function [T_i,rho_i,P_i,R,P_d,f,x_o,t_step,gamma,x_d,method,n,V_cv2_i,L_piston_max,g,D_piston,x_stroke,V_max,w_d,V_dead_valve,eta_motor,Ap,P_electric,x_max,leakage_on,mass_on,vibration_order,vibration_on,valve_dynamics_on,D_cv2,L_cv2,heat_transfer_on,h_2_s,brent,T_s_2] = givens(num_inputs,txt_inputs)
%A function that intializes all given parameters

% input_data = xlsread('inputs.xls');
% txt_inputs{8} = input('Type of Solver Method, Euler or RK4:  ');
% txt_inputs{56} = input('Name for Study:  ');
% num_inputs = input_data(1+2,2:61);

%Compression Process ODE Solver
%method='RK4';  %Options are 'RK4', and 'Euler'
method=txt_inputs(8);

%Turn model options on and off
leakage_on=num_inputs(19);               %Determines if leakage model is on or not. 1 it is open, 0 is no leakage
mass_on=num_inputs(20);                  %Determines if valves open or not, 0 valves stay closed, 1 they open normally
vibration_order=num_inputs(21);          %Order of Vibration model, 1=1st order, 2=2nd order, all else is first order
vibration_on=num_inputs(22);             %Determines if vibration model is used or not (0 is off, 1 is on)  If not the piston is assumed to move in a sine wave.
valve_dynamics_on=num_inputs(23);        %Determines if valve dynamics are on or off (1 is on, 0 is a digital valve)
heat_transfer_on=1;                      %Heat transfer in cylinder on (1) or off (2)
brent = 1;                               %Turn Brent's method on or off

%Inlet Properties
T_i=num_inputs(1);
T_i=T_i+273.15;                 %K
%P_i=572;                    %kPa
P_i=num_inputs(2);
rho_i=EOS(T_i,P_i,'rho','P');   %kg/m^3
%R=0.08149;                  %kJ/kg-K
R=num_inputs(3);

%Discharge Pressure
%P_d=715;                    %in kPa
P_d=num_inputs(4);

%Input Parameters
%f=44;           %Frequency of operation, Hz
f=num_inputs(5);
%n=1000;
n=num_inputs(9);
w_d=f*2*pi;                 %Damped natural frequency, assumed to be the driving frequency to drive resonance.

t_step=1/(f*n);       %time step in seconds

%gamma=1.239;
gamma=num_inputs(6);

%Piston Parameters
%D_piston=0.49*(0.0254);         % Piston Diameter in meters
D_piston=num_inputs(13);
%x_d=0.4*(0.0254); %free volume height in meters (dead space)
x_d=num_inputs(7);
%x_max=1*(0.0254);        %Max Stroke in meters
x_max=num_inputs(18);
%x_stroke=0.1*(0.0254);         %Inital Stroke Guess in meters
x_stroke=num_inputs(14);
Ap=(pi*D_piston^2)/4;                   %Piston Area in m^2
V_max=Ap*x_max;                  %Maximum Cylinder volume, m^3
x_o=(x_max)/2;                  %Piston Starting Position
%x_o=0.002523;
%V_dead_valve=V_max/25;          %estimate of dead volume in valves (FIX, from PRO-E)
V_dead_valve=num_inputs(15);


% M_mov=312/1000;  %moving mass in kg
% L_1=1.75*(0.0254);  %Length from cylinder head to COM in meters
% L_2=1.61*(0.0254);  %Length from rear bearing to COM in meters
% epsilon=0.1*(0.0254); %Eccentricy moment arm, for SIDE LOADS
% k_2=2134*100; %Torsional Spring rate in N/m, obtained from FEA simulation

%Control Volume 2 Givens
%V_cv2_i=1e-4;      %Nominal Volume of CV2 (WRONG, FIND FROM PRO-E)
V_cv2_i=num_inputs(10);
D_cv2=2*(0.0254);   %Internal Diameter of CV2
L_cv2=10*(0.0254);  %Length of CV2

%L_piston_max=2.5*(0.0254);      % Max length of piston for leakage path, in meters
L_piston_max=num_inputs(11);
%g=0.0005*(0.0254);              % Gap between piston and cylinder in meters
g=num_inputs(12);


%Motor Givens
%eta_motor=0.5;              %Motor efficiency
eta_motor=num_inputs(16);
%P_electric=20;              %Motor input power, in Watts
P_electric=num_inputs(17);

%Isentropic relations calculations, using REFROP
% s_1=refpropm('S','T',T_i,'D',rho_i,'R134a');
% h_2_s=refpropm('H','P',P_d,'S',s_1,'R134a');
% h_2_s=h_2_s-refpropm('H','T',(-40+273.15),'P',P_d,'R134a');
% h_2_s=h_2_s/1000;

%Isentropic Relations, using EOS
%s_ref = EOS(233.15,1417.64,'entropy','rho');
s_1=EOS(T_i,rho_i,'entropy','rho')


%Secant Method to Find Exit Temp for Isentropic Compression
dT = 1;     %Initalize change in Temperature
T(1) = 320; %Guess Values
T(2) = 325;
i=2;

%Secant Loop
while abs(dT)>1e-6
    T(i+1)=T(i)-(EOS(T(i),P_d,'entropy','P')-s_1)*(T(i)-T(i-1))/(EOS(T(i),P_d,'entropy','P')-EOS(T(i-1),P_d,'entropy','P'))
    
    dT=T(i+1)-T(i);
    i=i+1;
    if i>500
        error('stuck in interations')
    end
end

%Last Temperature in Iteration is Temperature at Isentropic Input
T_s_2 = T(i);
%Using Temperature and Discharge Pressure, enthalpy is Calcualted
h_2_s = EOS(T_s_2,P_d,'enthalpy','P');
