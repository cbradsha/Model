function [d_valve_suction,d_valve_discharge,C_d,A_suction_valve,A_discharge_valve,I_discharge,I_suction,k_discharge,k_suction,m_eff_discharge,m_eff_suction,A_suction,A_discharge,x_tr_suction,x_tr_discharge,x_stop,a_suction,a_discharge,d_suction,d_discharge] = valve_inputs(num_inputs,txt_inputs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Valve Flapper Geometry Information
d_valve_suction= num_inputs(24);  %Suction valve flapper diameter in meters
d_valve_discharge=num_inputs(25);  %Discharge valve flapper diameter in meters
E=num_inputs(26);  % (Pa) Young's Modulus of Spring Steel for valves, from Marks Std. Handbook
h_valve=num_inputs(27);
%h_valve=0.008*(0.0254); %height of valve (thickness) in meters
rho_valve=num_inputs(28); % kg/m^3, density of reed valve material
C_d=num_inputs(29);  %Drag coefficient for the valve
%C_d=90;

l_suction_valve=num_inputs(30);  %length of suction valve, meters
l_discharge_valve=num_inputs(31);  %length of discharge valve, meters

A_suction_valve=num_inputs(32);  % Suction Valve area,m^2 (From Pro-E Model)
A_discharge_valve=num_inputs(33);  %Discharge Valve area,m^2 (From Pro-E Model)

a_discharge=num_inputs(34);  % Distance from anchor to force,meters (from PRO_E)
a_suction=num_inputs(35);

%Valve flapper strength of materials values
I_discharge=(d_valve_suction*h_valve^3)/12;  %Moment of Intertia for discharge valve,m^4
I_suction=(d_valve_discharge*h_valve^3)/12;  % Moment of Intertia for suction valve,m^4

k_discharge=(6*E*I_discharge)/(a_discharge^2*(3*l_discharge_valve-a_discharge));
k_suction=(6*E*I_suction)/(a_suction^2*(3*l_suction_valve-a_suction));
%k_suction=599.41;       %N/m from Pro-E
%k_discharge=5216.48;        %N/m from Pro-E

m_eff_discharge=(1/3)*rho_valve*A_discharge_valve*h_valve;
m_eff_suction=(1/3)*rho_valve*A_suction_valve*h_valve;

w_n_discharge=sqrt(k_discharge/m_eff_discharge);
w_n_suction=sqrt(k_suction/m_eff_suction);

%Suction Valve Port Geometry Information
d_suction=num_inputs(36);  %suction oval diameter, meters
A_suction=((pi*((d_suction^2)/4))+0.12*0.0254*d_suction);  %Suction Port Area, m^2
d_h_suction=(4*A_suction)/(pi*d_suction+2*(0.12*0.0254));


%Discharge Valve Port Geometry Information
d_discharge=num_inputs(37);  %discharge diameter in meters
A_discharge=(pi*((d_discharge^2)/4));  %Discharge Port area,m^2

%Transitional valve lift
x_tr_suction=0.25*(d_h_suction^2/d_valve_suction);  % Transitional Valve Lift,meters
x_tr_discharge=0.25*(d_discharge^2/d_valve_discharge);
%Discharge Valve stopper distance
x_stop=num_inputs(38);      %Height of discharge stopper in meters
      
%keyboard

%end

