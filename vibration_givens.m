function [ M_mov,J,L_1,L_2,k_mech,f_friction,x_piston_i,x_dot_piston_i,theta_i,theta_dot_i,x_piston_m_i,J_a,ecc_1,ecc_2,L_load_1,L_load_2 ] = vibration_givens(x_o,vibration_order,num_inputs,txt_inputs)
%Inputs to vibration model

% global x_o vibration_order

L=num_inputs(40);        %Total length of piston
L_1=num_inputs(41);                %Distance from CG to front of piston
L_2=L-L_1;                %Distance from CG to back of piston
R=num_inputs(42);        %Radius of Piston Spring Seats
M_mov=num_inputs(39);          %Piston Equivalent moving mass, in kg
%M_mov=0.28;
M_piston=num_inputs(43);          %Mass of just the piston

J=((1/12)*0.8*M_piston*L^2)+((1/12)*0.2*M_piston*(0.98*0.0254)^2);      %Moment of Inertia from CG (kg-m^2)
J_a=((1/3)*0.8*M_piston*L^2)+((0.25)*0.2*M_piston*R^2)+(0.2*M_piston*L_2^2);    %Moment of Inertia from back of Piston (kg-m^2)

if vibration_order==2
    %ecc=0.02*(0.0254);       %Eccentricity from spring forces, in meters
    ecc_1=num_inputs(50);
    ecc_2=num_inputs(51);
elseif vibration_order==1
    ecc_1=0;
    ecc_2=0;
else
    ecc_1=0;
    ecc_2=0;
end

L_load_1=num_inputs(52);             %Distance between spring seats, 1 is forward spring, 2 is toward motor.
L_load_2=num_inputs(53);

%k_mech=8647.5;          %Mechanical Stiffness of Springs (N/m)
k_mech=num_inputs(44);
%k_mech=8000;

model_exist = exist('friction_model','file');
if model_exist == 2
    %new friction model
    f_friction=friction_model(num_inputs);
else
    f_friction=num_inputs(45);
end

% % Initial Conditions for vibration model
x_piston_i=x_o;
x_dot_piston_i=num_inputs(46);
theta_i=num_inputs(47);
theta_dot_i=num_inputs(48);
x_piston_m_i=num_inputs(49);

end

