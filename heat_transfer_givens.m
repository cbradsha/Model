
function [alpha,k_r,A,B,C,T_w_i,R_shell,T_amb] = heat_transfer_givens(T_i,num_inputs,txt_inputs)
%Inputs to the heat transfer model

alpha=0.00000055;         %Thermal diffusivity in m^2
k_alum=160;               %Thermal conductivity of aluminum, W/m-K
A=0.25;
B=0.65;
C=0.25;                   %Coefficients from correlation
T_w_i=22.2+273.15;                %Inital wall temperature, in K
T_amb=22+273.15;                %Ambient Temperature

k_r=15*10^-3;              %Thermal conductivity of refrigerant, W/m-K

rho_air=1.2;                %Density of air, kg/m^3
Vel_air=1;                  %Velocity of airflow across compressor, m/s
D_shell=2*0.0254;           %Outer Diameter of Piston Cylinder, m
L_shell=10*0.0254;          %Length of Compressor, m
A_shell=pi*D_shell*L_shell;       %Surface area of compressor shell, m^2
k_air=30e-3;                 %W/m-K , conductivity of air
mu_air=2*10^-5;             %Dynamic Viscosity of air, kg/m-s

nu_air=mu_air/rho_air;      %Kinematic Viscosity of air, m^2/s

alpha_air=2.2160*10^-5;      %Thermal diffusivity of air, m^2/s

Pr=nu_air/alpha_air;        %Prandtl Number of air flow
Re_air=(rho_air*Vel_air*D_shell)/mu_air;    %Reynolds number of air

C_2=0.683;                    %Constants from Table 7.2 in Incropera DeWitt
m=0.466;

Nu_D=C_2*Re_air^m*Pr^(1/3);

h_shell=(Nu_D*k_air)/D_shell;      %Heat transfer coefficient for outside heat transfer

R_shell=1/(h_shell*A_shell);        %Overall thermal resistance from shell to air