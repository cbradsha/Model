function [ dx ] = dx_vibration( t,x_1,x_2,theta_1,theta_2,c_eff,M_mov,J,P_electric,w_d,eta_motor,J_a,ecc_1,ecc_2,L_load_1,L_load_2,k_mech,k_gas,L_1,L_2,t_force_adjust )
%this function calculates the differentials for the vibration model


[ F_drive,Q_motor ] = motor( w_d,t,eta_motor,P_electric,t_force_adjust );

num2=(1/M_mov)*(k_mech*ecc_1*theta_1+F_drive-c_eff*x_2-(k_gas+k_mech)*x_1);
num4=(k_mech/J)*(x_1-ecc_1*theta_1)*ecc_1;

dx=[x_2,num2,theta_2,num4];


end

