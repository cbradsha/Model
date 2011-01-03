function [a] = x_RK_pressure_dom(x_1,x_2,rho,A_valve,k_valve,V,m_eff,C_d,P_high,P_low)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%C_d=1.17;  %Add to valve_dynamics

%keyboard
a=[x_2;(-k_valve*x_1+0.5*C_d*rho*V^2*A_valve+A_valve*(P_high-P_low)*1000-0.5*C_d*rho*x_2^2*A_valve)/m_eff];

%x_dot_1=a(1);
%x_dot_2=a(2);

%keyboard

end

