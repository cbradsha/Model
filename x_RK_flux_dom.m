function [a] = x_RK_flux_dom(x_1,x_2,rho,A_valve,k_valve,V,m_eff,C_d,A_port)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%C_d=1.17;  %Add to valve_dynamics

a=[x_2;(rho*A_port*(V-x_2)^2-k_valve*x_1+0.5*C_d*rho*V^2*A_valve-0.5*C_d*rho*x_2^2*A_valve)/m_eff];

%x_dot_1=a(1);
%x_dot_2=a(2);

%keyboard

end

