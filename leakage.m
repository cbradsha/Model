function [ dm_leak_in,dm_leak_out,Ma_cv2 ] = leakage( P,P_cv2,x_dot_piston,T,rho,T_cv2,x_piston,p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% global gamma L_piston_max g D_piston leakage_on R
%global d_valve_suction d_valve_discharge C_d A_suction_valve A_discharge_valve I_discharge I_suction k_discharge k_suction m_eff_discharge m_eff_suction A_suction A_discharge x_tr_suction x_tr_discharge x_stop
%global M_mov J ecc L_1 L_2 k_2 k_mech f_friction x_piston_i x_dot_piston_i theta_i theta_dot_i x_piston_m_i

if p.leakage_on==1

    %Couette-Posille Flow, for when Ma<0.3

    mu=5.87509602E-07+3.79308232E-08*T;  %approximation for mu between 290 and 390K, in kg/m-s

    L_piston=p.L_piston_max-x_piston;

    dP_dx=(P_cv2-P)/L_piston;
     
    u_bar=(-x_dot_piston/2)+(p.g^2/(4*mu))*(-dP_dx)+(p.g^3/(6*mu))*dP_dx;

    A=(pi/4)*(p.D_piston+2*p.g)^2-(pi/4)*(p.D_piston^2);

    dm_leak=rho*u_bar*A;

    if dm_leak<0
        dm_leak_out=0;
        dm_leak_in=abs(dm_leak);
    elseif dm_leak>0
        dm_leak_out=abs(dm_leak);
        dm_leak_in=0;
    else
        dm_leak_in=0;
        dm_leak_out=0;
    end

    Ma_1=abs(u_bar)/sqrt(p.gamma*p.R*T);

    Ma_2=abs(u_bar)/sqrt(p.gamma*p.R*T_cv2);

    if Ma_1>Ma_2
        Ma_cv2=Ma_1;
    else
        Ma_cv2=Ma_2;
    end

elseif p.leakage_on==0
    
    Ma_cv2=0;
    dm_leak_in=0;
    dm_leak_out=0;
    
else
    
    Ma_cv2=0;
    dm_leak_in=0;
    dm_leak_out=0;
end

end

