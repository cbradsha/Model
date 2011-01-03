function [ dT,drho ] = DT_RK_cv2( t,T_cv2,P,P_cv2,x_dot_piston_1,T,rho,rho_cv2,h,T_w_cv2,x_piston_1,theta_1,theta_dot_1,x_piston_m_1,dP,dx,theta_dot_dot_1,p)
%the function for RK analysis
%   Detailed explanation goes here


   %global T_i rho_i P_i R P_d f x_o t_step gamma r_p x_d method n V_cv2_i L_piston_max g D_piston x_stroke V_max w_d V_dead_valve eta_motor Ap P_electric
   %global d_valve_suction d_valve_discharge C_d A_suction_valve A_discharge_valve I_discharge I_suction k_discharge k_suction m_eff_discharge m_eff_suction A_suction A_discharge x_tr_suction x_tr_discharge x_stop
   %global M_mov J ecc L_1 L_2 k_2 k_mech f_friction x_piston_i x_dot_piston_i theta_i theta_dot_i x_piston_m_i
%     global V_cv2_i
    
    %Volume from geometry
    %[dV,V,x] = volume(x_o,t,f,r_p,x_d);
    [ dV,V,x_piston,x_dot_piston,theta,theta_dot,x_piston_m,theta_dot_dot ] = vibration( t,dP,dx,x_piston_1,x_dot_piston_1,theta_1,theta_dot_1,x_piston_m_1,theta_dot_dot_1,p );
    dV_cv2=dV;
    V_cv2=p.V_cv2_i-V;

    %Leakage Massflow
    [dm_leak_in dm_leak_out,Ma_cv2] = leakage( P,P_cv2,x_dot_piston,T,rho,T_cv2,x_piston,p );
    
    %Total flowrates - in and out are relative to compression chamber
    dm=dm_leak_out-dm_leak_in;
    dm_in=dm_leak_out;
    dm_out=dm_leak_in;
    
    [ F_drive,Q_motor ] = motor( p.w_d,t,p.eta_motor,p.P_electric,p.t_force_adjust );
    
    %Heat Transfer
    [ Q_cv2 ] = Ins_HT_cv2( T_cv2,rho_cv2,T_w_cv2,V_cv2,dV_cv2,x_dot_piston,Q_motor,p);
    
    %Mass of the gas
    m=rho_cv2*V_cv2;    

%Fluid Properties for approximation
    h_in_cv2=h;
    h_cv2=EOS(T_cv2,rho_cv2,'enthalpy','rho');
    dP_dT_v= dP_dT(T_cv2,rho_cv2);
    Cv_cv2=EOS(T_cv2,rho_cv2,'Cv','rho');
    
%%Function
    dT=(Q_cv2+T_cv2*dP_dT_v*(1/1000)*((1/rho_cv2)*dm-dV_cv2)-h*dm+dm_in*h_in_cv2-dm_out*h_cv2)/(m*Cv_cv2);
    
    drho=(dm-rho_cv2*dV_cv2)/V_cv2;
    
    if length(dT)>1
        keyboard
    end

end

