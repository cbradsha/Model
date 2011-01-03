function [ dT,drho ] = DT_RK( t,T,rho,P,T_w,x_valve_1,x_dot_valve_1,x_piston_1,x_dot_piston_1,theta_1,theta_dot_1,x_piston_m_1,dP,dx,P_cv2,T_cv2,theta_dot_dot_1,p)
%the function for RK analysis
%   Detailed explanation goes here
    
    %Volume from Vibration Model
    %[dV,V,x] = volume(x_o,t,f,r_p,x_d);
    [ dV,V,x_piston,x_dot_piston,theta,theta_dot,x_piston_m,theta_dot_dot ] = vibration( t,dP,dx,x_piston_1,x_dot_piston_1,theta_1,theta_dot_1,x_piston_m_1,theta_dot_dot_1,p );

    %Valve Massflow
    %[dm_valve,Ma]= valve(P_i,P_d,P,rho_i,rho,T_i,T,R,gamma);
    [ x_valve,x_dot_valve,dm_valve, Ma ] = valve_dynamics(P,rho,T,x_valve_1,x_dot_valve_1,p);
    
    %Leakage Flowrates (Leakage Model)
    [dm_leak_in dm_leak_out,Ma_cv2] = leakage( P,P_cv2,x_dot_piston_1,T,rho,T_cv2,x_piston_1,p );
    

    %Total flowrates 
    if dm_valve<0     %flow IN
        dm_in=abs(dm_valve)+dm_leak_in;
        dm_out=dm_leak_out;
        dm=dm_in-dm_out;

    elseif dm_valve>0 %flow OUT
        dm_in=dm_leak_in;
        dm_out=dm_valve+dm_leak_out;
        dm=dm_in-dm_out;
    else
        dm_in=dm_leak_in;
        dm=dm_leak_in-dm_leak_out;
        dm_out=dm_leak_out;
    end
    
        %Heat Transfer
   [ Q ]=Ins_HT( T,rho,T_w,V,dV,x_piston_1,x_dot_piston_1,p);
         
    
    %Mass of the gas
    m=rho*V;    

%Fluid Properties for approximation
    %h_in=EOS(T_i,rho_i,'enthalpy','rho');
    h=EOS(T,rho,'enthalpy','rho');
    dP_dT_v= dP_dT(T,rho);
    Cv=EOS(T,rho,'Cv','rho');
    
%%Function
    dT=(Q+T*dP_dT_v*(1/1000)*((1/rho)*dm-dV)-h*dm+dm_in*p.h_in-dm_out*h)/(m*Cv);
    
    drho=(dm-rho*dV)/V;
    
    %keyboard

end

