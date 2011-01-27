function [ dm_leak_in,dm_leak_out,Ma_cv2 ] = leakage( P,P_cv2,x_dot_piston,T,rho,T_cv2,x_piston,p )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Pr_limit=0.556999;          %Sonic limit to the pressure ratio, for R-134a

if p.leakage_on==1

    %Couette-Posille Flow, for when Ma<0.3

    mu=5.87509602E-07+3.79308232E-08*T;  %approximation for mu between 290 and 390K, in kg/m-s

    L_piston=p.L_piston_max-x_piston;

    dP_dx=(P_cv2-P)/L_piston;
    dP_dx = dP_dx*1000;  %convert to N/m^3
     
    u_bar=(-x_dot_piston/2)+(p.g^2/(4*mu))*(-dP_dx)+(p.g^3/(6*mu))*dP_dx;

    A=(pi/4)*(p.D_piston+2*p.g)^2-(pi/4)*(p.D_piston^2);
    
    
    Ma_1=abs(u_bar)/sqrt(p.gamma*p.R*T);

    Ma_2=abs(u_bar)/sqrt(p.gamma*p.R*T_cv2);

    if Ma_1>Ma_2
        Ma_cv2=Ma_1;
    else
        Ma_cv2=Ma_2;
    end
    
    % Isentropic compressible flow when Ma>=0.3
    
    if Ma_cv2 >= 0.3    %Compressible flow
        
        Pr=min(P/P_cv2,P_cv2/P);  %grab the PR smaller than one
        
        
        if P >= P_cv2  %flow into CV2
            Ma_cv2=sqrt((2/(p.gamma-1))*((Pr^((1-p.gamma)/p.gamma))-1));
            if Ma_cv2>1
                Ma_cv2=1;
            end
            u_bar=sqrt(p.gamma*p.R*T*1000)*Ma_cv2;   %T because this flow comes from compression chamber
                
                
        elseif P < P_cv2  %flow into compression chamber
            Ma_cv2=sqrt((2/(p.gamma-1))*((Pr^((1-p.gamma)/p.gamma))-1));
            if Ma_cv2>1
                Ma_cv2=1;
            end
            u_bar=-sqrt(p.gamma*p.R*T_cv2*1000)*Ma_cv2;   %T because this flow comes from CV2   
            
        else
            error('your logic sucks in the leakage sub-model')
        end
        
    end

    dm_leak=rho*u_bar*A;

    if dm_leak<0   %flow into compression chamber
        dm_leak_out=0;
        dm_leak_in=abs(dm_leak);
    elseif dm_leak>0   %flow out of compression chamber
        dm_leak_out=abs(dm_leak);
        dm_leak_in=0;
    else
        dm_leak_in=0;
        dm_leak_out=0;
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

