function [ Q_cv2 ] = Ins_HT_cv2( T_cv2,rho_cv2,T_w_cv2,V_cv2,dV_cv2,x_dot_piston,Q_motor,p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if p.heat_transfer_on==1 && p.leakage_on == 1

%%%%%%%%
%Approximate Viscosity
%%%%%%%%

mu_cv2=5.87509602E-07+3.79308232E-08*T_cv2;  %approximation for mu between 290 and 390K, in kg/m-s

%%%%%%%
%Fagotti Correlation, 1998 Purdue Compressor conference

%%%%%%%
%Compressibility Number
%%%%%%%


L_cv2=((p.gamma-1)/V_cv2)*dV_cv2*sqrt(p.D_cv2^3/(p.alpha*abs(p.x_dot_ave)));

if isnan(L_cv2)
    L_cv2=0;
end

%%%%%%%
%Reynolds Numbers for the process
%%%%%%%


Re_cv2=(rho_cv2*abs(p.x_dot_ave)*p.D_cv2)/mu_cv2;

%%%%%%%
%Heat Transfer Area
%%%%%%%


A_ht_cv2=p.D_cv2*pi*p.L_cv2;

%%%%%%
%Calculate Heat Transfer
%%%%%%


Q_cv2=((A_ht_cv2*p.k_r*(T_cv2-T_w_cv2))/p.D_cv2)*(p.A*(Re_cv2^p.B)+p.C*L_cv2*(T_w_cv2/(T_cv2-T_w_cv2)));

if isnan(Q_cv2)
    Q_cv2=0;
end

Q_cv2=((-Q_cv2/1000)+((-Q_motor)/1000));

else
    Q_cv2=0;
end

end

