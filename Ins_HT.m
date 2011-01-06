function [ Q ] = Ins_HT( T,rho,T_w,V,dV,x_piston,x_dot_piston,p)
%Instanteneous heat transfer model for heat transfer between the cylinder
%wall and the compressor gas.  Used in all control volumes.

if p.heat_transfer_on==1

%%%%%%%%
%Approximate Viscosity
%%%%%%%%
mu=5.87509602E-07+3.79308232E-08*T;  %approximation for mu between 290 and 390K, in kg/m-s



%%%%%%%
%Fagotti Correlation, 1998 Purdue Compressor conference

%%%%%%%
%Compressibility Number
%%%%%%%

L=((p.gamma-1)/V)*dV*sqrt(p.D_piston^3/(p.alpha*abs(p.x_dot_ave)));

if isnan(L)
    L=0;
end

%%%%%%%
%Reynolds Numbers for the process
%%%%%%%

Re=(rho*abs(p.x_dot_ave)*p.D_piston)/mu;


%%%%%%%
%Heat Transfer Area
%%%%%%%

A_ht=p.D_piston*pi*x_piston;


%%%%%%
%Calculate Heat Transfer
%%%%%%

Q=((A_ht*p.k_r*(T-T_w))/p.D_piston)*(p.A*(Re^p.B)+p.C*L*(T_w/(T-T_w)));

W_dot_friction = (p.c_eff - p.c_gas)*p.x_max*p.f_list^2;

if isnan(Q)
    Q=0;
end

Q=-Q/1000+W_dot_friction/1000;


else
    
    Q=0;
end



end

