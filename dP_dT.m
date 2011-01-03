function [ dP_dT ] = dP_dT( T,rho )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%T=400;
%rho=1/0.04;

P=EOS(T,rho,'pressure','rho');

dT=0.0001;

dP_dT=((EOS(T+dT,rho,'pressure','rho')-P)/dT)*1000;

end

