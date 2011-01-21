function [ f_friction ] = friction_model( num_inputs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

g=num_inputs(12)*1e6; %leakage gap (microns)

%alpha = 3 (g - 0 to 30 microns)
%f_friction=0.000979*g^2 - 0.028555*g + 0.377921;     %friction factor for dry friction between teflon piston and aluminum cylinder

%alpha = 5 (g - 0 to 40 microns)
%f_friction = 0.000353*g^2 - 0.017133*g + 0.377921;

%alpha = 6 (g - 0 to 50 microns)
%f_friction = 0.000245*g^2 - 0.014277*g + 0.377921;

%alpha = 8 (g - 0 to 60 microns)
%f_friction = 0.000138*g^2 - 0.010708*g + 0.377921;

%alpha = 12 (g - 0 to 100 microns)
f_friction = 0.0000612*g^2 - 0.0071386*g + 0.3779207;
end

