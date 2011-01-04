function [ f_friction ] = friction_model( num_inputs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

g=num_inputs(12)*1e6; %leakage gap (microns)
f_friction=0.000979*g^2 - 0.028555*g + 0.377921;     %friction factor for dry friction between teflon piston and aluminum cylinder

end

