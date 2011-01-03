function [ F_drive,Q_motor ] = motor( w_d,t,eta_motor,P_electric,t_force_adjust )
%A simple motor model for linear compressor

%global eta_motor P_electric

if isempty(t_force_adjust) == 1
    t_force_adjust=0;
end

%keyboard
P_motor=eta_motor*P_electric;
Q_motor=P_electric-P_motor;             %Heat transfer

%Converts power input to force output by motor, based on motor literature
k_m = -0.0008*P_motor^2 - 0.0047*P_motor + 4.6831;

F_drive_max=sqrt(P_motor)*(k_m);       %Maximum Motor force, in Newtons.

F_drive=F_drive_max*(cos(w_d*(t-t_force_adjust)+pi));                %Driving force v. Time, Newtons
   


end

