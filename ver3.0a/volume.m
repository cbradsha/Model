function [dV,V,x,x_dot_piston] = volume(x_o,t,f,r,x_d)
%Volume is a function that calculate the volume of the piston/cylinder
%as a function of time.


%V_o=0;
A_piston=((pi*r^2)/4); %in square meters

w=2*pi*f;

x=x_o*(cos(w*t+pi)+1)+x_d;

V=A_piston*x;

dV=-A_piston*x_o*w*sin(w*t+pi);

x_dot_piston=-x_o*w*(sin(w*t+pi));


end

