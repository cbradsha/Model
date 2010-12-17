function [ J ] = jacobian( T_1,rho_1,T_1_cv2,rho_1_cv2,x_stroke_save,dT_res,drho_res,dT_res_cv2,drho_res_cv2,dx_res,k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Temp for CV1 components

den(1)=((T_1(k)-T_1(k-1))*(1/T_1(k)));
den(2)=((rho_1(k)-rho_1(k-1))*(1/rho_1(k)));
den(3)=((T_1_cv2(k)-T_1_cv2(k-1)*(1/T_1_cv2(k))));
den(4)=(rho_1_cv2(k)-rho_1_cv2(k-1))*(1/rho_1_cv2(k));
den(5)=(x_stroke_save(k)-x_stroke_save(k-1))*(1/x_stroke_save(k));

% num(1)=(dT_res(k)-dT_res(k-1));
% num(2)=(drho_res(k)-drho_res(k-1));
% num(3)=(dT_res_cv2(k)-dT_res_cv2(k-1));
% num(4)=(drho_res_cv2(k)-drho_res_cv2(k-1));
% num(5)=(dx_res(k)-dx_res(k-1));

num(1)=(dT_res(k));
num(2)=(drho_res(k));
num(3)=(dT_res_cv2(k));
num(4)=(drho_res_cv2(k));
num(5)=(dx_res(k));

J=zeros(5);

for i=1:length(num)
    for j=1:length(den)
        J(i,j)=num(i)/den(j);
    end
end

% J(1,1)=(dT_res(k)-dT_res(k-1))/den_1;
% J(1,2)=(dT_res(k)-dT_res(k-1))/den_2;
% J(1,3)=(dT_res(k)-dT_res(k-1))/den_3;
% J(1,4)=(dT_res(k)-dT_res(k-1))/den_4;
% J(1,5)=(dT_res(k)-dT_res(k-1))/den_5;
% 
% 
% J(2,1)=(drho_res(k)-drho_res(k-1))/(T_1(k)-T_1(k-1));
% J(2,2)=(drho_res(k)-drho_res(k-1))/(rho_1(k)-rho_1(k-1));
% J(2,3)=(drho_res(k)-drho_res(k-1))/(T_1_cv2(k)-T_1_cv2(k-1));
% J(2,4)=(drho_res(k)-drho_res(k-1))/(rho_1_cv2(k)-rho_1_cv2(k-1));
% J(2,5)=(drho_res(k)-drho_res(k-1))/(x_stroke_save(k)-rho_1_cv2(k-1));
% 
% 
% J(3,1)=(dT_res_cv2(k)-dT_res_cv2(k-1))/(T_1(k)-T_1(k-1));
% J(3,2)=(dT_res_cv2(k)-dT_res_cv2(k-1))/(rho_1(k)-rho_1(k-1));
% J(3,3)=(dT_res_cv2(k)-dT_res_cv2(k-1))/(T_1_cv2(k)-T_1_cv2(k-1));
% J(3,4)=(dT_res_cv2(k)-dT_res_cv2(k-1))/(rho_1_cv2(k)-rho_1_cv2(k-1));
% J(3,5)=(dT_res_cv2(k)-dT_res_cv2(k-1))/(x_stroke_save(k)-rho_1_cv2(k-1));
% 
% 
% J(4,1)=(drho_res_cv2(k)-drho_res_cv2(k-1))/(T_1(k)-T_1(k-1));
% J(4,2)=(drho_res_cv2(k)-drho_res_cv2(k-1))/(rho_1(k)-rho_1(k-1));
% J(4,3)=(drho_res_cv2(k)-drho_res_cv2(k-1))/(T_1_cv2(k)-T_1_cv2(k-1));
% J(4,4)=(drho_res_cv2(k)-drho_res_cv2(k-1))/(rho_1_cv2(k)-rho_1_cv2(k-1));
% J(4,5)=(drho_res_cv2(k)-drho_res_cv2(k-1))/(x_stroke_save(k)-rho_1_cv2(k-1));
% 
% 
% J(5,1)=(dx_res(k)-dx_res(k-1))/(T_1(k)-T_1(k-1));
% J(5,2)=(dx_res(k)-dx_res(k-1))/(rho_1(k)-rho_1(k-1));
% J(5,3)=(dx_res(k)-dx_res(k-1))/(T_1_cv2(k)-T_1_cv2(k-1));
% J(5,4)=(dx_res(k)-dx_res(k-1))/(rho_1_cv2(k)-rho_1_cv2(k-1));
% J(5,5)=(dx_res(k)-dx_res(k-1))/(x_stroke_save(k)-rho_1_cv2(k-1));


end

