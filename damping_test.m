
%A= [1 0 0 0;0 p.M_mov 0 0;0 0 1 0;0 0 0 p.J];
c_eff_iter = linspace(1,15);
M = [p.M_mov 0;0 p.J];
%K = [k_eff 0;0 p.k_mech*p.ecc_1^2];
K = [k_eff -p.k_mech*p.ecc_1;-p.k_mech*p.ecc_1 p.k_mech*p.ecc_1^2];


for j = 1:1:length(c_eff_iter)
    
    %B = [0 -1 0 0;k_eff c_eff_iter(j) -p.k_mech*p.ecc_1 0;0 0 0 -1;-p.k_mech*p.ecc_1 0 p.k_mech*p.ecc_1^2 0];
    %C = inv(A)*B;
    C = [c_eff 0;0 0];
    A = [zeros(2) eye(2);-(M\K) -(M\C)];
    [v,d]= eig(A);
    w_d(j) = imag(d(1,1));
    f_d(j) = w_d(j)/(2*pi);
        
end