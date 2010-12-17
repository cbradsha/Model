function [p]=create_struct(num_inputs,txt_inputs)
%This function creates a structure of input parameters that will be used
%throughout the program. The title of this structure will be 'p' for
%parameter.

[T_i,rho_i,P_i,R,P_d,f,x_o,t_step,gamma,x_d,method,n,V_cv2_i,L_piston_max,g,D_piston,x_stroke,V_max,w_d,V_dead_valve,eta_motor,Ap,P_electric,x_max,leakage_on,mass_on,vibration_order,vibration_on,valve_dynamics_on,D_cv2,L_cv2,heat_transfer_on,h_2_s] = givens(num_inputs,txt_inputs);
[d_valve_suction,d_valve_discharge,C_d,A_suction_valve,A_discharge_valve,I_discharge,I_suction,k_discharge,k_suction,m_eff_discharge,m_eff_suction,A_suction,A_discharge,x_tr_suction,x_tr_discharge,x_stop,a_suction,a_discharge,d_suction,d_discharge] = valve_inputs(num_inputs,txt_inputs);
[ M_mov,J,L_1,L_2,k_mech,f_friction,x_piston_i,x_dot_piston_i,theta_i,theta_dot_i,x_piston_m_i,J_a,ecc_1,ecc_2,L_load_1,L_load_2 ] = vibration_givens(x_o,vibration_order,num_inputs,txt_inputs);
[alpha,k_r,A,B,C,T_w_i,R_shell,T_amb] = heat_transfer_givens(T_i,num_inputs,txt_inputs);

%%%%%%
%From General Givens
%%%%%%

p.T_i=T_i;
p.rho_i=rho_i;
p.P_i=P_i;
p.R=R;
p.P_d=P_d;
p.f=f;
p.x_o=x_o;
p.t_step=t_step;
p.gamma=gamma;
p.x_d=x_d;
p.method=method;
p.n=n;
p.V_cv2_i=V_cv2_i;
p.L_piston_max=L_piston_max;
p.g=g;
p.D_piston=D_piston;
p.x_stroke=x_stroke;
p.V_max=V_max;
p.w_d=w_d;
p.V_dead_valve=V_dead_valve;
p.eta_motor=eta_motor;
p.Ap=Ap;
p.P_electric=P_electric;
p.x_max=x_max;
p.leakage_on=leakage_on;
p.mass_on=mass_on;
p.vibration_order=vibration_order;
p.vibration_on=vibration_on;
p.valve_dynamics_on=valve_dynamics_on;
p.heat_transfer_on=heat_transfer_on;
p.D_cv2=D_cv2;
p.L_cv2=L_cv2;
p.h_2_s=h_2_s;


%%%%%%%%
%From valve givens
%%%%%%%%

p.d_valve_suction=d_valve_suction;
p.d_valve_discharge=d_valve_discharge;
p.C_d=C_d;
p.A_suction_valve=A_suction_valve;
p.A_discharge_valve=A_discharge_valve;
p.I_discharge=I_discharge;
p.I_suction=I_suction;
p.k_discharge=k_discharge;
p.k_suction=k_suction;
p.m_eff_discharge=m_eff_discharge;
p.m_eff_suction=m_eff_suction;
p.A_suction=A_suction;
p.A_discharge=A_discharge;
p.x_tr_suction=x_tr_suction;
p.x_tr_discharge=x_tr_discharge;
p.x_stop=x_stop;
p.a_suction=a_suction;
p.a_discharge=a_discharge;
p.d_suction=d_suction;
p.d_discharge=d_discharge;


%%%%%%%
%From Vibration Givens
%%%%%%%

p.M_mov=M_mov;
p.J=J;
p.L_1=L_1;
p.L_2=L_2;
p.k_mech=k_mech;
p.f_friction=f_friction;
p.x_piston_i=x_piston_i;
p.x_dot_piston_i=x_dot_piston_i;
p.theta_i=theta_i;
p.theta_dot_i=theta_dot_i;
p.x_piston_m_i=x_piston_m_i;
p.J_a=J_a;
p.ecc_1=ecc_1;
p.ecc_2=ecc_2;
p.L_load_1=L_load_1;
p.L_load_2=L_load_2;

%%%%%%%%
%Heat Transfer Inputs
%%%%%%%%

p.alpha=alpha;
p.k_r=k_r;
p.A=A;
p.B=B;
p.C=C;
p.T_w_i=T_w_i;
p.T_w_cv2_i=T_w_i;
p.R_shell=R_shell;
p.T_amb=T_amb;


%%%%%%
% Miscellaneous inputs
%%%%%%
p.n_period=num_inputs(54);
p.loop_error=num_inputs(55);
p.save_name=txt_inputs(56);
p.save_all=num_inputs(57);

%%%%%%%%%%%
% Frequency Sweep Parameters
%%%%%%%%%%%

p.f_begin=num_inputs(58);
p.f_increment=num_inputs(59);
p.f_end=num_inputs(60);



