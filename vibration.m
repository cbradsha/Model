function [ dV,V,x_piston,x_dot_piston,theta,theta_dot,x_piston_m,theta_dot_dot,p ] = vibration( t,dP,dx,x_piston_1,x_dot_piston_1,theta_1,theta_dot_1,x_piston_m_1,theta_dot_dot_1,p )
%Vibration Model Component of Linear Compressor Model, added ver. 1.4.

% 
% global P_i P_d x_o t_step gamma x_stroke V_max V_dead_valve eta_motor Ap P_electric w_d x_d
% %global d_valve_suction d_valve_discharge C_d A_suction_valve A_discharge_valve I_discharge I_suction k_discharge k_suction m_eff_discharge m_eff_suction A_suction A_discharge x_tr_suction x_tr_discharge x_stop
% global M_mov J L_1 L_2 k_mech f_friction vibration_on J_a ecc_1 ecc_2 L_load_1 L_load_2

if p.vibration_on==1

    V=x_piston_1*p.Ap+p.V_dead_valve;           %Piston volume, when x_piston=0 volume is dead volume in valves.

    %dV=-x_dot_piston_1*p.Ap;
    dV=-x_dot_piston_1*p.Ap;        %x_piston and x_piston_m are opposites, x_dot_piston follows x_piston_m but I want it to follow x_piston, hence the negative.
        
    %k_gas=abs((dP*p.Ap*1000)/dx);           %Calculated Spring Rate of the Gas in the cylinder (N/m), 1000 converts from 
%     if p.dP_max>=(p.P_d-p.P_i)
%     
%         k_gas=((p.P_d-p.P_i)*p.Ap*1000)/p.x_stroke;  %Gas spring rate, linear estimate
%         
%     elseif p.dP_max<(p.P_d-p.P_i)
%         
%         k_gas=(p.dP_max*p.Ap*1000)/p.x_stroke;  %Gas spring rate, linear estimate
%         
%     end

    p.theta=theta_1;
    
    %Linearized vibration model flag
    linear = 0;
    
    if linear == 1
    
        theta_tmp = atan(p.g/p.L_1);
        p.F_wall=(1/p.L_1)*(p.k_mech*(p.x_stroke_ref-p.ecc_1*theta_tmp)*p.ecc_1);
        %k_gas adjusted for max stroke
        k_gas=((p.P_d-p.P_i)*p.Ap*1000)/p.x_stroke_ref;  %Gas spring rate, linear estimate
        p.V_max=p.Ap*p.x_stroke_ref;
        W_gas=((p.gamma*p.P_i*p.V_max)/(p.gamma-1))*((p.P_d/p.P_i)^((p.gamma-1)/p.gamma)-1);        %Work done on gas in one cycle, kJ
        c_friction=(4*p.f_friction*p.F_wall)/(p.w_d*p.x_stroke_ref*pi);     %effective damping due to dry friction in cylinder   
        c_gas=(W_gas*1000)/(p.w_d*p.x_stroke_ref^2*pi);                %effecitve damping due to work done on gas, 1000 to convert work to J, (N-s)/m

    else
        
        if abs(theta_1) >= atan(p.g/p.L_1)  %if piston is in contact with wall, there is friction.
            
            p.F_wall=(1/p.L_1)*(p.k_mech*abs((abs(x_piston_m_1)-p.ecc_1*abs(theta_1)))*p.ecc_1); %Newtons

        else
            p.F_wall=0;
            %keyboard
            p.F_wall=(1/p.L_1)*(p.k_mech*abs((abs(x_piston_m_1)-p.ecc_1*abs(theta_1)))*p.ecc_1); %Newtons
        end
        
        %p.V_max=p.Ap*(x_piston_m_1+p.x_o);
        p.V_max = p.Ap*p.x_stroke;
        if p.dP_max>=(p.P_d-p.P_i)
        
            k_gas=((p.P_d-p.P_i)*p.Ap*1000)/p.x_stroke;  %Gas spring rate, linear estimate
            %W_gas=((p.gamma*p.P_i*p.V_max)/(p.gamma-1))*((p.P_d/p.P_i)^((p.gamma-1)/p.gamma)-1);        %Work done on gas in one cycle, kJ
            W_gas = (-p.P_current*p.dV + 0.5*p.dV*dP)*1000;
            
        elseif p.dP_max<(p.P_d-p.P_i)
            
            k_gas=(p.dP_max*p.Ap*1000)/p.x_stroke;  %Gas spring rate, linear estimate
            %W_gas=((p.gamma*p.P_i*p.V_max)/(p.gamma-1))*(((p.dP_max/p.P_i)+1)^((p.gamma-1)/p.gamma)-1);        %Work done on gas in one cycle, kJ
            W_gas = (-p.P_current*p.dV + 0.5*p.dV*dP)*1000;
        end
        
        %Testing F_gas instead of k_gas
        k_gas = 0;
        
        
        c_friction=(4*p.f_friction*p.F_wall)/(p.w_d*p.x_stroke*pi);
        %%effective damping due to dry friction in cylinder
        c_gas=(W_gas)/(p.w_d*p.x_stroke^2*pi);
        %%effecitve damping due to work done on gas, 1000 to convert work to J, (N-s)/m
    end
    
    
    k_eff=k_gas+p.k_mech;         %Effective Spring Rate, gas + mechanical springs
    p.k_eff=k_eff;
    
    c_eff=c_gas+c_friction;                         %Total effective damping, (N-s)/m
    p.c_eff = c_eff;
    p.c_friction = c_friction;
    p.c_gas = c_gas;
    %keyboard

%     k1=p.t_step* dx_vibration( t,x_piston_m_1,x_dot_piston_1,theta_1,theta_dot_1,c_eff,p.M_mov,p.J,p.P_electric,p.w_d,p.eta_motor,p.J_a,p.ecc_1,p.ecc_2,p.L_load_1,p.L_load_2,p.k_mech,k_gas,p.L_1,p.L_2,p.t_force_adjust );
%     k2=p.t_step* dx_vibration( t+0.5*p.t_step,x_piston_m_1+0.5*k1(1),x_dot_piston_1+0.5*k1(2),theta_1+0.5*k1(3),theta_dot_1+0.5*k1(4),c_eff,p.M_mov,p.J,p.P_electric,p.w_d,p.eta_motor,p.J_a,p.ecc_1,p.ecc_2,p.L_load_1,p.L_load_2,p.k_mech,k_gas,p.L_1,p.L_2,p.t_force_adjust );
%     k3=p.t_step* dx_vibration( t+0.5*p.t_step,x_piston_m_1+0.5*k2(1),x_dot_piston_1+0.5*k2(2),theta_1+0.5*k2(3),theta_dot_1+0.5*k2(4),c_eff,p.M_mov,p.J,p.P_electric,p.w_d,p.eta_motor,p.J_a,p.ecc_1,p.ecc_2,p.L_load_1,p.L_load_2,p.k_mech,k_gas,p.L_1,p.L_2,p.t_force_adjust );
%     k4=p.t_step* dx_vibration( t+p.t_step,x_piston_m_1+k3(1),x_dot_piston_1+k3(2),theta_1+k3(3),theta_dot_1+k3(4),c_eff,p.M_mov,p.J,p.P_electric,p.w_d,p.eta_motor,p.J_a,p.ecc_1,p.ecc_2,p.L_load_1,p.L_load_2,p.k_mech,k_gas,p.L_1,p.L_2,p.t_force_adjust );
% 
%     dx_total=[x_piston_m_1,x_dot_piston_1,theta_1,theta_dot_1]+(1/6)*(k1+2*k2+2*k3+k4);

    %%%%%%%%%%%%%%%%%%
    %% Start - dx_vibration replacement
    %%%%%%%%%%%%%%%%%%

    %Vector of Timesteps
    t_combined = [t,t+0.5*p.t_step,t+0.5*p.t_step,t+p.t_step];
    %k_m = -0.0008*(p.eta_motor*p.P_electric)^2 - 0.0047*(p.eta_motor*p.P_electric) + 4.6831;
    k_m = 4.6831;
    F_drive_max=sqrt(p.eta_motor*p.P_electric)*(k_m);       %Maximum Motor force, in Newtons.
    F_drive=F_drive_max*(cos(p.w_d*(t-p.t_force_adjust)+pi));                %Driving force v. Time
    
    if linear == 0
        F_drive = F_drive - p.dP_piston*p.Ap*1000;
    end
    
    %Renaming displacements and angles for Runge-Kutta Calculations
    x_1 = x_piston_m_1;
    x_2 = x_dot_piston_1;
    theta_1_calc = theta_1;
    theta_2 = theta_dot_1;

    %Replacing dx_vibration%
    %[ F_drive,Q_motor ] = motor( p.w_d,t_combined(1),p.eta_motor,p.P_electric,p.t_force_adjust );
    F_drive=F_drive_max*(cos(p.w_d*(t_combined(1)-p.t_force_adjust)+pi));                %Driving force v. Time
    
    if linear == 0
        F_drive = F_drive - p.dP_piston*p.Ap*1000;
    end
    
    num2=(1/p.M_mov)*(p.k_mech*p.ecc_1*theta_1_calc+F_drive-c_eff*x_2-(k_gas+p.k_mech)*x_1);
    num4=(p.k_mech/p.J)*(x_1-p.ecc_1*theta_1)*p.ecc_1;
    dx=[x_2,num2,theta_2,num4];
    k1 = p.t_step*dx;
    
    
    %Step 2, update displacements and angles
    x_1 = x_1+0.5*k1(1);
    x_2 = x_2+0.5*k1(2);
    theta_1_calc = theta_1_calc+0.5*k1(3);
    theta_2 = theta_2+0.5*k1(4);
    
    %Replacing dx_vibration%
    %[ F_drive,Q_motor ] = motor( p.w_d,t_combined(2),p.eta_motor,p.P_electric,p.t_force_adjust );
    F_drive=F_drive_max*(cos(p.w_d*(t_combined(2)-p.t_force_adjust)+pi));

    if linear == 0
        F_drive = F_drive - p.dP_piston*p.Ap*1000;
    end
    
    num2=(1/p.M_mov)*(p.k_mech*p.ecc_1*theta_1_calc+F_drive-c_eff*x_2-(k_gas+p.k_mech)*x_1);
    num4=(p.k_mech/p.J)*(x_1-p.ecc_1*theta_1)*p.ecc_1;
    dx=[x_2,num2,theta_2,num4];
    k2 = p.t_step*dx;
    
    %Step 3, update displacements and angles
    x_1 = x_1+0.5*k2(1);
    x_2 = x_2+0.5*k2(2);
    theta_1_calc = theta_1_calc+0.5*k2(3);
    theta_2 = theta_2+0.5*k2(4);
    
    %Replacing dx_vibration%
    num2=(1/p.M_mov)*(p.k_mech*p.ecc_1*theta_1_calc+F_drive-c_eff*x_2-(k_gas+p.k_mech)*x_1);
    num4=(p.k_mech/p.J)*(x_1-p.ecc_1*theta_1)*p.ecc_1;
    dx=[x_2,num2,theta_2,num4];
    k3 = p.t_step*dx;
    
    %Step 4, update displacements and angles
    x_1 = x_1+k2(1);
    x_2 = x_2+k2(2);
    theta_1_calc = theta_1_calc+k2(3);
    theta_2 = theta_2+k2(4);
    
    %Replacing dx_vibration%
    %[ F_drive,Q_motor ] = motor( p.w_d,t_combined(4),p.eta_motor,p.P_electric,p.t_force_adjust );
    F_drive=F_drive_max*(cos(p.w_d*(t_combined(4)-p.t_force_adjust)+pi));
    
    if linear == 0
        F_drive = F_drive - p.dP_piston*p.Ap*1000;
    end
    
    num2=(1/p.M_mov)*(p.k_mech*p.ecc_1*theta_1_calc+F_drive-c_eff*x_2-(k_gas+p.k_mech)*x_1);
    num4=(p.k_mech/p.J)*(x_1-p.ecc_1*theta_1)*p.ecc_1;
    dx=[x_2,num2,theta_2,num4];
    k4 = p.t_step*dx;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %% End, dx_vibration replacement
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    dx_total=[x_piston_m_1,x_dot_piston_1,theta_1,theta_dot_1]+(1/6)*(k1+2*k2+2*k3+k4);
    
    x_piston_m=dx_total(1);
    x_dot_piston=dx_total(2);
    theta=dx_total(3);
    theta_dot=dx_total(4);
    theta_dot_dot=k1(4);
    p.x_dot_dot_piston=k1(2);

%     if x_piston_m <= 0
%         x_piston=x_piston_m+p.x_o;
%     else 
%         x_piston=p.x_o-x_piston_m;
%     end

    x_piston=-x_piston_m+p.x_o;

    if x_piston < 0                 %This second condition does not allow the piston to travel beyond the top of the cylinder
        x_piston = 0;
        x_piston_m=p.x_o;
        x_dot_piston=0;
    end
    
    if abs(theta) > atan(p.g/p.L_1)          %Restricts the rotation to approx 5 deg.
        if theta < 0
            theta=-atan(p.g/p.L_1);
            theta_dot=0;
        elseif theta > 0
            theta=atan(p.g/p.L_1);
            theta_dot=0;
        end
    end
    
    %damping_ratio=c_eff/(2*k_eff*p.M_mov);
    %p.w_d=p.w_d*sqrt(1-damping_ratio^2);
    
elseif p.vibration_on==0
    
    %V_o=0;
    %Ap=((pi*r^2)/4); %in square meters

    %w_d=2*pi*f;

    %x_piston=p.x_o*(cos(p.w_d*t+pi)+1)+p.x_d;
    x_piston=p.x_o*(cos(p.w_d*t+pi)+1);

    V=p.Ap*x_piston+p.V_dead_valve;

    dV=-p.Ap*p.x_o*p.w_d*sin(p.w_d*t+pi);

    x_dot_piston=-p.x_o*p.w_d*(sin(p.w_d*t+pi));
    theta=0;
    theta_dot=0;
    %x_piston_m=0;
    theta_dot_dot=0;
    x_piston_m=x_piston-p.x_o;
    p.F_wall=0;
    p.theta=0;
    p.k_eff=0;
    p.c_eff = 0;
    p.c_friction = 0;
    p.c_gas = 0;
    p.x_dot_dot_piston=0;
end


%keyboard

end

