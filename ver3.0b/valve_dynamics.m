function [ x_valve,x_dot_valve,dm,Ma ] = valve_dynamics(P,rho,T,x_valve_1,x_dot_valve_1,p)
%Valve flow model

Pr_limit=0.556999;          %Sonic limit to the pressure ratio, for R-134a

if p.mass_on==1
    if p.valve_dynamics_on==1

        %keyboard

        %If there is flowrate coming in...
        if p.P_i>P      
            Pr=P/p.P_i;
            
            if Pr > Pr_limit
                Ma=sqrt((2/(p.gamma-1))*((Pr^((1-p.gamma)/p.gamma))-1));
                if Ma>1
                    Ma=1;
                end
                V=sqrt(p.gamma*p.R*p.T_i*1000)*Ma;   %T_i because this flow comes from suction valve
            else
                %error('Sonic Limit Reached at Suction Valve')
                Ma=1;
                V=sqrt(p.gamma*p.R*p.T_i*1000)*Ma;   %T_i because this flow comes from suction valve
            end
            
            x_1=x_valve_1;
            x_2=x_dot_valve_1;

            if x_valve_1>=p.x_tr_suction


                k1 = p.t_step*x_RK_flux_dom(x_1,x_2,p.rho_i,p.A_suction_valve,p.k_suction,V,p.m_eff_suction,p.C_d,p.A_suction);

                k2 = p.t_step*x_RK_flux_dom(x_1+0.5*k1(1),x_2+0.5*k1(2),p.rho_i,p.A_suction_valve,p.k_suction,V,p.m_eff_suction,p.C_d,p.A_suction);
                k3 = p.t_step*x_RK_flux_dom(x_1+0.5*k2(1),x_2+0.5*k2(2),p.rho_i,p.A_suction_valve,p.k_suction,V,p.m_eff_suction,p.C_d,p.A_suction);
                k4 = p.t_step*x_RK_flux_dom(x_1+k3(1),x_2+k3(2),p.rho_i,p.A_suction_valve,p.k_suction,V,p.m_eff_suction,p.C_d,p.A_suction);

                x_total=[x_1;x_2]+(1/6)*(k1+2*k2+2*k3+k4);


            else 
                %keyboard
                k1 = p.t_step*x_RK_pressure_dom(x_1,x_2,p.rho_i,p.A_suction_valve,p.k_suction,V,p.m_eff_suction,p.C_d,p.P_i,P);
                k2 = p.t_step*x_RK_pressure_dom(x_1+0.5*k1(1),x_2+0.5*k1(2),p.rho_i,p.A_suction_valve,p.k_suction,V,p.m_eff_suction,p.C_d,p.P_i,P);
                k3 = p.t_step*x_RK_pressure_dom(x_1+0.5*k1(1),x_2+0.5*k1(2),p.rho_i,p.A_suction_valve,p.k_suction,V,p.m_eff_suction,p.C_d,p.P_i,P);
                k4 = p.t_step*x_RK_pressure_dom(x_1+k1(1),x_2+k1(2),p.rho_i,p.A_suction_valve,p.k_suction,V,p.m_eff_suction,p.C_d,p.P_i,P);


                x_total=[x_1;x_2]+(1/6)*(k1+2*k2+2*k3+k4);


            end

            if x_total(1)<0
                x_total(1)=0;              %This is the hard stop limit, valve cannot have a negative position
            end   
    %         if x_total(1)>1.4*x_stop
    %             x_total(1)=1.4*x_stop;      %Hard stop limit for discharge valve
    %         end

            if x_total(1)==0 %|| x_total(1)==1.4*x_stop
                x_total(2)=0;
            end
            %keyboard
            x_valve=x_total(1);         
            x_dot_valve=x_total(2);
            
            x_right=0.5*p.d_suction+p.a_suction;
            x_left=p.a_suction-0.5*p.d_suction;
            
            x_ave_suction=(1/(x_right-x_left))*(((3*x_valve)/p.a_suction^2)*(p.a_suction^3-x_left^3)-((3*x_valve)/(4*p.a_suction^3))*(p.a_suction^4-x_left^4)...               
                +((9*x_valve)/(2*p.a_suction))*(x_right^2-p.a_suction^2)-3*x_valve*(x_right-p.a_suction));
            
            A_port=p.A_suction;
            
            if pi*p.d_suction*x_ave_suction < A_port
                dm=-p.rho_i*V*pi*p.d_suction*x_ave_suction;
            else
                dm=-p.rho_i*V*A_port;
            end
            %keyboard    

        %if there is flowrate coming out...
        elseif P>=p.P_d
           Pr=p.P_d/P;
           if Pr > Pr_limit
               Ma=sqrt((2/(p.gamma-1))*((Pr^((1-p.gamma)/p.gamma))-1));
                if Ma>1
                    Ma=1;
                end
               V=sqrt(p.gamma*p.R*T*1000)*Ma;
           else
               %error('Sonic Limit Reached at Discharge Port')
               Ma=1;
               V=sqrt(p.gamma*p.R*T*1000)*Ma;
           end
           x_1=x_valve_1;
           x_2=x_dot_valve_1;

            if x_valve_1>=p.x_tr_discharge


                k1 = p.t_step*x_RK_flux_dom(x_1,x_2,rho,p.A_discharge_valve,p.k_discharge,V,p.m_eff_discharge,p.C_d,p.A_discharge);

                k2 = p.t_step*x_RK_flux_dom(x_1+0.5*k1(1),x_2+0.5*k1(2),rho,p.A_discharge_valve,p.k_discharge,V,p.m_eff_discharge,p.C_d,p.A_discharge);
                k3 = p.t_step*x_RK_flux_dom(x_1+0.5*k2(1),x_2+0.5*k2(2),rho,p.A_discharge_valve,p.k_discharge,V,p.m_eff_discharge,p.C_d,p.A_discharge);
                k4 = p.t_step*x_RK_flux_dom(x_1+k3(1),x_2+k3(2),rho,p.A_discharge_valve,p.k_discharge,V,p.m_eff_discharge,p.C_d,p.A_discharge);

                x_total=[x_1;x_2]+(1/6)*(k1+2*k2+2*k3+k4);

            else 

                k1 = p.t_step*x_RK_pressure_dom(x_1,x_2,rho,p.A_discharge_valve,p.k_discharge,V,p.m_eff_discharge,p.C_d,P,p.P_d);
                %keyboard
                k2 = p.t_step*x_RK_pressure_dom(x_1+0.5*k1(1),x_2+0.5*k1(2),rho,p.A_discharge_valve,p.k_discharge,V,p.m_eff_discharge,p.C_d,P,p.P_d);
                k3 = p.t_step*x_RK_pressure_dom(x_1+0.5*k1(1),x_2+0.5*k1(2),rho,p.A_discharge_valve,p.k_discharge,V,p.m_eff_discharge,p.C_d,P,p.P_d);
                k4 = p.t_step*x_RK_pressure_dom(x_1+k1(1),x_2+k1(2),rho,p.A_discharge_valve,p.k_discharge,V,p.m_eff_discharge,p.C_d,P,p.P_d);


                x_total=[x_1;x_2]+(1/6)*(k1+2*k2+2*k3+k4);

            end

            if x_total(1)<0
                x_total(1)=0;              %This is the hard stop limit, valve cannot have a negative position
            end
            if x_total(1)>p.x_stop
                x_total(1)=p.x_stop;      %Hard stop limit for discharge valve
            end

            if x_total(1)==0 || x_total(1)==p.x_stop
                x_total(2)=0;
            end

            %keyboard
            x_valve=x_total(1);         
            x_dot_valve=x_total(2);   
            
            x_right=0.5*p.d_discharge+p.a_discharge;
            x_left=p.a_discharge-0.5*p.d_discharge;
            
            x_ave_discharge=(1/(x_right-x_left))*(((3*x_valve)/p.a_discharge^2)*(p.a_discharge^3-x_left^3)-((3*x_valve)/(4*p.a_discharge^3))*(p.a_discharge^4-x_left^4)...               
                +((9*x_valve)/(2*p.a_discharge))*(x_right^2-p.a_discharge^2)-3*x_valve*(x_right-p.a_discharge));
            
            A_port=(pi*p.d_discharge^2)/4;
            
            if pi*p.d_valve_discharge*x_ave_discharge < A_port
                dm=rho*V*pi*p.d_valve_discharge*x_ave_discharge;  
            else
                dm=rho*V*A_port;
            end

        else
           dm=0;
           Ma=0;
           x_valve=0;
           x_dot_valve=0;


        end
        
        
    elseif p.valve_dynamics_on==0
        
        
        %If there is flowrate coming in...
        if p.P_i>P      
            Pr=P/p.P_i;
            Ma=sqrt((2/(p.gamma-1))*((Pr^((1-p.gamma)/p.gamma))-1));
            V=sqrt(p.gamma*p.R*(p.T_i)*1000)*Ma;
            dm=-p.rho_i*V*p.A_suction;  %negative means flow in!!
            x_valve=0;
            x_dot_valve=0;

        %dm=0;


        %if there is flowrate coming out...
        elseif P>=p.P_d
           Pr=p.P_d/P;
           Ma=sqrt((2/(p.gamma-1))*((Pr^((1-p.gamma)/p.gamma))-1));
           V=sqrt(p.gamma*p.R*(T)*1000)*Ma;
           dm=rho*V*p.A_discharge;  %positive means flow out!!
           x_valve=0;
           x_dot_valve=0;

        %dm=0;

        else
           dm=0;
           Ma=0;
           x_valve=0;
           x_dot_valve=0;

        end
        
    end

elseif p.mass_on==0
    
    x_valve=0;
    dm=0;
    x_dot_valve=0;
    Ma=0;
    
else

    x_valve=0;
    dm=0;
    x_dot_valve=0;
    Ma=0;
end

    

end

