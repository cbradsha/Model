clc
clear all
close all


f=42;
p.w_d=2*pi*f;
t_step=1/(100*f);
Period=1/f;
%t=0:t_step:Period;
i=1;
p.n_period=1;
t_cross=0;
p.P_electric=20;
p.eta_motor=0.5;
t_force_adjust=0;
t=0;

for k=1:1:20
    
    t_cross=0;
    t=0;
    F_drive=0;
    x_piston=0;
    


    while max(size(t_cross))<p.n_period*2+1


        [ F_drive(i) ] = motor( p.w_d,t(i),p.eta_motor,p.P_electric,t_force_adjust );

        x_piston(i)=10*(cos(p.w_d*(1-(1/(50*k)))*(t(i))+(pi))); 
        
        if i>1
            t_cross=find_cross(t',x_piston',0);
        end
        
        t(i+1)=t(i)+t_step;
        i=i+1;
        
        %keyboard

    end
    t=t(1:end-1);
    %disp(k)
    i=1;
    temp(k)=max(t_cross)-max(find_cross(t',F_drive',0));
    t_force_adjust=t_force_adjust+temp(k);
    t_change(k)=t_force_adjust;
    
    if k==1
        plot(t,F_drive);hold on; plot(t,x_piston,'k');
    elseif k==2
        plot(t,F_drive,'r');plot(t,x_piston,'ko');
    elseif k==3
        plot(t,F_drive,'c');plot(t,x_piston,'k+-');
    elseif k==4
        plot(t,F_drive,'g');plot(t,x_piston,'k+');
    elseif k==5
        plot(t,F_drive,'k--');plot(t,x_piston,'ko-');
    end
    hold on

end

figure(2)
plot(t_change)            
%plotyy(t,F_drive,t,x_piston)
