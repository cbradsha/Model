
%Main Model Program
%Craig Bradshaw
%June 2009
%Herrick Labs, cbradsha@purdue.edu
%
%Version History
%
%Version 1.0
% - Inital Compression process model
%
%Version 1.1
% - Improved speed by eliminating EOS calls
% - Fixed unit errors in DT_RK sub-function
% - Fixed unit errors in Valve model
% - Added diagnostic plots
% - 
%
%Verson 1.2
% - Added Valve Dynamics
%
%Version 1.2.1
% - Added option of utilizing Property information from REFPROP.
% Information comes from REFPROP v.7 through a link in MATLAB.
% REFPROP slightly faster than EOS coded in
% MATLAB.
%
% - Added option to utilize ideal gas properties in EOS of model.  Operator
% changed in EOS.m instead of in givens.m function
%
% - Changed time step initialization.  The user now supplies the number of
% steps and the frequency of operation and the time step is calculated
% based on this information.
%
%Version 1.2.2
% - Found several bugs in the valve model, fixed.
% - Added a damping term to both the flux dominated and pressure dominated
% regions of valve operation.
% - Added total flowrate change to residuals list.
%
%Version 1.3
% - Major Revision to Add Leakage model
% - Found a sign error in the DT-RK subfunction that generated impossible
% results.  Fixed in all previous model versions.
%
%Version 1.3.1
% - Made some revisions to the valve model to make it more accurate,
% including changing some of the geometric values so they were closer to
% reality.
%
%Version 1.4
% - Added  2 degree of freedom Vibration Model
% - Added options to turn on and off each part of the model for diagnostics
% - Fixed a bug in the valve model, an incorrect area was being used in the
% flux dominant equation
% - Updated general model to incorporate global variables 
%
%Version 1.5
% - Added heat tranfer model
% - Changed constant input parameters from global variables to a single
% struct that is passed to each function in an effort to inprove speed.
% - Added the ability to batch process many different parametric runs using
% the excel spreadsheet 'inputs.xls'.
%
%Version 2.0
% - Attempted to reduce convergence time by adding Newtons Method, failed.
% - Changed integration techniquest for heat transfer and mass flow from 
%   rectangular method to trapezoidal method.
% - Added ability to calculate compressor losses (valve, motor, heat
%   transfer, friction, gas work) 

%%%%%%%%%%%%
%TO - DO
%%%%%%%%%%%%

% 1 - Add viscosity to EOS
% 2 - Add diffusivity to EOS
% 3 - Add dP_dT_v to EOS
% 4 - Finish Adding Losses
% 5 - Add compressibility to leakage model


try

clc
clear all
close all



%% Initializations

pc_flag = input('Is this run on a server or PC without Excel COM Server? (1 - yes, 0 - no)');


if pc_flag == 0
    %Note to shut down Excel
    disp('**************************')
    disp('Please shut down all Excel Notebooks')
    disp('**************************')
    disp('   ')
    disp('Once notebooks are closed, press any key to continue')
    pause

    %How long do I need to run this simulation for? Obtained from first two
    %numbers on the inputs.xls file.
    num_size=xlsread('inputs.xls','A3:A4');
    analy_length=(num_size(2)-num_size(1))+1;


    %Open Excel Activex Server to use xlsread/write for prelim stuff
    Excel = actxserver ('Excel.Application');
    File=strcat(pwd,'\data_log.xls');
    if ~exist(File,'file')
        ExcelWorkbook = Excel.workbooks.Add;
        ExcelWorkbook.SaveAs(File,1);
        ExcelWorkbook.Close(false);
    end
    invoke(Excel.Workbooks,'Open',File);


    %Now using modified xlsread/write to work on data_log.xls
    [num,txt,raw]=xlsread1('data_log.xls');
    [rows,col]=size(raw);


    %Adding a note to the batch log so I know what I was doing
    batch_on=input('Is this a batch file? 1 for yes, 0 for no ');
    if batch_on==1
        batch_note=input('Please input a string to describe why you are running this batch.');
        batch_note=cellstr(batch_note);
        batch_range=strcat('A',num2str(rows+2),':A',num2str(rows+2));
        xlswrite1('data_log.xls',batch_note,batch_range);
    end

    %Shutting down activex server to excel
    invoke(Excel.ActiveWorkbook,'Save');
    Excel.Quit
    Excel.delete
    clear Excel
    
elseif pc_flag == 1
    
    input_data = xlsread('inputs.xls');
    
    %How long do I need to run this simulation for? Obtained from first two
    %numbers on the inputs.xls file. 
    analy_length = (input_data(4,1) - input_data(3,1))+1; 
    
        
else
    error('Invalid pc_flag.  Please use either 0 or 1.');
    clear all
    
end

%turn annoying warning off
warning('off','MATLAB:Print:SavingToDifferentName')

%%%%%%%%%%%%%
% For loop used for batch operation, line by line from input file
%%%%%%%%%%%%%

for z=1:1:analy_length
    
    if pc_flag == 0
    
        %pulling input data from inputs.xls
        range=strcat('B',num2str(z+2),':','BI',num2str(z+2));
        [num_inputs,txt_inputs]=xlsread('inputs.xls',range);

        %creating my variable structure
        [p]=create_struct(num_inputs,txt_inputs);

        %name that I use to generate a save folder
        p.save_name=strcat(p.save_name,'_',num2str(z));


        %Open Excel Activex Server to use xlsread/write for the remainder of
        %the analysis
        Excel = actxserver ('Excel.Application');
        File=strcat(pwd,'\data_log.xls');
        if ~exist(File,'file')
            ExcelWorkbook = Excel.workbooks.Add;
            ExcelWorkbook.SaveAs(File,1);
            ExcelWorkbook.Close(false);
        end
        invoke(Excel.Workbooks,'Open',File);


        %Open data log
        [num,txt,raw]=xlsread1('data_log.xls');
        [rows,col]=size(raw);
        
    elseif pc_flag == 1
        
        if z == 1
            txt_inputs{8} = 'Euler';
            txt_inputs{56} = input('Name for Study:  ');
        end
        
        num_inputs = input_data(z+2,2:61);
        
        %Create structure of constant variables
        [p]=create_struct(num_inputs,txt_inputs);

        %name that I use to generate a save folder
        p.save_name=strcat(p.save_name,'_',num2str(z));
        
        
    end
    

    %Time/Volume Calculations
    error=0.0001;
    f_list=p.f_begin:p.f_increment:p.f_end;
    w_d_list=f_list*2*pi;
    p.x_stroke_initial = p.x_stroke;

%%%%%%%%%
% Loop used for frequency sweeps
%%%%%%%%%

p.resonant_loop = 0;        %initialize variable that tells frequency sweep loop if we have reached the maximum frequency.
l=1;                        %loop variable for the following loop.
    
while p.resonant_loop<1 && l<=length(w_d_list)
%for l=1:1:length(w_d_list)

    close all
    p.w_d=w_d_list(l);
    p.Period=1/f_list(l);     %Time of one cycle
    p.t_step=1/(f_list(l)*p.n);       %time step in seconds
    %t=0:p.t_step:p.n_period*p.Period; 
    p.transient=0;
    t_adjust_force=0;
    
    
    k=1;
    dT_loop_k=1;
    drho_loop_k=1;
    dT_cv2_loop_k=1;
    drho_cv2_loop_k=1;
    dx_piston=1;

    %% Initalize Vectors
    %These are the vectors that are iterated on and are thus not included in
    %the parameters structure

    P=zeros(1,p.n_period*p.n);
    dT=zeros(1,p.n_period*p.n);
    h=zeros(1,p.n_period*p.n-1);
    V=zeros(1,p.n_period*p.n-1);
    dV=zeros(1,p.n_period*p.n-1);
    dm=zeros(1,p.n_period*p.n-1);
    dm_in=zeros(1,p.n_period*p.n-1);
    dm_out=zeros(1,p.n_period*p.n-1);
    drho=zeros(1,p.n_period*p.n-1);
    m=zeros(1,p.n_period*p.n-1);
    rho=zeros(1,p.n_period*p.n);
    x_piston=zeros(1,p.n_period*p.n-1);
    x_dot_piston=zeros(1,p.n_period*p.n-1);
    Cv=zeros(1,p.n_period*p.n-1);
    iter=zeros(1,p.n_period*p.n);
    dP_dT_v=zeros(1,p.n_period*p.n-1);
    Ma=zeros(1,p.n_period*p.n-1);
    dT_check=zeros(1,p.n_period*p.n-1);
    T=zeros(1,p.n_period*p.n);
    x_valve=zeros(1,p.n_period*p.n);
    x_dot_valve=zeros(1,p.n_period*p.n);
    dm_leak_in=zeros(1,p.n_period*p.n-1);
    dm_leak_out=zeros(1,p.n_period*p.n-1);
    dm_out_valve=zeros(1,p.n_period*p.n-1);
    dV_cv2=zeros(1,p.n_period*p.n-1);
    V_cv2=zeros(1,p.n_period*p.n-1);
    Ma_cv2=zeros(1,p.n_period*p.n-1);
    h_cv2=zeros(1,p.n_period*p.n-1);
    Cv_cv2=zeros(1,p.n_period*p.n-1);
    dP_dT_v_cv2=zeros(1,p.n_period*p.n);
    P_cv2=zeros(1,p.n_period*p.n);
    drho_cv2=zeros(1,p.n_period*p.n-1);
    rho_cv2=zeros(1,p.n_period*p.n-1);
    dT_cv2=zeros(1,p.n_period*p.n-1);
    T_cv2=zeros(1,p.n_period*p.n-1);
    F_drive=zeros(1,p.n_period*p.n);
    Q_motor=zeros(1,p.n_period*p.n);
    theta=zeros(1,p.n_period*p.n);
    theta_dot=zeros(1,p.n_period*p.n);
    dP=zeros(1,p.n_period*p.n);
    dx=zeros(1,p.n_period*p.n);
    theta_dot_dot=zeros(1,p.n_period*p.n);
    Q=zeros(1,p.n_period*p.n);
    Q_cv2=zeros(1,p.n_period*p.n);
    m_cv2=zeros(1,p.n_period*p.n);
    F_wall=zeros(1,p.n_period*p.n);
    theta_1_c=zeros(1,p.n_period*p.n);
    k_eff=zeros(1,p.n_period*p.n);
    %F_drive = zeros(1,p.n_period*p.n);
   
     x_piston_m=0;
     
    T_w = zeros(1,p.loop_error);
    T_w_cv2 = zeros(1,p.loop_error);
    x_stroke_save = zeros(1,p.loop_error);
    loop = zeros(1,p.loop_error);
    dT_res = zeros(1,p.loop_error);
    drho_res = zeros(1,p.loop_error);
    dx_res = zeros(1,p.loop_error);
    dT_res_cv2 = zeros(1,p.loop_error);
    drho_res_cv2 = zeros(1,p.loop_error);
    dT_loop = zeros(1,p.loop_error);
    drho_loop = zeros(1,p.loop_error);
    dT_cv2_loop = zeros(1,p.loop_error);
    drho_cv2_loop = zeros(1,p.loop_error);    

    mu=1;       %placeholder mu spot, should add viscosity to EOS.


    m_dot=0;
    m_dot_leak_out=0;
    m_dot_leak_in=0;
    W_dot=0;
    eta_o=0;
    eta_vol=0;
    y=1;        %loop variable for piston 'bumps'
    loop=1;     %loop counter variable, re-initialize
    p.dP_max=p.P_d-p.P_i;
     

    
    
    
 %keyboard


         %%%%%%%%%%%%
         %  while loop to enforce convergence
         %%%%%%%%%%%%
     
     
        while (abs(dT_loop_k)>error || abs(drho_loop_k)>error) || (abs(dT_cv2_loop_k)>error || abs(drho_cv2_loop_k)>error) || abs(dx_piston)>error
            

            

            %% Inital Cylinder Conditions
            if k==1

                i=1;
                T(i)=p.T_i;
                rho(i)=p.rho_i;

                P(i)=EOS(T(i),rho(i),'pressure','rho');

                T_cv2(1)=p.T_i;
                rho_cv2(1)=p.rho_i;

                P_cv2(1)=P(1);

                x_piston(1)=p.x_piston_i;
                x_dot_piston(1)=p.x_dot_piston_i;
                V(1)=(p.x_piston_i)*p.Ap+p.V_dead_valve;
                dV(1)=-p.x_dot_piston_i*p.Ap;
                theta(1)=p.theta_i;
                theta_dot(1)=p.theta_dot_i;
                x_piston_m(1)=p.x_piston_m_i;       %%%%%%IMPORTANT, AVERAGE SPEED
                p.x_dot_ave=1;

                dP(1)=0;
                dx(1)=1;
                
                T_w(k)=p.T_w_i;
                T_w_cv2(k)=p.T_w_i;
                
                x_stroke_save(k)=p.x_stroke_initial;
                p.x_stroke = p.x_stroke_initial;
                p.dP_max = p.P_d-p.P_i;
                
                t=0;

            else
                                        %This forces the new inlet conditions to
                                        %equal the previous outlet conditions
                %%%%%%%%%%%%%%%%%%%%%
                %Update without Newtons Method
                %%%%%%%%%%%%%%%%%%%%%
                
%                 T(1)=T(i);
%                 rho(1)=rho(i);
%                 P(1)=P(i);
% 
%                 T_cv2(1)=T_cv2(i);
%                 rho_cv2(1)=rho_cv2(i);
%                 P_cv2(1)=P_cv2(i);
%                p.x_stroke=max(x_piston_m)-min(x_piston_m); %max stroke

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Update values, using Newtons Method
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                T(1)=X_n_1(1);
                rho(1)=X_n_1(2);
                %P(1)=pressure(T(1),rho(1));
                P(1)=EOS( T(1),rho(1),'pressure','rho' );

                T_cv2(1)=X_n_1(3);
                rho_cv2(1)=X_n_1(4);
                %P_cv2(1)=pressure(T_cv2(1),rho_cv2(1));
                P_cv2(1)=EOS( T_cv2(1),rho_cv2(1),'pressure','rho' );

                %p.x_stroke=X_n_1(5); %max stroke
                p.x_stroke=max(x_piston_m)-min(x_piston_m); %max stroke
                p.x_dot_ave=0.707*max(x_dot_piston);            %%%%%%IMPORTANT, AVERAGE SPEED
                
                
%                 if k==20*y                      %check every 20 iterations
%                     p.x_stroke=1.2*p.x_stroke;   %20% higher value, bump.
%                     p.x_dot_ave=1.2*p.x_dot_ave;
%                     y=y+1;
%                 end

                
                
                %disp(strcat('x_dot_ave: ',num2str(p.x_dot_ave)))
                x_stroke_save(k)=p.x_stroke;
                %disp(strcat('x_stroke: ',num2str(x_stroke_save(k))))
    
                
                x_piston(1)=x_piston(i);
                %x_piston(1)=p.x_piston_i;
                x_dot_piston(1)=x_dot_piston(i);

                V(1)=V(i);

                dV(1)=dV(i);
                theta(1)=theta(i);
                theta_dot(1)=theta_dot(i);
                x_piston_m(1)=x_piston_m(i);
                t=t(i);

                dP(1)=dP(i);
                dx(1)=dx(i);
                %disp(i)
                i=1;                %resets inner loop variable
                %disp(i)
                
                
            end
            
            %% Initalize inner loop variables
            
            t_cross=0;
            p.t_force_adjust=0;
            %t=0;
            x_piston_m=0;
            F_drive=0;
            %F_drive = zeros(1,p.n_period*p.n);
            %x_piston_m = zeros(1,p.n_period*p.n);

            %%%%%%%%%%
            %Main Loop for calculations
            %%%%%%%%%%
  
            %for i=1:1:length(t)
            while size(t_cross)<p.n_period*2-1
                
       

                %% Motor Model
                [ F_drive(i),Q_motor(i) ] = motor( p.w_d,t(i),p.eta_motor,p.P_electric,p.t_force_adjust );            %Motor Model

                %% Vibration Model and Volume and Piston Position  

                [ dV(i+1),V(i+1),x_piston(i+1),x_dot_piston(i+1),theta(i+1),theta_dot(i+1),x_piston_m(i+1),theta_dot_dot(i+1),p ] = vibration( t(i),dP(i),dx(i),x_piston(i),x_dot_piston(i),theta(i),theta_dot(i),x_piston_m(i),theta_dot_dot(i),p );

                F_wall(i)=p.F_wall;
                theta_1_c(i)=p.theta;
                k_eff(i)=p.k_eff;
                
                %% Volume for Second Control Volume
                dV_cv2(i)=dV(i);
                V_cv2(i)=p.V_cv2_i-V(i);

                %% Massflow (Valve Model)
                [ x_valve(i+1),x_dot_valve(i+1),dm_valve, Ma(i) ] = valve_dynamics(P(i),rho(i),T(i),x_valve(i),x_dot_valve(i),p);

                %% Leakage Flowrates (Leakage Model)
                [dm_leak_in(i),dm_leak_out(i),Ma_cv2(i)] = leakage( P(i),P_cv2(i),x_dot_piston(i),T(i),rho(i),T_cv2(i),x_piston(i),p );


                %Total flowrates (Massflow Summations)
                if dm_valve<0     %flow IN
                    dm_in(i)=abs(dm_valve)+dm_leak_in(i);
                    dm_out(i)=dm_leak_out(i);
                    dm(i)=dm_in(i)-dm_out(i);
                    dm_out_valve(i)=0;

                elseif dm_valve>0 %flow OUT
                    dm_in(i)=dm_leak_in(i);
                    dm_out(i)=dm_valve+dm_leak_out(i);
                    dm(i)=dm_in(i)-dm_out(i);
                    dm_out_valve(i)=dm_valve;
                else
                    dm_in(i)=dm_leak_in(i);
                    dm(i)=dm_leak_in(i)-dm_leak_out(i);
                    dm_out(i)=dm_leak_out(i);
                    dm_out_valve(i)=0;
                end

                %% Instantaneous Heat Transfer in control volumes

                [ Q(i) ] = Ins_HT( T(i),rho(i),T_w(k),V(i),dV(i),x_piston(i),x_dot_piston(i),p);

                [ Q_cv2(i) ] = Ins_HT_cv2( T_cv2(i),rho_cv2(i),T_w_cv2(k),V_cv2(i),dV_cv2(i),x_dot_piston(i),Q_motor(i),p);
                
                
                

                %% Mass of the gas in the cylinder
                m(i)=rho(i)*V(i);
                m_cv2(i)=rho_cv2(i)*V_cv2(i);

                %% Fluid Properties for approximation
                if i==1
                    p.h_in=EOS(p.T_i,p.rho_i,'enthalpy','rho');
                end
                h_out=EOS(T(i),rho(i),'enthalpy','rho');
                h(i)=h_out;
                dP_dT_v(i)= dP_dT(T(i),rho(i));
                Cv(i)=EOS(T(i),rho(i),'Cv','rho');

                h_cv2(i)=EOS(T_cv2(i),rho_cv2(i),'enthalpy','rho');
                dP_dT_v_cv2(i)=dP_dT(T_cv2(i),rho_cv2(i));
                Cv_cv2(i)=EOS(T_cv2(i),rho_cv2(i),'Cv','rho');

                if strcmp(p.method,'RK4')==1
                    %% RK4 Approximations

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Runge-Kutta 4th Order Approximation for Temperature and
                    %Density
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    %%%%%
                    %1st Control Volume, Piston Cylinder
                    %%%%%

                    [dT(i),drho(i)]=DT_RK(t(i),T(i),rho(i),P(i),T_w(k),x_valve(i),x_dot_valve(i),x_piston(i),x_dot_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),P_cv2(i),T_cv2(i),theta_dot_dot(i),p);



                    k1=[p.t_step*dT(i),p.t_step*drho(i)];
                    [k2(1),k2(2)]=DT_RK( t(i)+0.5*p.t_step,T(i)+0.5*k1(1),rho(i)+0.5*k1(2),P(i),T_w(k),x_valve(i),x_dot_valve(i),x_piston(i),x_dot_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),P_cv2(i),T_cv2(i),theta_dot_dot(i),p);
                    k2=[p.t_step*k2(1),p.t_step*k2(2)];
                    [k3(1),k3(2)]=DT_RK( t(i)+0.5*p.t_step,T(i)+0.5*k2(1),rho(i)+0.5*k2(2),P(i),T_w(k),x_valve(i),x_dot_valve(i),x_piston(i),x_dot_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),P_cv2(i),T_cv2(i),theta_dot_dot(i),p);
                    k3=[p.t_step*k3(1),p.t_step*k3(2)];
                    [k4(1),k4(2)]=DT_RK( t(i)+p.t_step,T(i)+k3(1),rho(i)+k3(2),P(i),T_w(k),x_valve(i),x_dot_valve(i),x_piston(i),x_dot_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),P_cv2(i),T_cv2(i),theta_dot_dot(i),p);
                    k4=[p.t_step*k4(1),p.t_step*k4(2)];

                    T(i+1)=T(i)+(1/6)*(k1(1)+2*k2(1)+2*k3(1)+k4(1));
                    rho(i+1)=rho(i)+(1/6)*(k1(2)+2*k2(2)+2*k3(2)+k4(2));
                    %keyboard

                    %%%%%%
                    %2nd Control Volume, Temperature and Density
                    %%%%%%

                    %[ dT,drho ] = DT_RK_cv2( t,T_cv2,P,P_cv2,x_dot_piston_1,T,rho,rho_cv2,h,Q_cv2,x_piston_1,theta_1,theta_dot_1,x_piston_m_1)

                    [dT_cv2(i),drho_cv2(i)]=DT_RK_cv2(t(i),T_cv2(i),P(i),P_cv2(i),x_dot_piston(i),T(i),rho(i),rho_cv2(i),h(i),T_w_cv2(k),x_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),theta_dot_dot(i),p);

                    k1_cv2=[p.t_step*dT_cv2(i),p.t_step*drho_cv2(i)];
                    [k2_cv2(1),k2_cv2(2)]=DT_RK_cv2(t(i)+0.5*p.t_step,T_cv2(i)+0.5*k1_cv2(1),P(i),P_cv2(i),x_dot_piston(i),T(i),rho(i),rho_cv2(i)+0.5*k1_cv2(2),h(i),T_w_cv2(k),x_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),theta_dot_dot(i),p);
                    k2_cv2=[p.t_step*k2_cv2(1),p.t_step*k2_cv2(2)];
                    [k3_cv2(1),k3_cv2(2)]=DT_RK_cv2(t(i)+0.5*p.t_step,T_cv2(i)+0.5*k2_cv2(1),P(i),P_cv2(i),x_dot_piston(i),T(i),rho(i),rho_cv2(i)+0.5*k2_cv2(2),h(i),T_w_cv2(k),x_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),theta_dot_dot(i),p);
                    k3_cv2=[p.t_step*k3_cv2(1),p.t_step*k3_cv2(2)];
                    [k4_cv2(1),k4_cv2(2)]=DT_RK_cv2(t(i)+p.t_step,T_cv2(i)+k3_cv2(1),P(i),P_cv2(i),x_dot_piston(i),T(i),rho(i),rho_cv2(i)+k3_cv2(2),h(i),T_w_cv2(k),x_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),theta_dot_dot(i),p);
                    k4_cv2=[p.t_step*k4_cv2(1),p.t_step*k4_cv2(2)];

                    T_cv2(i+1)=T_cv2(i)+(1/6)*(k1_cv2(1)+2*k2_cv2(1)+2*k3_cv2(1)+k4_cv2(1));
                    rho_cv2(i+1)=rho_cv2(i)+(1/6)*(k1_cv2(2)+2*k2_cv2(2)+2*k3_cv2(2)+k4_cv2(2));



                elseif strcmp(p.method,'Euler')==1

                    %% Simple Euler Approximation

                    %%%%%%%%%%%%%%%%%
                    %Euler Approximation, Temperature and Density
                    %%%%%%%%%%%%%%%%%

                    %%%%%
                    %1st Control Volume, Piston Cylinder
                    %%%%%
                    
                    
                    %dT=(Q+T*dP_dT_v*(1/1000)*((1/rho)*dm-dV)-h*dm+dm_in*p.h_in-dm_out*h)/(m*Cv);
                    %drho=(dm-rho*dV)/V;
                    [dT(i),drho(i)]=DT_RK(t(i),T(i),rho(i),P(i),T_w(k),x_valve(i),x_dot_valve(i),x_piston(i),x_dot_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),P_cv2(i),T_cv2(i),theta_dot_dot(i),p);

                    %dT(i)=(Q(i)+T(i)*dP_dT_v(i)*(1/1000)*((1/rho(i))*dm(i)-dV(i))-h(i)*dm(i)+dm_in(i)*p.h_in-dm_out(i)*h(i))/(m(i)*Cv(i));
                    %drho(i)=(dm(i)-rho(i)*dV(i))/V(i);
                    
                    T(i+1)=T(i)+p.t_step*dT(i);
                    rho(i+1)=rho(i)+p.t_step*drho(i);

                    %%%%%%
                    %2nd Control Volume
                    %%%%%%
                    
                    dm_cv2=dm_leak_out(i)-dm_leak_in(i);
                    dm_in_cv2=dm_leak_out(i);
                    dm_out_cv2=dm_leak_in(i);
                    

                    %[dT_cv2(i),drho_cv2(i)]=DT_RK_cv2(t(i),T_cv2(i),P(i),P_cv2(i),x_dot_piston(i),T(i),rho(i),rho_cv2(i),h(i),T_w_cv2(k),x_piston(i),theta(i),theta_dot(i),x_piston_m(i),dP(i),dx(i),theta_dot_dot(i),p);

                    dT_cv2(i)=(Q_cv2(i)+T_cv2(i)*dP_dT_v_cv2(i)*(1/1000)*((1/rho_cv2(i))*dm_cv2-dV_cv2(i))-h_cv2(i)*dm_cv2+dm_in_cv2*h(i)-dm_out_cv2*h_cv2(i))/(m_cv2(i)*Cv_cv2(i));
                    drho_cv2(i)=(dm_cv2-rho_cv2(i)*dV_cv2(i))/V_cv2(i);
                    
                    T_cv2(i+1)=T_cv2(i)+p.t_step*dT_cv2(i);
                    rho_cv2(i+1)=rho_cv2(i)+p.t_step*drho_cv2(i);


                else
                    i=floor(p.Period/p.t_step)-1; 
                    disp('Method error, please select a correct ODE Solution Method')
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Calculate next step for pressure
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                P(i+1)=EOS(T(i+1),rho(i+1),'pressure','rho');
                P_cv2(i+1)=EOS(T_cv2(i+1),rho_cv2(i+1),'pressure','rho');

                %%%%%%%%%%%%%
                %% Misc other properties
                %%%%%%%%%%%%
                dP(i+1)=P(i+1)-P(i);                    %Change in pressure
                dx(i+1)=x_piston(i+1)-x_piston(i);      %change in piston position
 

                %%%%%%%%%%%%
                %% Loop Adjustments
                %%%%%%%%%%%%
                if k>1
                    if i>1
                        %t_cross=find_cross(t',x_piston_m(1:end-1)',(abs(max(x_piston_m))-abs(min(x_piston_m)))/2);
                        t_cross=find_cross(t',x_piston_m(1:end-1)',0);
                    end
                else
                    if i>1
                        t_cross=find_cross(t',x_piston_m(1:end-1)',0);
                    end
                end
                
                t(i+1)=t(i)+p.t_step;
                i=i+1;
                
                if length(t)>=p.n*p.n_period*1.5
                    t_cross=ones(2*p.n_period);
                    disp('Convergence Error, piston did not cross zero twice')
                end

            end     %End of While loop checking t_cross
            
            
            %%%%%%%%%%%%%%%%%
            %% Adjust loop vectors
            %%%%%%%%%%%%%%%%%
            
            %F_drive = F_drive(1:length(t));
            %x_piston_m = x_piston_m(1:length(t));
            
            %%%%%%%%%%%%%%%
            %% Adjustments for Time
            %%%%%%%%%%%%%%%
            
            i=i-1;
            t=t(1:end-1);
            temp(k)=max(t_cross)-max(find_cross(t',F_drive',0));
            p.t_force_adjust=p.t_force_adjust+temp(k);
            t_change(k)=p.t_force_adjust;
            p.Period=t(end)-t(1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Estimate Massflow, Power Input,Heat Transfer, and Efficiencies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            

%             m_dot(k)=(sum(dm_out_valve))/((p.n_period-1)*length(t)); %#ok<SAGROW>
%             m_dot_leak_out(k)=(sum(dm_leak_out))/(length(t)*(p.n_period-1)); %#ok<SAGROW>
%             m_dot_leak_in(k)=(sum(dm_leak_in))/(length(t)*(p.n_period-1)); %#ok<SAGROW>
            
            %P_applied=(2*max(F_drive)*p.x_stroke)/p.Period;     %Calculated Applied Force

            m_dot(k)=(trapz(dm_out_valve)*(t(2)-t(1)))/p.Period; %#ok<SAGROW>
            m_dot_leak_out(k)=(trapz(dm_leak_out)*(t(2)-t(1)))/p.Period; %#ok<SAGROW>
            m_dot_leak_in(k)=(trapz(dm_leak_in)*(t(2)-t(1)))/p.Period; %#ok<SAGROW>
            
            
            W_dot(k)=m_dot(k)*(max(h)-p.h_in); %#ok<SAGROW>
            eta_o(k)=(m_dot(k)*(p.h_2_s-p.h_in)*1000)/p.P_electric;%#ok<SAGROW>      %W_dot is in kW
            eta_vol(k)=(m_dot(k))/(p.rho_i*f_list(l)*(max(x_piston_m)-min(x_piston_m))*p.Ap);%#ok<SAGROW>
            
            
%             Q_dot(k)=(sum(Q))/((p.n_period-1)*length(t));%#ok<SAGROW>
%             Q_dot_cv2(k)=(sum(Q_cv2))/((p.n_period-1)*length(t));%#ok<SAGROW>
            
            Q_dot(k)=(trapz(Q)*(t(2)-t(1)))/p.Period;%#ok<SAGROW>
            Q_dot_cv2(k)=(trapz(Q_cv2)*(t(2)-t(1)))/p.Period;%#ok<SAGROW>
            
            T_w(k+1)=(Q_dot(k)+Q_dot_cv2(k))*1000*p.R_shell+p.T_amb;
            T_w_cv2(k+1)=T_w(k);    %one lumped shell temp assumption
            

            %% Iteration Counter
            count1=num2str(k);
            count2='Iteration Number ';
            count=[count2 count1];
            disp(' ')
            disp(count)
            loop(k)=k;
            
            %%%%%%%%%%%%%%
            %% Residuals
            %%%%%%%%%%%%%%
            
            dT_res(k)=(T(1)-T(i));
            drho_res(k)=(rho(1)-rho(i));
            %T_1(k)=T(1);
            %rho_1(k)=rho(1);
            
            if k>1
            dx_res(k)=x_stroke_save(k)-x_stroke_save(k-1);
            end

            
            dT_loop(k)=dT_res(k)/min(T(1:length(t)));
            dP_loop=(P(1)-P(i))/min(P(1:length(t)));
            drho_loop(k)=(rho(1)-rho(i))/min(rho(1:length(t)));
            dT_loop_k=dT_loop(k);
            drho_loop_k=drho_loop(k);
            CV1_res=strcat('CV1 Res.:| ',' DT| ',num2str(dT_loop(k)),' Drho| ',num2str(drho_loop(k)),' dP| ',num2str(dP_loop));
            %disp(CV1_res);

            if p.leakage_on==1
                dT_res_cv2(k)=(T_cv2(1)-T_cv2(i));
                drho_res_cv2(k)=(rho_cv2(1)-rho_cv2(i));
                %T_1_cv2(k)=T_cv2(1);
                %rho_1_cv2(k)=rho_cv2(1);
                
                dT_cv2_loop(k)=(T_cv2(1)-T_cv2(i))/min(T_cv2(1:length(t)));
                drho_cv2_loop(k)=(rho_cv2(1)-rho_cv2(i))/min(rho_cv2(1:length(t)));
                
                dT_cv2_loop_k=dT_cv2_loop(k);
                drho_cv2_loop_k=drho_cv2_loop(k);
                
                dP_cv2_loop=(P_cv2(1)-P_cv2(i))/min(P_cv2(1:length(t)));
                CV2_res=strcat('CV2 Res.:| ',' DT| ',num2str(dT_cv2_loop(k)),' Drho| ',num2str(drho_cv2_loop(k)),' dP| ',num2str(dP_cv2_loop));
                %disp(CV2_res);

            else
                dT_cv2_loop=0;
                drho_cv2_loop=0;
                dP_cv2_loop=0;
                dT_res_cv2(k)=0;
                drho_res_cv2(k)=0;
                dT_cv2_loop_k=0;
                drho_cv2_loop_k=0;
            end     
            
            if k>1
                dx_piston=(x_stroke_save(k)-x_stroke_save(k-1))/max(x_stroke_save);
                %disp(dx_piston)
                dm_dot_cycle=(m_dot(k)-m_dot(k-1))/max(m_dot);
                %disp(dm_dot_cycle)
                if isnan(dm_dot_cycle)
                    dm_dot_cycle=1;
                end
                dT_w=(T_w(k)-T_w(k-1))/max(T_w);
                Other_res=strcat('RES:|','dx:| ',num2str(dx_piston),' dm:| ',num2str(dm_dot_cycle),' dT_w:| ',num2str(dT_w));
                %disp(Other_res)
            end
            
            
            %%%%%%%%%%%%%%%%
            %Plot Residuals
            %%%%%%%%%%%%%%%%
%             if k>2
%             figure(1)
%             set(gcf,'Visible','off')
%             handles.plot(1)=plot(loop,dT_loop,'bo-'); 
%             hold on
%             plot(k,dP_loop,'ro-');
%             plot(loop,drho_loop,'go-');
%             plot(loop,dT_cv2_loop,'co-');
%             plot(loop,drho_cv2_loop,'rx-');
%             plot(k,dx_piston,'ko-');
%             
%             
%             title('Residuals');xlabel('Iterations');ylabel('Residual Value');
%             
%             if k==2
%                 legend('Temperature','Pressure','density','Temp CV2','Density CV2','Piston','Location','SouthEastOutside');
%             end
%             handles.figure(1)=gcf;
%             handles.axis(1)=gca;
%             %set(handles.figure(1),'Position',get(0,'ScreenSize'))
%             set(handles.axis(1),'YLim',[-0.01,0.01])
%             
%             end
%             
            if pc_flag == 0
                figure(30)
                plot(t,x_piston_m(1:end-1))
                hold on
            end

            if k>p.loop_error+5
                


            else
                X_n=[T(i),rho(i),T_cv2(i),rho_cv2(i),x_stroke_save(k)];
                F_n=[dT(i),drho(i),dT_cv2(i),drho_cv2(i),x_dot_piston(i)];
                J=0;
                DELTAX_n=0;
                X_n_1=X_n+DELTAX_n;
            end
            
            %Error Loop Counter
            if k==p.loop_error
                dT_loop_k=0;
                drho_loop_k=0;
                count='Convergence Error, too many iterations';
                disp(count)
                dT_cv2_loop_k=0;
                drho_cv2_loop_k=0;
                dx_piston=0;
                dm_dot_cycle=0;
            end
            
            
            iter(k)=k;
            k=k+1;
            
            p.dP_max=max(P)-p.P_i;        %Max change in pressure over cycle relative to cylinder
            
            if k == 30
                error = error + 0.00005;
            elseif k == 60
                error = error + 0.00005;
            end
            
            %% Adjustments to power to try and maintain a correct stroke
            
            %If stroke is outside of the range of 0.001 and 0.0254 m either
            %add or substract 1W of power input to try and get it into the
            %appropriate range
            if k > 6
                if p.x_stroke < 0.001
                    p.P_electric = p.P_electric +1;
                    disp('Increase P_electric =')
                    disp(p.P_electric)
                elseif p.x_stroke > 0.0254
                    p.P_electric = p.P_electric -1;
                end
            end
            
            

        end                                 %End of k while loop
        
        %%%%%%%%%%%%%%
        %% Post-Allocation, adjust for proper length
        %%%%%%%%%%%%%%

        P=P(1:length(t));
        dT=dT(1:length(t));
        h=h(1:length(t));
        V=V(1:length(t));
        dV=dV(1:length(t));
        dm=dm(1:length(t));
        dm_in=dm_in(1:length(t));
        dm_out=dm_out(1:length(t));
        drho=drho(1:length(t));
        m=m(1:length(t));
        rho=rho(1:length(t));
        x_piston=x_piston(1:length(t));
        x_dot_piston=x_dot_piston(1:length(t));
        Cv=Cv(1:length(t));
        %iter=iter(1:length(t));
        dP_dT_v=dP_dT_v(1:length(t));
        Ma=Ma(1:length(t));

        T=T(1:length(t));
        x_valve=x_valve(1:length(t));
        x_dot_valve=x_dot_valve(1:length(t));
        dm_leak_in=dm_leak_in(1:length(t));
        dm_leak_out=dm_leak_out(1:length(t));
        dm_out_valve=dm_out_valve(1:length(t));
        dV_cv2=dV_cv2(1:length(t));
        V_cv2=V_cv2(1:length(t));
        Ma_cv2=Ma_cv2(1:length(t));
        h_cv2=h_cv2(1:length(t));
        Cv_cv2=Cv_cv2(1:length(t));
        dP_dT_v_cv2=dP_dT_v_cv2(1:length(t));
        P_cv2=P_cv2(1:length(t));
        drho_cv2=drho_cv2(1:length(t));
        rho_cv2=rho_cv2(1:length(t));
        dT_cv2=dT_cv2(1:length(t));
        T_cv2=T_cv2(1:length(t));
        F_drive=F_drive(1:length(t));
        Q_motor=Q_motor(1:length(t));
        theta=theta(1:length(t));
        theta_dot=theta_dot(1:length(t));
        x_piston_m=x_piston_m(1:length(t));
        dP=dP(1:length(t));
        dx=dx(1:length(t));
        theta_dot_dot=theta_dot_dot(1:length(t));
        Q=Q(1:length(t));
        Q_cv2=Q_cv2(1:length(t));
        F_wall=F_wall(1:length(t));
        theta_1_c=theta_1_c(1:length(t));
        k_eff=k_eff(1:length(t));
        
        T_w = T_w(1:k-1);
        T_w_cv2 = T_w_cv2(1:k-1);
        x_stroke_save = x_stroke_save(1:k-1);
        loop = loop(1:k-1);
        dT_res = dT_res(1:k-1);
        drho_res = drho_res(1:k-1);
        dx_res = dx_res(1:k-1);
        dT_res_cv2 = dT_res_cv2(1:k-1);
        drho_res_cv2 = drho_res_cv2(1:k-1);
        dT_loop = dT_loop(1:k-1);
        drho_loop = drho_loop(1:k-1);
        dT_cv2_loop = dT_cv2_loop(1:k-1);
        drho_cv2_loop = drho_cv2_loop(1:k-1);    

        
        x_mag(l)=p.x_stroke/p.x_max; %#ok<SAGROW>
        x_stroke_resonant_save(l) = (max(x_piston_m)-min(x_piston_m)); %#ok<SAGROW>
        f_list_resonant_save(l) = f_list(l); %#ok<SAGROW>
        
        %Converts power input to force output by motor, based on motor
        %literature
        k_m = -0.0008*(p.P_electric*p.eta_motor)^2 - 0.0047*(p.P_electric*p.eta_motor) + 4.6831;
        F_drive_max=sqrt(p.P_electric*p.eta_motor)*(k_m);       %Maximum Motor force, in Newtons.
        
        W_dot_estimated = F_drive_max*x_stroke_resonant_save(l)*f_list(l)*(2/pi);
        
    %% Calculating resonant frequency
    
    if l>1
        
        change_in_stroke = x_stroke_resonant_save(l) - x_stroke_resonant_save(l-1);
        
        if change_in_stroke < 0
            number_change_stroke = number_change_stroke +1;
        end
        
        if number_change_stroke < 3                 
            
            p.resonant_loop = 0;
            f_resonant_calculated = 0;
            w_d_resonant_calculated = 0;
            
        elseif number_change_stroke >= 3      % is stroke decreasing
            
            p.resonant_loop = 1;
            
            if length(x_stroke_resonant_save) < 5
                
                f_resonant_calculated = 99;
                w_d_resonant_calculated = 99;
                
            elseif length(x_stroke_resonant_save) <= 15 && length(x_stroke_resonant_save) >=5
                
                coefficients = polyfit(f_list_resonant_save,x_stroke_resonant_save,2);            
                f_resonant_calculated = -coefficients(2)/(2*coefficients(1));
                w_d_resonant_calculated = f_resonant_calculated*2*pi;                
                
            else
                
                x_stroke_resonant_save = x_stroke_resonant_save(end-14:end);
                f_list_resonant_save = f_list_resonant_save(end-14:end);
                coefficients = polyfit(f_list_resonant_save,x_stroke_resonant_save,2);            
                f_resonant_calculated = -coefficients(2)/(2*coefficients(1));
                w_d_resonant_calculated = f_resonant_calculated*2*pi;  
                
            end
                                      
        end
        
    else
        f_resonant_calculated = 0;
        w_d_resonant_calculated = 0;
        number_change_stroke = 0;           %initialize variable
        p.resonant_loop = 0;
        
    end
    
        

    %% Saving Important parameters
    
    if pc_flag == 0 %&& p.resonant_loop == 1
        
        d=date;
        save_data={d,char(p.save_name),z,p.w_d,f_list(l),(max(x_piston_m)-min(x_piston_m)),W_dot(k-1),...
            m_dot(k-1),m_dot_leak_out(k-1),p.P_i,p.P_d,p.T_i,length(t),char(p.method),p.P_electric,...
            p.eta_motor,p.k_mech,p.f_friction,p.ecc_1,x_mag(l),p.g,eta_o(k-1),eta_vol(k-1),max(T),...
            T_w(k-1),k-1,Q_dot(k-1),Q_dot_cv2(k-1),p.P_d/p.P_i,dT_loop(k-1),drho_loop(k-1),dT_cv2_loop(k-1),drho_cv2_loop(k-1),...
            dx_piston,dm_dot_cycle,dT_w,p.M_mov,(max(x_piston_m)-min(x_piston_m)),p.P_electric,f_resonant_calculated,w_d_resonant_calculated,p.k_eff};

        range_2=strcat('A',num2str(l+rows),':','AP',num2str(l+rows));
        xlswrite1('data_log.xls',save_data,range_2)
    
    elseif pc_flag == 1
        
        save_data_tmp = {z,p.w_d,f_list(l),(max(x_piston_m)-min(x_piston_m)),W_dot(k-1),...
            m_dot(k-1),m_dot_leak_out(k-1),p.P_i,p.P_d,p.T_i,length(t),p.P_electric,...
            p.eta_motor,p.k_mech,p.f_friction,p.ecc_1,x_mag(l),p.g,eta_o(k-1),eta_vol(k-1),max(T),...
            T_w(k-1),k-1,Q_dot(k-1),Q_dot_cv2(k-1),p.P_d/p.P_i,dT_loop(k-1),drho_loop(k-1),dT_cv2_loop(k-1),drho_cv2_loop(k-1),...
            dx_piston,dm_dot_cycle,dT_w,p.M_mov,(max(x_piston_m)-min(x_piston_m)),p.P_electric,f_resonant_calculated,w_d_resonant_calculated,p.k_eff};
        
        
        save_data(z,:) = cell2mat(save_data_tmp);
        
        save_name = strcat('data_log_',p.save_name);
        save_name = char(save_name);
        xlswrite(save_name,save_data);
        
    end
    
    k=1;
    dT_loop=1;
    drho_loop=1;
    dT_cv2_loop=1;
    drho_cv2_loop=1;
    
    l=l+1;              %Resonant loop counter update
    
end     %End of freq sweep while loop, l variable

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Diagrams
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    if pc_flag == 0
        
    %%%%%%%%%%%%%%
    %Cylinder Plots
    %%%%%%%%%%%%%%
    
%     figure(2)
%     set(gcf,'Visible','off')
%     handles.plot(2)=plot(t,T); title('Temperature v. Time - Cylinder');
%     ylabel('Temperature (K)');
%     handles.figure(2)=gcf;
%     handles.axis(2)=gca;
%     xlabel('Time (sec)')
%     
%     figure(3)
%     set(gcf,'Visible','off')
%     handles.plot(3)=plot(t,P); title('Pressure v. Time - Cylinder');
%     ylabel('Pressure (kPa)');
%     handles.figure(3)=gcf;
%     handles.axis(3)=gca;
%     xlabel('Time (sec)')
%     
%     figure(4)
%     set(gcf,'Visible','off')
%     handles.plot(4)=plot(t,dm); title('Massflow v. Time - Cylinder');
%     ylabel('Massflowrate (kg/sec)');
%     handles.figure(4)=gcf;
%     handles.axis(4)=gca;
%     xlabel('Time (sec)')
%     
%     figure(5)
%     set(gcf,'Visible','off')
%     handles.plot(5)=plot(t,Q); title('Heat Transfer v. Time - Cylinder');
%     ylabel('Heat Transfer (kW)');
%     handles.figure(5)=gcf;
%     handles.axis(5)=gca;
%     xlabel('Time (sec)')
%     
%     figure(6)
%     set(gcf,'Visible','off')
%     handles.plot(6)=plot(t,rho); title('Density v. Time - Cylinder');
%     ylabel('Density (kg/m^3)');
%     handles.figure(6)=gcf;
%     handles.axis(6)=gca;
%     xlabel('Time (sec)')
%     
%     figure(7)
%     set(gcf,'Visible','off')
%     handles.plot(7)=plot(t,h); title('Enthalpy v. Time - Cylinder');
%     ylabel('Enthalpy (kJ/kg)');
%     handles.figure(7)=gcf;
%     handles.axis(7)=gca;
%     xlabel('Time (sec)')
%     
%     figure(8)
%     set(gcf,'Visible','off')
%     handles.plot(8)=plot(t,dV); title('Change in Volume v. Time - Cylinder');
%     ylabel('Change in Volume (m^3/sec)');
%     handles.figure(8)=gcf;
%     handles.axis(8)=gca;
%     xlabel('Time (sec)')
%     
%     figure(9)
%     %set(gcf,'Visible','off')
%     handles.plot(9)=plot(t,x_piston); title('Piston Position v. Time - Cylinder');
%     ylabel('Piston Position (m)');
%     handles.figure(9)=gcf;
%     handles.axis(9)=gca;
%     xlabel('Time (sec)')
%     
%     figure(10)
%     set(gcf,'Visible','off')
%     handles.plot(10)=plot(t,drho); title('Change in Density v. Time - Cylinder');
%     ylabel('Change in Density (kg/m^3/sec)');
%     handles.figure(10)=gcf;
%     handles.axis(10)=gca;
%     xlabel('Time (sec)')
%     
%     %%%%%%%%%%%%%%%%%
%     %Valve and Leakage Plots
%     %%%%%%%%%%%%%%%%%
%     
%     figure(11)
%     set(gcf,'Visible','off')
%     handles.plot(11)=plot(t,x_valve); title('Valve Lift v. Time');
%     ylabel('Valve Lift (m)');
%     handles.figure(11)=gcf;
%     handles.axis(11)=gca;
%     xlabel('Time (sec)')
%     
%     figure(12)
%     set(gcf,'Visible','off')
%     handles.plot(12)=plot(t,Ma); title('Mach Number v. Time');
%     ylabel('Mach Number (-)');
%     handles.figure(12)=gcf;
%     handles.axis(12)=gca;
%     xlabel('Time (sec)')
%     
%     figure(13)
%     set(gcf,'Visible','off')
%     handles.plot(13)=plot(t,x_valve); title('Valve Lift v. Time');
%     ylabel('Valve Lift (m)');
%     handles.figure(13)=gcf;
%     handles.axis(13)=gca;
%     xlabel('Time (sec)')
%     
%     figure(14)
%     set(gcf,'Visible','off')
%     handles.plot(14)=plot(t,dm_leak_out); title('Leakage Massflow out of Cylinder v. Time');
%     ylabel('Flowrate (kg/sec)');
%     handles.figure(14)=gcf;
%     handles.axis(14)=gca;
%     xlabel('Time (sec)')
%     
%     figure(15)
%     set(gcf,'Visible','off')
%     handles.plot(15)=plot(t,Ma_cv2); title('Mach Number Past Cylinder v. Time');
%     ylabel('Mach Number (-)');
%     handles.figure(15)=gcf;
%     handles.axis(15)=gca;
%     xlabel('Time (sec)')
%   
%     
%   
% 
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Plots for Full Iterations
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     figure(16)
%     set(gcf,'Visible','off')
%     handles.plot(16)=plot(m_dot); title('Total Massflow v. Iteration');
%     ylabel('Massflow Rate (kg/sec)');
%     handles.figure(16)=gcf;
%     handles.axis(16)=gca;
%     xlabel('Iterations (-)')
%     
%     figure(17)
%     set(gcf,'Visible','off')
%     handles.plot(17)=plot(W_dot); title('Work Done on Gas v. Iteration');
%     ylabel('Work (kW)');
%     handles.figure(17)=gcf;
%     handles.axis(17)=gca;
%     xlabel('Iterations (-)')
%     
%     figure(18)
%     set(gcf,'Visible','off')
%     handles.plot(18)=plot(m_dot_leak_out); title('Leakage Out of Cylinder v. Iteration');
%     ylabel('Massflow Rate (kg/sec)');
%     handles.figure(18)=gcf;
%     handles.axis(18)=gca;
%     xlabel('Iterations (-)')
%     
%     figure(19)
%     set(gcf,'Visible','off')
%     handles.plot(19)=plot(eta_o); title('Efficienies v. Iteration');
%     ylabel('Efficiency (-)');
%     handles.figure(19)=gcf;
%     handles.axis(19)=gca;
%     hold on
%     handles.plot(20)=plot(eta_vol,'r');
%     xlabel('Iterations (-)')
%     handles.figure(20)=gcf;
%     handles.axis(20)=gca;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % P-V Plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    figure(21)
    set(gcf,'Visible','off')
    handles.plot(21)=plot(V,P); title('Pressure v. Volume');
    ylabel('Pressure (kPa)');
    xlabel('Volume (m^3)');
    handles.figure(21)=gcf;
    handles.axis(21)=gca;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Diagnostic Plots
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(22)
    set(gcf,'Visible','off')
    handles.plot(22)=plot(Q_dot); title('Heat Transfer v. Iteration');
    ylabel('Heat Transfer (kW)');
    xlabel('Iteration');
    handles.figure(22)=gcf;
    handles.axis(22)=gca;
    
    figure(23)
    set(gcf,'Visible','off')
    handles.plot(23)=plot(Q_dot_cv2); title('Heat Transfer v. Iteration - CV2');
    ylabel('Heat Transfer (kW)');
    xlabel('Iteration');
    handles.figure(23)=gcf;
    handles.axis(23)=gca;
    
    figure(24)
    set(gcf,'Visible','off')
    handles.plot(24)=plot(T_w); title('Shell Wall Temp v. Iteration - CV2');
    ylabel('Wall Temp (K)');
    xlabel('Iteration');
    handles.figure(24)=gcf;
    handles.axis(24)=gca;
    

    %Main subplot diagram
    figure(25)
    set(gcf,'Visible','off')
    subplot(3,3,1),plot(t,T);title('Temperature v. Time');
    subplot(3,3,2),plot(t,P);title('Pressure v. Time');
    subplot(3,3,3),plot(t,dm);title('Massflow v. Time');
    subplot(3,3,4),plot(t,V);title('Volume v. Time');
    subplot(3,3,5),plot(t,V);title('Volume v. Time');
    subplot(3,3,6),plot(t,rho);title('Density v. Time');
    subplot(3,3,7),plot(t,dV);title('Change in Volume v. Time');
    subplot(3,3,8),plot(t,x_piston);title('Displacement v. Time');
    subplot(3,3,9),plot(t,drho);title('Change in Rho v. Time');
    
% 
    figure(26)
    set(gcf,'Visible','off')
    %plot(t,dT_check);hold on;plot(t,dT,'r');legend('dT_c_h_e_c_k','dm')
    subplot(2,2,1),plot(t,Ma);title('Mach Number v. Time')
    subplot(2,2,2),plot(t,x_valve);title('Valve Lift v. Time')
    subplot(2,2,3),plot(m_dot_leak_out);title('Leakage Outflow v. Time')
    subplot(2,2,4),plot(m_dot);title('Massflow v. Iteration')
% 
% 
    % CV2 Plots
    figure(27)
    set(gcf,'Visible','off')
    subplot(3,3,1),plot(t,T_cv2);title('Temperature v. Time for CV2');
    subplot(3,3,2),plot(t,P_cv2);title('Pressure v. Time for CV2');
    subplot(3,3,3),plot(t,dm_leak_in);title('Leakage Massflow In v. Time for CV2');
    subplot(3,3,4),plot(t,V_cv2);title('Volume v. Time for CV2');
    subplot(3,3,5),plot(t,rho_cv2);title('Density v. Time');
    subplot(3,3,6),plot(t,h_cv2);title('Enthalpy v. Time');
    subplot(3,3,7),plot(t,dT_cv2);title('Change in Temperature v. Time for CV2');
    subplot(3,3,8),plot(t,Ma_cv2);title('Mach Number v. Time');
    subplot(3,3,9),plot(t,drho_cv2);title('Change in Rho v. Time for CV2');
    

    %%%%%%%%%%%%%%%%%%%%%%%%
    %% Format and Save all Plots
    %%%%%%%%%%%%%%%%%%%%%%%%
    
%  
%     if p.save_all==1
%        for w=1:1:length(handles.plot)
% 
%             set(handles.axis(w),'FontSize',14)
%             handles.title(w)=get(handles.axis(w),'title');
%             handles.xlabel(w)=get(handles.axis(w),'xlabel');
%             handles.ylabel(w)=get(handles.axis(w),'ylabel');
%             
%             set(handles.title(w),'FontSize',16)
%             set(handles.xlabel(w),'FontSize',14)
%             set(handles.ylabel(w),'FontSize',14)
% 
%             set(handles.plot(w),'LineWidth',2)
% 
%             print (handles.figure(w),'-dpng','-r300')
%             print (handles.figure(w),'-depsc','-r300')
% 
%         end
%     end
%     
%     
    
    
    
        %% Shutting down activex server to excel
        invoke(Excel.ActiveWorkbook,'Save');

        Excel.Quit
        Excel.delete
        clear Excel   
        
    end
    
    %% other stuff
    
%     if p.save_all==1
%         %Save the current workspace
%         workspace_name=strcat(p.save_name,'_workspace.mat');
%         save (char(workspace_name))
%         
%         
%         if pc_flag == 0
%             %Make a directory based on the current save name
%             path_name=strcat(pwd,'\batch\',date,'\',p.save_name);
%             temp_2=char(path_name);
%         
%         if pc_flag == 1
%             
%             path_name=strcat(pwd,'/batch/',date,'/',p.save_name);
%             temp_2=char(path_name);
%         
%         mkdir('batch')    
%         cd batch
%         mkdir(date)
%         cd (date)
%         temp=char(p.save_name);
%         mkdir(temp)
%         cd ..
%         cd ..
%         
%         %zip everything up and save it
%         zip_name=strcat(p.save_name,'.zip');
%         temp=char(zip_name);
%         zip(temp,{'*.m','*.fig','*.mat','*.png','*.eps','*.xls'})
%         
%         
%         cd batch
%         cd (date)
%         cd (char(p.save_name))
%         unzip(temp,temp_2)
%         
%         cd ..
%         cd ..
%         cd ..
%         
%         delete(temp)
        delete('*.png','*.eps','*.mat','*.asv')
 %   end
    
%    disp(strcat('Batch Line No: ',num2str(z)))
   
end

%turn annoying warning back on
warning('on','MATLAB:Print:SavingToDifferentName')

disp('finished')

catch ME

    disp('finished w/ error')
    disp(ME.message)
    disp(ME.cause)
    disp(ME.stack)
    disp(ME)
    
    if pc_flag == 0
        %% Shutting down activex server to excel
        invoke(Excel.ActiveWorkbook,'Save');

        Excel.Quit
        Excel.delete
        clear Excel   
        
    end
    
    if pc_flag == 1

        save_name = strcat('data_log_',p.save_name);
        save_name = char(save_name);
        xlswrite(save_name,save_data);

    end


%turn annoying warning back on
warning('on','MATLAB:Print:SavingToDifferentName')



end














































