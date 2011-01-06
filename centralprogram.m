
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

% 5 - Add compressibility to leakage model


try

    clc
    clear all
    close all

    %% Initializations

    %pc_flag = input('Is this run on a server or PC without Excel COM Server? (1 - yes, 0 - no)');
    %pc_flag = 1;
    %pc_flag = 0;
    
    pc_flag=getflag();
    
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
        %analy_length = (input_data(4,1) - input_data(3,1))+1; 
        analy_length = (input_data(2,1) - input_data(1,1))+1; 

    else
        error('Invalid pc_flag.  Please use either 0 or 1.');
        clear all

    end

    %turn annoying warning off
    warning('off','MATLAB:Print:SavingToDifferentName')

    %%%%%%%%%%%%%
    % For loop used for batch operation, line by line from input file, z
    % loop
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
                d_name = strcat(txt_inputs{56},'.txt');
                diary(d_name)
            end

            %num_inputs = input_data(z+2,2:61);
            num_inputs = input_data(z,2:61);

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
    % Loop used for frequency sweeps, l loop
    %%%%%%%%%

    p.resonant_loop = 0;        %initialize variable that tells frequency sweep loop if we have reached the maximum frequency.
    l=1;                        %loop variable for the following loop.

    while p.resonant_loop<1 && l<=length(w_d_list)
    %for l=1:1:length(w_d_list)

        close all
        p.w_d=w_d_list(l);
        p.Period=1/f_list(l);     %Time of one cycle
        p.t_step=1/(f_list(l)*p.n);       %time step in seconds
        p.transient=0;
        t_adjust_force=0;
%         k=1;
%         dT_loop=1;
%         drho_loop=1;
%         dT_cv2_loop=1;
%         drho_cv2_loop=1;
%         dT_loop_k=1;
%         drho_loop_k=1;
%         dT_cv2_loop_k=1;
%         drho_cv2_loop_k=1;
%         dx_piston=1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%    
        %Brents Method estimations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        a = 1;              %loop iteration variable
        p.P_brent(1) = 10;   %lower power guess
        p.P_brent(2) = 15; %second lower power guess
        error_brents = 0.00000005;    %0.05 microns
        error_stroke = 1;
        
    while error_stroke > error_brents   %a loop, brents method to back out power for full stroke
            
        %Initialze residuals and convergence loop iteration counter
        k=1;
        dT_loop=1;
        drho_loop=1;
        dT_cv2_loop=1;
        drho_cv2_loop=1;
        dT_loop_k=1;
        drho_loop_k=1;
        dT_cv2_loop_k=1;
        drho_cv2_loop_k=1;
        dx_piston=1;

        if p.brent == 1    
            p.P_electric = p.P_brent(a);
        end
     

        %% Initalize Vectors
        %These are the vectors that are iterated on and are thus not included in
        %the parameters structure, i loop variables

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
        c_eff = zeros(1,p.n_period*p.n);
        c_friction = zeros(1,p.n_period*p.n);
        c_gas = zeros(1,p.n_period*p.n);
        x_dot_dot_piston = zeros(1,p.n_period*p.n);
        F_gas = zeros(1,p.n_period*p.n);

        % Initialize convergence loop variables, k loop
        x_piston_m=0;
        T_w = zeros(1,p.loop_error);
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
        
        m_dot = zeros(1,p.loop_error);
        m_dot_leak_in = zeros(1,p.loop_error);
        m_dot_leak_out = zeros(1,p.loop_error);
        W_dot = zeros(1,p.loop_error);
        eta_o = zeros(1,p.loop_error);
        eta_vol = zeros(1,p.loop_error);
        Q_dot = zeros(1,p.loop_error);
        Q_dot_cv2 = zeros(1,p.loop_error);
        m_dot_in = zeros(1,p.loop_error);
        

        loop=1;     %loop counter variable, re-initialize
        p.dP_max=p.P_d-p.P_i;
        P_electric_temp = p.P_electric;         %save input value
        

             %%%%%%%%%%%%
             %%  while loop to enforce convergence, k loop
             %%%%%%%%%%%%

            while (abs(dT_loop_k)>error || abs(drho_loop_k)>error) || (abs(dT_cv2_loop_k)>error || abs(drho_cv2_loop_k)>error) || abs(dx_piston)>(error/10)

%                 if k < 2
%                     p.P_electric = 100;
%                 elseif k >= 10 && k < 35
%                     delta_power = (100 - P_electric_temp)/25;
%                     p.P_electric = -delta_power*(k-10)+100;
%                 elseif k >= 35
%                     p.P_electric = P_electric_temp;
%                 end
                
                W_dot_in(k) = p.P_electric;

                %% Inital Cylinder Conditions
                if k==1

                    i=1;                                                    %initialize i loop counter
                    T(i)=p.T_i;                                             %initial temperature
                    rho(i)=p.rho_i;                                         %initial density
                    P(i)=EOS(T(i),rho(i),'pressure','rho');                 %calculate initial pressure
                    T_cv2(1)=p.T_i;                                         %initial temperature, CV2
                    rho_cv2(1)=p.rho_i;                                     %initial density, CV2
                    P_cv2(1)=P(1);
                    %P_cv2(1)=(p.P_d + p.P_i)/2;                             %initial pressure, CV2
                    %rho_cv2(1) = EOS(T_cv2(1),P_cv2(1),'rho','P');
                    x_piston(1)=p.x_piston_i;                               %initial piston position, default is (max stroke/2)
                    x_dot_piston(1)=p.x_dot_piston_i;                       %initial piston velocity, default is 0
                    V(1)=(p.x_piston_i)*p.Ap+p.V_dead_valve;                %initial volume
                    dV(1)=-p.x_dot_piston_i*p.Ap;                           %initial change in volume
                    theta(1)=p.theta_i;                                     %initial piston rotation
                    theta_dot(1)=p.theta_dot_i;                             %initial change in piston rotation                   
                    p.x_dot_ave=1;                                          %%%%%%IMPORTANT, AVERAGE SPEED
                    dP(1)=0;                                                %Change in pressure
                    dx(1)=1;                                                %Change in stroke
                    T_w(k)=p.T_w_i;                                         %initial compressor lump temperature (wall temp)
                    x_stroke_save(k)=p.x_stroke_initial;                    %saved value of stroke
                    p.x_stroke = p.x_stroke_initial;                        %guess stroke
                    p.dP_max = p.P_d-p.P_i;
                    t=0;
                    p.dV = 0;                                               %initialize dV
                    p.P_current = P(i);

                    %Updated guess values
%                     i=1;                                                    %initialize i loop counter
%                     T(i)=p.T_s_2;                                           %initial temperature
%                     P(i)=(EOS(T(i),rho(i),'pressure','rho')+p.P_d)/2;       %calculate initial pressure
%                     rho(i)=EOS(T(i),P(i),'rho','P');                        %initial density
%                     T_cv2(1)=T(i);                                         %initial temperature, CV2
%                     rho_cv2(1)=rho(i);                                     %initial density, CV2
%                     P_cv2(1)=P(1);                                          %initial pressure, CV2, same as CV1
%                     x_piston(1)=p.x_piston_i;                               %initial piston position, default is (max stroke/2)
%                     x_dot_piston(1)=p.x_max*f_list(l);                       %initial piston velocity, default is 0
%                     V(1)=(p.x_piston_i)*p.Ap+p.V_dead_valve;                %initial volume
%                     dV(1)=-x_dot_piston(i)*p.Ap;                           %initial change in volume
%                     theta(1)=p.theta_i;                                     %initial piston rotation
%                     theta_dot(1)=p.theta_dot_i;                             %initial change in piston rotation                   
%                     p.x_dot_ave=x_dot_piston(i)*0.707;                                          %%%%%%IMPORTANT, AVERAGE SPEED
%                     dP(1)=0;                                                %Change in pressure
%                     dx(1)=1;                                                %Change in stroke
%                     T_w(k)=p.T_w_i;                                         %initial compressor lump temperature (wall temp)
%                     x_stroke_save(k)=p.x_stroke_initial;                    %saved value of stroke
%                     p.x_stroke = p.x_stroke_initial;                        %guess stroke
%                     p.dP_max = p.P_d-p.P_i;
%                     t=0;                    
                    
%                 elseif k < 2    
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %Update values for intermediate iterations
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     
%                     T(1)=X_n_1(1);                                          %Temperatures and pressures equal to the end values of the previous iteration
%                     rho(1)=X_n_1(2);
%                     P(1)=EOS( T(1),rho(1),'pressure','rho' );
%                     T_cv2(1)=X_n_1(3);
%                     rho_cv2(1)=X_n_1(4);
%                     P_cv2(1)=EOS( T_cv2(1),rho_cv2(1),'pressure','rho' );
%                                      
%                     x_piston(1)=x_piston(i);
%                     x_dot_piston(1)=p.x_max*f_list(l);
%                     p.x_dot_ave=0.707*max(x_dot_piston);                    %%%%%%IMPORTANT, AVERAGE SPEED, RMS value of speed
%                     V(1)=V(i);
%                     dV(1)=-x_dot_piston(i)*p.Ap;
%                     theta(1)=theta(i);
%                     theta_dot(1)=theta_dot(i);                   
%                     t=t(i);
%                     
%                     p.x_stroke=max(x_piston_m)-min(x_piston_m);             %max stroke, current value
%                     x_stroke_save(k)=p.x_stroke;                            % stroke variable saved to see progression
%                     dP(1)=dP(i);
%                     dx(1)=dx(i);                                            %rest of above resets remaining variables to the ending previous values of the previous k loop
%                     i=1;                %resets inner loop variable
                else
  
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Update values for iterations beyond first
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    T(1)=X_n_1(1);                                          %Temperatures and pressures equal to the end values of the previous iteration
                    rho(1)=X_n_1(2);
                    P(1)=EOS( T(1),rho(1),'pressure','rho' );
                    T_cv2(1)=X_n_1(3);
                    rho_cv2(1)=X_n_1(4);
                    P_cv2(1)=EOS( T_cv2(1),rho_cv2(1),'pressure','rho' );
                    p.x_stroke=max(x_piston_m)-min(x_piston_m);             %max stroke, current value
                    p.x_dot_ave=0.707*max(x_dot_piston);                    %%%%%%IMPORTANT, AVERAGE SPEED, RMS value of speed
                    x_stroke_save(k)=p.x_stroke;                            % stroke variable saved to see progression
                    x_piston(1)=x_piston(i);
                    x_dot_piston(1)=x_dot_piston(i);
                    V(1)=V(i);
                    dV(1)=dV(i);
                    theta(1)=theta(i);
                    theta_dot(1)=theta_dot(i);                   
                    t=t(i);
                    dP(1)=dP(i);
                    dx(1)=dx(i);                                            %rest of above resets remaining variables to the ending previous values of the previous k loop
                    i=1;                %resets inner loop variable
  
                end

                %% Initalize inner loop variables

                t_cross=0;                                                  %number of times piston crosses zero
                p.t_force_adjust=0;                                         %adjustment in the time vector to account for variable frequency response
                x_piston_m=0;                                               %initialize x_piston
                F_drive=0;                                                  %initalize driving force vector
                F_drive_clean=0;
                %%%%%%%%%%
                %% Main Loop for calculations, i loop
                %%%%%%%%%%

                while size(t_cross)<p.n_period*2-1

                    %% Motor Model                  
                    if isempty(p.t_force_adjust) == 1
                        p.t_force_adjust=0;
                    end

                    if i == 1                                               %only want to perform these calculations on first time through
                        p.P_motor=p.eta_motor*p.P_electric;                 %power actually delivered to piston
                        p.Q_motor=p.P_electric-p.P_motor;                   %Heat transfer from motor

                        %Converts power input to force output by motor, based on motor literature
                        %p.k_m = -0.0008*p.P_motor^2 - 0.0047*p.P_motor + 4.6831;
                        p.k_m = 4.6831;

                        p.F_drive_max=sqrt(p.P_motor)*(p.k_m);              %Maximum Motor force, in Newtons.
                    end

                    F_drive_clean(i)=p.F_drive_max*(cos(p.w_d*(t(i)-p.t_force_adjust)+pi));                %Driving force v. Time, Newtons
                    
                    p.dP_piston = P(i) - P_cv2(i);
                    F_drive(i)=F_drive_clean(i) - p.dP_piston*p.Ap*1000;

                    %% Vibration Model and Volume and Piston Position  

                    p.dP_piston = P(i) - P_cv2(i);
                    
                    [ dV(i+1),V(i+1),x_piston(i+1),x_dot_piston(i+1),theta(i+1),theta_dot(i+1),x_piston_m(i+1),theta_dot_dot(i+1),p ] = vibration( t(i),dP(i),dx(i),x_piston(i),x_dot_piston(i),theta(i),theta_dot(i),x_piston_m(i),theta_dot_dot(i),p );
                    
                    c_eff(i) = p.c_eff;
                    c_friction(i) = p.c_friction;
                    c_gas(i) = p.c_gas;
                    x_dot_dot_piston(i) = p.x_dot_dot_piston;
                    k_eff(i) = p.k_eff;
                    F_gas(i) = (P(i) - P_cv2(i))*p.Ap;
                    F_wall(i) = p.F_wall;
                    p.dV = V(i+1) - V(i);
                    p.P_current = P(i);
                    
                    %% Volume for Second Control Volume
                    dV_cv2(i)=dV(i);
                    V_cv2(i)=p.V_cv2_i-V(i);

                    %% Massflow (Valve Model)
                    [ x_valve(i+1),x_dot_valve(i+1),dm_valve, Ma(i) ] = valve_dynamics(P(i),rho(i),T(i),x_valve(i),x_dot_valve(i),p);

                    %% Leakage Flowrates (Leakage Model)
                    [dm_leak_in(i),dm_leak_out(i),Ma_cv2(i)] = leakage( P(i),P_cv2(i),x_dot_piston(i),T(i),rho(i),T_cv2(i),x_piston(i),p );
                  

                    %% Total flowrates (Massflow Summations)
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

                    p.f_list=f_list(l);
                    [ Q(i) ] = Ins_HT( T(i),rho(i),T_w(k),V(i),dV(i),x_piston(i),x_dot_piston(i),p);

                    [ Q_cv2(i) ] = Ins_HT_cv2( T_cv2(i),rho_cv2(i),T_w(k),V_cv2(i),dV_cv2(i),x_dot_piston(i),p.Q_motor,p);

                    %% Mass of the gas in each control volume
                    m(i)=rho(i)*V(i);
                    m_cv2(i)=rho_cv2(i)*V_cv2(i);

                    %% Fluid Properties for approximation
                    if k==1                                                 %only need to compute this value once for each batch line
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

                        error('RK4 not available in this version');

                    elseif strcmp(p.method,'Euler')==1

                        %% Simple Euler Approximation

                        %%%%%%%%%%%%%%%%%
                        %Euler Approximation, Temperature and Density
                        %%%%%%%%%%%%%%%%%

                        %%%%%
                        %1st Control Volume, Piston Cylinder
                        %%%%%

                        dT(i)=(Q(i)+T(i)*dP_dT_v(i)*(1/1000)*((1/rho(i))*dm(i)-dV(i))-h(i)*dm(i)+dm_in(i)*p.h_in-dm_out(i)*h(i))/(m(i)*Cv(i));
                        drho(i)=(dm(i)-rho(i)*dV(i))/V(i);
                        
                        T(i+1)=T(i)+p.t_step*dT(i);
                        rho(i+1)=rho(i)+p.t_step*drho(i);

                        %%%%%%
                        %2nd Control Volume
                        %%%%%%
                        
                        dm_cv2=dm_leak_out(i)-dm_leak_in(i);
                        dm_in_cv2=dm_leak_out(i);
                        dm_out_cv2=dm_leak_in(i);

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
                    
                    %Section determines when the piston crosses back
                    %through the midpoint, signifying one cycle.
                    
                    if k>1
                        if i>1
                            t_cross=find_cross(t',x_piston_m(1:end-1)',0);
                        end
                    else
                        if i>1
                            t_cross=find_cross(t',x_piston_m(1:end-1)',0);
                        end
                    end

                    t(i+1)=t(i)+p.t_step;                                   %initially frequency has to float, transient frequency is different than driving frequency.
                    i=i+1;

                    if length(t)>=p.n*p.n_period*1.5
                        t_cross=ones(2*p.n_period);
                        disp('Convergence Error, piston did not cross zero twice')
                    end

                end     %End of While loop checking t_cross, i loop

                %%%%%%%%%%%%%%%
                %% Adjustments for Time
                %%%%%%%%%%%%%%%

                i=i-1;                                                      %step i back one for sanity
                t=t(1:end-1);                                               %step time vector back one
                temp(k)=max(t_cross)-max(find_cross(t',F_drive_clean',0));
                p.t_force_adjust=p.t_force_adjust+temp(k);                  %when t_force_adjust is greater than 3, stop i loop
                t_change(k)=p.t_force_adjust;
                p.Period=t(end)-t(1);                                       %floating period

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Estimate Massflow, Power Input,Heat Transfer, and Efficiencies
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
                %Mass flows calculated using numerical integration
                m_dot(k)=(trapz(dm_out_valve)*(t(2)-t(1)))/p.Period; 
                m_dot_leak_out(k)=(trapz(dm_leak_out)*(t(2)-t(1)))/p.Period; 
                m_dot_leak_in(k)=(trapz(dm_leak_in)*(t(2)-t(1)))/p.Period; 
                m_dot_in(k) = (trapz(dm_in)*(t(2) - t(1)))/p.Period;
                
                %W_dot calculated assuming max(h) is similar to h_avg
                %(proper calculation of h_avg done after the loop is over,
                %this value is just for diagnostic purposes)
                W_dot(k)=m_dot(k)*(max(h)-p.h_in); 
                
                %Efficiency estimations, diagnostic values
                eta_o(k)=(m_dot(k)*(p.h_2_s-p.h_in)*1000)/(p.P_electric*2);      %W_dot is in kW
                eta_vol(k)=(m_dot(k))/(p.rho_i*f_list(l)*(max(x_piston_m)-min(x_piston_m))*p.Ap);
                
                %Power consumed by friction
                W_dot_friction = (c_eff - c_gas)*p.x_max*f_list(l)^2;
                W_dot_friction_ave = sum(W_dot_friction)/(1000*length(W_dot_friction)); %converted to kW
                
                %Heat Transfer estimations, numerically integrated.
                Q_dot(k)=((trapz(Q)*(t(2)-t(1)))/p.Period)+W_dot_friction_ave;
                Q_dot_cv2(k)=(trapz(Q_cv2)*(t(2)-t(1)))/p.Period;

                T_w(k+1)=(Q_dot(k)+Q_dot_cv2(k))*1000*p.R_shell+p.T_amb;
                
                %Diagnostic values
                m_change(k) = m(end) - m(1);
                m_change_cv2(k) = m_cv2(end) - m_cv2(1);
                dE(k) = -Q_dot(k) + m_dot(k)*max(h) - (2*p.P_electric/1000) - m_dot_in(k)*p.h_in;
                resonant_freq(k) = sqrt(p.k_eff/p.M_mov)/2/pi;
 

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

                if k>1
                dx_res(k)=x_stroke_save(k)-x_stroke_save(k-1);
                end


                dT_loop(k)=dT_res(k)/min(T(1:length(t)));
                dP_loop=(P(1)-P(i))/min(P(1:length(t)));
                drho_loop(k)=(rho(1)-rho(i))/min(rho(1:length(t)));
                dT_loop_k=dT_loop(k);
                drho_loop_k=drho_loop(k);
                CV1_res=strcat('CV1 Res.:| ',' DT| ',num2str(dT_loop(k)),' Drho| ',num2str(drho_loop(k)),' dP| ',num2str(dP_loop));
                disp(CV1_res);

                if p.leakage_on==1
                    dT_res_cv2(k)=(T_cv2(1)-T_cv2(i));
                    drho_res_cv2(k)=(rho_cv2(1)-rho_cv2(i));


                    dT_cv2_loop(k)=(T_cv2(1)-T_cv2(i))/min(T_cv2(1:length(t)));
                    drho_cv2_loop(k)=(rho_cv2(1)-rho_cv2(i))/min(rho_cv2(1:length(t)));

                    dT_cv2_loop_k=dT_cv2_loop(k);
                    drho_cv2_loop_k=drho_cv2_loop(k);

                    dP_cv2_loop=(P_cv2(1)-P_cv2(i))/min(P_cv2(1:length(t)));
                    CV2_res=strcat('CV2 Res.:| ',' DT| ',num2str(dT_cv2_loop(k)),' Drho| ',num2str(drho_cv2_loop(k)),' dP| ',num2str(dP_cv2_loop));
                    disp(CV2_res);

                else
                    %dT_cv2_loop=0;
                    %drho_cv2_loop=0;
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
                    disp(Other_res)
                    add_info=strcat('Line: ',num2str(z),' | x_stroke: ',num2str(p.x_stroke),'| a: ',num2str(a),'| P_electric: ',num2str(p.P_electric));
                    disp(add_info)
                end

                %Plots during iterations
                if pc_flag == 0
                    figure(30)
                    plot(t,x_piston_m(1:end-1))
                    hold on
                    
                    W_dot_plot = W_dot_in(k)*ones(1,length(t));
                    figure(31)
                    plot(t,W_dot_plot)
                    hold on
                    
%                     figure(1)
%                     plot(k, dE(k),'bo-')
%                     hold on
%                     ylabel('Differential in Energy, (KW)')
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

                if k == 150
                    error = error + 0.00005;
                elseif k == 200
                    error = error + 0.00005;
                end

                %% Adjustments to power to try and maintain a correct stroke

                %If stroke is outside of the range of 0.001 and 0.0254 m either
                %add or substract 1W of power input to try and get it into the
                %appropriate range
                if k > 6
                    if p.x_stroke < 0.0001
                        p.P_electric = p.P_electric +1;
                        disp('Increase P_electric =')
                        disp(p.P_electric)
                    elseif p.x_stroke > 0.0254
                        %p.P_electric = p.P_electric -1;
                    end
                end

            end                                 %End of k while loop

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Brents Method: Attempting to find the correct power input to
            %% obtain full stroke
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [p,error_stroke,a] = brents(p,a,error_brents);
            if p.brent == 1
                P_electric_save(l) = p.P_brent(a);
                p.P_electric = p.P_brent(a);      
                
                if a == 50
                    error_brents= 0.000005; % 5 microns
                elseif a==60
                    error_brents= 0.00001; %10 microns
                elseif a==80
                    error_brents= 0.00005; %50 microns
                elseif a>100
                    error_brents= 1;
                end
                
            end
            
        end     %end a loop
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
            c_eff = c_eff(1:length(t));
            c_friction = c_friction(1:length(t));
            c_gas = c_gas(1:length(t));
            x_dot_dot_piston = x_dot_dot_piston(1:length(t));
            F_gas = F_gas(1:length(t));

            T_w = T_w(1:k-1);
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
            
            m_dot = m_dot(1:k-1);
            m_dot_leak_in = m_dot_leak_in(1:k-1);
            m_dot_leak_out = m_dot_leak_out(1:k-1);
            W_dot = W_dot(1:k-1);
            eta_o = eta_o(1:k-1);
            eta_vol = eta_vol(1:k-1);
            Q_dot = Q_dot(1:k-1);
            Q_dot_cv2 = Q_dot_cv2(1:k-1);
            m_dot_in = m_dot_in(1:k-1);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Estimate Massflow, Power Input,Heat Transfer, and Efficiencies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            % Estimating exit temperature and enthalpy
            if m_dot(k-1) > 0
                %finding key volumes
                V_1=max(V);                         %BDC
                V_3=min(V);                         %TDC
                V_4=min(find_cross(V',P',p.P_i));   %Suction valve open
                V_2=max(find_cross(V',P',p.P_d));   %discharge valve open

                %Find indexes of important volumes
                i_v1=find(V == V_1);
                i_v3=find(V == V_3);
                
                if i_v1 > i_v3
                    %Discharge peak is on the left
                    if max(crossing(V,[],V_2))&& min(crossing(V,[],V_2)) < i_v1
                        %i_v2 is a low value
                        i_v2=min(crossing(V,[],V_2));
                        i_v4=max(crossing(V,[],V_4));
                        T_d=mean(T(i_v2:i_v3));
                        h_2=mean(h(i_v2:i_v3));
                    else
                        %i_v2 is high
                        i_v2=max(crossing(V,[],V_2));
                        i_v4=max(crossing(V,[],V_4));
                        T_d=mean(cat(2,T(i_v2:end),T(1,i_v3)));
                        h_2=mean(cat(2,h(i_v2:end),h(1,i_v3)));
                    end
                else
                    %Discharge peak is on the right
                    if max(crossing(V,[],V_4))&& min(crossing(V,[],V_4)) < i_v3
                        %i_v4 is low
                        i_v2=max(crossing(V,[],V_2));
                        i_v4=min(crossing(V,[],V_4));
                        T_d=mean(T(i_v2:i_v3));
                        h_2=mean(h(i_v2:i_v3));                        
                    else
                        %i_v4 is high
                        i_v2=max(crossing(V,[],V_2));
                        i_v4=min(crossing(V,[],V_4));
                        T_d=mean(T(i_v2:i_v3));
                        h_2=mean(h(i_v2:i_v3));                        
                    end
                end
 
                %estimate exit temperature and enthalpy
                %T_d = sum(T(i_v2:i_v3))/(i_v3 - i_v2 + 1);
                %h_2 = sum(h(i_v2:i_v3))/(i_v3 - i_v2 + 1);
                h_2_calc = EOS(T_d,p.P_d,'enthalpy','P');
                
                %Boundary work calculations
                W_2_3 = mean(P(i_v2:i_v3))*(V_2 - V_3);
                W_4_1 = mean(P(i_v1:i_v4))*(V_1 - V_4);
                W_1_2 = mean(mean(P(i_v1:end)))*(V_1 - V_2);
                W_3_4 = mean(P(i_v3:i_v4))*(V_4 - V_3);
                %W_4_1 = sum(P(i_v4:end).*V(i_v4:end))/(length(P) - i_v4 + 1);
                %W_1_2 = sum(P(i_v1:i_v2).*V(i_v1:i_v2))/(i_v2 - i_v1 + 1);
                %W_3_4 = sum(P(i_v3:i_v4).*V(i_v3:i_v4))/(i_v4 - i_v3 + 1);
                
                %Actual boundary work of real compression process
                W_boundary = W_1_2 + W_2_3 - W_3_4 - W_4_1;
                W_dot_boundary = W_boundary*f_list(l);
                
                %Idealized suction and discharge processes
                W_gas_dis_ideal = p.P_d*(V_2 - V_3);
                W_gas_suc_ideal = p.P_i*(V_1 - V_4);

                %Hybrid idealized boundary work (compression and expansion
                %are real)
                W_gas_net_ideal = W_gas_dis_ideal + W_1_2 - W_3_4 - W_gas_suc_ideal;
                W_dot_net_ideal = W_gas_net_ideal*f_list(l);
                
                %Valve Losses, based on differences between ideal valves and
                %real valves, in W
                W_dot_suc_loss = 1000*(-W_4_1 + W_gas_suc_ideal)*f_list(l);
                W_dot_dis_loss = 1000*(W_2_3 - W_gas_dis_ideal)*f_list(l);

            else

                T_d = 0;
                h_2 = 0;
                h_2_calc = 0;
                W_dot_suc_loss = 0;
                W_dot_dis_loss = 0;


            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Losses and power calculations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Power in the piston shaft, energy integrated over one period
            W_dot_shaft = 0.5*p.M_mov*x_dot_piston.^2*f_list(l);
            W_dot_shaft_ave = sum(W_dot_shaft)/length(W_dot_shaft);
            
            %Power in the rotation of piston shaft, energy integrated over
            %one period
            W_dot_shaft_rotation = 0.5*p.J*theta_dot.^2*f_list(l);
            W_dot_shaft_rotation_ave = sum(W_dot_shaft_rotation)/length(W_dot_shaft_rotation);
            
            %Power consumed by friction
            W_dot_friction = (c_eff - c_gas)*p.x_max*f_list(l)^2;
            W_dot_friction_ave = sum(W_dot_friction)/length(W_dot_friction);
            
            %Total power consumed by compressor, moving from compression
            %chamber to motor
            W_dot_total_ave = W_dot_shaft_ave + W_dot_shaft_rotation_ave + W_dot_friction_ave + p.Q_motor;
       
            %Maximum power stored in the mechanical springs
            W_dot_stored_max = 0.5*p.k_mech*p.x_stroke^2*f_list(l);

            %Net compressor power into the gas
            W_dot(k-1)=1000*m_dot(k-1)*(h_2-p.h_in); %W_dot is in W
            
            %Efficiency estimations
            eta_o(k-1)=(m_dot(k-1)*(p.h_2_s-p.h_in)*1000)/(p.P_electric*2);
            eta_o_test = (m_dot(k-1)*(p.h_2_s-p.h_in)*1000)/W_dot_total_ave;
            eta_is = (p.h_2_s - p.h_in)/(h_2 - p.h_in);

            %Net leakage losses
            W_dot_leakage_loss = W_dot_shaft_ave - W_dot(k-1) - Q_dot(k-1) - W_dot_suc_loss - W_dot_dis_loss;
            W_dot_leak = 1000*mean(m_dot_leak_out(end)*h - m_dot_leak_in(end)*h_cv2);
            
            %% Other Estimations
            
            x_mag(l)=p.x_stroke/p.x_max; %#ok<SAGROW>
            x_stroke_resonant_save(l) = (max(x_piston_m)-min(x_piston_m)); %#ok<SAGROW>
            f_list_resonant_save(l) = f_list(l); %#ok<SAGROW>
  
            
            W_gas=((p.gamma*p.P_i*(max(V) - p.V_dead_valve))/(p.gamma-1))*(((p.P_d/p.P_i)^((p.gamma-1)/p.gamma))-1); %kJ
            W_dot_gas = W_gas*f_list(l)*1000;
            
            
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
                m_dot(k-1),m_dot_leak_out(k-1),p.P_i,p.P_d,p.T_i,length(t),char(p.method),2*p.P_electric,...
                p.eta_motor,p.k_mech,p.f_friction,p.ecc_1,x_mag(l),p.g,eta_o(k-1),eta_vol(k-1),T_d,...
                T_w(k-1),k-1,Q_dot(k-1),Q_dot_cv2(k-1),p.P_d/p.P_i,dT_loop(k-1),drho_loop(k-1),dT_cv2_loop(k-1),drho_cv2_loop(k-1),...
                dx_piston,dm_dot_cycle,dT_w,p.M_mov,(max(x_piston_m)-min(x_piston_m)),p.P_electric,f_resonant_calculated,w_d_resonant_calculated,p.k_eff,eta_is,...
                eta_o_test,W_dot_shaft_ave,W_dot_shaft_rotation_ave,W_dot_friction_ave,W_dot_total_ave,W_dot_stored_max,W_dot_suc_loss,W_dot_dis_loss,...
                W_dot_leakage_loss,W_dot_leak,mean(x_piston_m)};

            range_2=strcat('A',num2str(l+rows),':','BB',num2str(l+rows));
            xlswrite1('data_log.xls',save_data,range_2)

            save_name = strcat('wksp_',p.save_name);
            save_name = char(save_name);  
            save(save_name)

        elseif pc_flag == 1

            save_data_tmp = {z,p.w_d,f_list(l),(max(x_piston_m)-min(x_piston_m)),W_dot(k-1),...
                m_dot(k-1),m_dot_leak_out(k-1),p.P_i,p.P_d,p.T_i,length(t),2*p.P_electric,...
                p.eta_motor,p.k_mech,p.f_friction,p.ecc_1,x_mag(l),p.g,eta_o(k-1),eta_vol(k-1),T_d,...
                T_w(k-1),k-1,Q_dot(k-1),Q_dot_cv2(k-1),p.P_d/p.P_i,dT_loop(k-1),drho_loop(k-1),dT_cv2_loop(k-1),drho_cv2_loop(k-1),...
                dx_piston,dm_dot_cycle,dT_w,p.M_mov,(max(x_piston_m)-min(x_piston_m)),p.P_electric,f_resonant_calculated,w_d_resonant_calculated,p.k_eff,eta_is,...
                eta_o_test,W_dot_shaft_ave,W_dot_shaft_rotation_ave,W_dot_friction_ave,W_dot_total_ave,W_dot_stored_max,W_dot_suc_loss,W_dot_dis_loss,...
                W_dot_leakage_loss,W_dot_leak,mean(x_piston_m)};


            save_data(z,:) = cell2mat(save_data_tmp);

            save_name = strcat('data_log_',p.save_name);
            wksp_name = strcat('wksp_',p.save_name);
            save_name = char(save_name);
            xlswrite(save_name,save_data);
            save(char(wksp_name))

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


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % P-V Plot
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            figure(21)
            %set(gcf,'Visible','off')
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
            %set(gcf,'Visible','off')
            subplot(3,3,1),plot(t,T);title('Temperature v. Time');
            subplot(3,3,2),plot(t,P);title('Pressure v. Time');
            subplot(3,3,3),plot(t,dm);title('Massflow v. Time');
            subplot(3,3,4),plot(t,V);title('Volume v. Time');
            subplot(3,3,5),plot(t,V);title('Volume v. Time');
            subplot(3,3,6),plot(t,rho);title('Density v. Time');
            subplot(3,3,7),plot(t,dV);title('Change in Volume v. Time');
            subplot(3,3,8),plot(t,x_piston);title('Displacement v. Time');
            subplot(3,3,9),plot(t,drho);title('Change in Rho v. Time');

            figure(26)
            %set(gcf,'Visible','off')
            %plot(t,dT_check);hold on;plot(t,dT,'r');legend('dT_c_h_e_c_k','dm')
            subplot(2,2,1),plot(t,Ma);title('Mach Number v. Time')
            subplot(2,2,2),plot(t,x_valve);title('Valve Lift v. Time')
            subplot(2,2,3),plot(m_dot_leak_out);title('Leakage Outflow v. Time')
            subplot(2,2,4),plot(m_dot);title('Massflow v. Iteration')

            % CV2 Plots
            figure(27)
            %set(gcf,'Visible','off')
            subplot(3,3,1),plot(t,T_cv2);title('Temperature v. Time for CV2');
            subplot(3,3,2),plot(t,P_cv2);title('Pressure v. Time for CV2');
            subplot(3,3,3),plot(t,dm_leak_in);title('Leakage Massflow In v. Time for CV2');
            subplot(3,3,4),plot(t,V_cv2);title('Volume v. Time for CV2');
            subplot(3,3,5),plot(t,rho_cv2);title('Density v. Time');
            subplot(3,3,6),plot(t,h_cv2);title('Enthalpy v. Time');
            subplot(3,3,7),plot(t,dT_cv2);title('Change in Temperature v. Time for CV2');
            subplot(3,3,8),plot(t,Ma_cv2);title('Mach Number v. Time');
            subplot(3,3,9),plot(t,drho_cv2);title('Change in Rho v. Time for CV2');


                %% Shutting down activex server to excel
                invoke(Excel.ActiveWorkbook,'Save');

                Excel.Quit
                Excel.delete
                clear Excel   

        end

        %% other stuff

        disp(strcat('Batch Line No: ',num2str(z)))

    end         %end of z loop (number of batch lines)

    %turn annoying warning back on
    warning('on','MATLAB:Print:SavingToDifferentName')

    disp('finished')
    
    diary

catch ME

    disp('finished w/ error')
    disp(ME.message)
    disp(ME.cause)
    disp(ME.stack)
    disp(ME)
    diary
    
    if pc_flag == 0
        %% Shutting down activex server to excel
        invoke(Excel.ActiveWorkbook,'Save');

        Excel.Quit
        Excel.delete
        clear Excel   
        
        crash_name = strcat('crash_',save_name);
        save(crash_name)
        
    end
    
    if pc_flag == 1

        save_name = strcat('data_log_',p.save_name);
        save_name = char(save_name);
        xlswrite(save_name,save_data);
        crash_name = strcat('crash_',save_name);
        save(crash_name)

    end

%turn annoying warning back on
warning('on','MATLAB:Print:SavingToDifferentName')

end














































