%% Clean up

clc
clear all


%% 1 - Add header to data_log

log_file = 'data_log_ecc_g_36.csv';
combined_file = 'data_log_ecc_g.xlsx';

[num,txt]=xlsread(log_file);
[num2,txt2]=xlsread('header.xlsx');

if exist(combined_file,'file') == 0
    
    xlswrite(combined_file,txt2)
    xlswrite(combined_file,num,'sheet1','A3')
    
end

%% Diagnostic Plots
close all

colors=char('b','k','r','bo-','ko-','ro-');
g=[18,15];    %Sweeping leakage gap (18), and mu (15)

for z=1:1:length(g)

    %Figure 1
    figure
    var=[12 15 6 19 20 7 45 49 51];
    % 12 - P_electric
    % 15 - mu
    % 6  - m_dot
    % 19 - eta_o
    % 20 - eta_vol
    % 7  - m_dot_leakage
    % 45 - W_dot_total
    % 49 - W_dot_leakage_loss
    % 51 - Equilibrium Position

    for i=1:1:length(var)
        vars=var(i); %y is m_dot
        n=1;        %initial starting point.
        subplot(3,3,i)
        for j=1:1:6
            plot(num(n:n+5,g(z)),num(n:n+5,vars),colors(j,:),'LineWidth',2)
            hold on
            n=n+6;
        end
        xlabel(txt2(1,g(z)))
        ylabel(txt2(1,vars))
    end

    %Figure 2
    figure
    var=[41 40 21 24 25 44];
    % 41 - eta_o_test
    % 40 - eta_is
    % 21 - T_d
    % 24 - Q_dot
    % 25 - Q_dot_cv2
    % 44 - W_dot_friction

    for i=1:1:length(var)
        vars=var(i); %y is m_dot
        n=1;        %initial starting point.
        subplot(3,3,i)
        for j=1:1:6
            plot(num(n:n+5,g(z)),num(n:n+5,vars),colors(j,:),'LineWidth',2)
            hold on
            n=n+6;
        end
        xlabel(txt2(1,g(z)))
        ylabel(txt2(1,vars))
    end

end

