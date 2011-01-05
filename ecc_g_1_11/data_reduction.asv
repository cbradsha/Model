%% Clean up

clc
clear all


%% 1 - Add header to data_log

log_file = 'data_log_ecc_g_36.csv';
combined_file = 'data_log_ecc_g.xlsx';

[num,txt]=xlsread(log_file);
[num2,txt2]=xlsread('header.xlsx');

if length(txt) == 0
    
    xlswrite(combined_file,txt2)
    xlswrite(combined_file,num,'sheet1','A3')
    
end

%% Diagnostic Plots
close all

colors=char('b','k','r','bo-','ko-','ro-');
g=18;    %Sweeping leakage gap

figure

% Power input vs. leakage gap
P_electric=12; %y is P_electric
n=1;        %initial starting point.
subplot(3,3,1)
for j=1:1:6
    plot(num(n:n+5,g),num(n:n+5,P_electric),colors(j,:),'LineWidth',2)
    hold on
    n=n+6;
end
xlabel(txt2(1,g))
ylabel(txt2(1,P_electric))


% Friction factor vs. leakage gap
mu=15; %y is mu
n=1;        %initial starting point.
subplot(3,3,2)
for j=1:1:6
    plot(num(n:n+5,g),num(n:n+5,mu),colors(j,:),'LineWidth',2)
    hold on
    n=n+6;
end
xlabel(txt2(1,g))
ylabel(txt2(1,mu))


% Mass Flow vs. leakage gap
m_dot=6; %y is m_dot
n=1;        %initial starting point.
subplot(3,3,3)
for j=1:1:6
    plot(num(n:n+5,g),num(n:n+5,m_dot),colors(j,:),'LineWidth',2)
    hold on
    n=n+6;
end
xlabel(txt2(1,g))
ylabel(txt2(1,m_dot))


% Leakage Mass Flow vs. leakage gap
m_dot_leak=7; %y is m_dot
n=1;        %initial starting point.
subplot(3,3,6)
for j=1:1:6
    plot(num(n:n+5,g),num(n:n+5,m_dot_leak),colors(j,:),'LineWidth',2)
    hold on
    n=n+6;
end
xlabel(txt2(1,g))
ylabel(txt2(1,m_dot_leak))

% eta_o vs. leakage gap
eta_o=19; %y is m_dot
n=1;        %initial starting point.
subplot(3,3,4)
for j=1:1:6
    plot(num(n:n+5,g),num(n:n+5,eta_o),colors(j,:),'LineWidth',2)
    hold on
    n=n+6;
end
xlabel(txt2(1,g))
ylabel(txt2(1,eta_o))

% eta_vol vs. leakage gap
eta_vol=20; %y is m_dot
n=1;        %initial starting point.
subplot(3,3,5)
for j=1:1:6
    plot(num(n:n+5,g),num(n:n+5,eta_vol),colors(j,:),'LineWidth',2)
    hold on
    n=n+6;
end
xlabel(txt2(1,g))
ylabel(txt2(1,eta_vol))


% W_dot_total vs. leakage gap
W_dot_total=45; %y is m_dot
n=1;        %initial starting point.
subplot(3,3,7)
for j=1:1:6
    plot(num(n:n+5,g),num(n:n+5,W_dot_total),colors(j,:),'LineWidth',2)
    hold on
    n=n+6;
end
xlabel(txt2(1,g))
ylabel(txt2(1,W_dot_total))

% W_dot_leakage vs. leakage gap
W_dot_leakage=49; %y is m_dot
n=1;        %initial starting point.
subplot(3,3,8)
for j=1:1:6
    plot(num(n:n+5,g),num(n:n+5,W_dot_leakage),colors(j,:),'LineWidth',2)
    hold on
    n=n+6;
end
xlabel(txt2(1,g))
ylabel(txt2(1,W_dot_leakage))


% Equilibrium Position vs. leakage gap
eq=51; %y is m_dot
n=1;        %initial starting point.
subplot(3,3,9)
for j=1:1:6
    plot(num(n:n+5,g),num(n:n+5,eq),colors(j,:),'LineWidth',2)
    hold on
    n=n+6;
end
xlabel(txt2(1,g))
ylabel(txt2(1,eq))
