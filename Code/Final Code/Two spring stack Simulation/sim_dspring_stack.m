%% Simulation of 2-spring stack - Mass system in MATLAB

% ****** RUN THIS IF YOU WANT TO SIMULATE FOR A PARTICULAR SET OF SPRING
% AND MASS PARAMETERS. CHANGE DAMPING PARAMETERS IN 'dyn_dspring_stack.m' ******

% The code below is used to simulate the dynamics of a 2-spring stack
% consisting of disc-springs with a mass attached to it. 
% The displacement and velocity profiles of the mass are plotted 
% and an animation of the spring-mass motion is also shown.

% Convention - Positive values of displacement and velocity are in the
% "downward direction" for convenience (since the springs are mostly in 
% compression).

% In the Problem Formulation Report, the positive values of displacement
% and velocity are in the upward direction. In order to use that, one needs
% to reverse the signs of x_st, x_1, delta_st, delta_1 and their derivatives 
% obtained in the simulation.



%% Clear and Close All windows

clear;
close all;

%% System Parameters
m1 = 11.2; % Mass (kg)
g = 9.81; % Acceleration due to gravity (m/s^2)
ht1_ratio = 1.41; % Ratio of height to thickness of spring 1 (the top spring)
tau1 = 0.5; % Thickness of spring 1 (mm)
ht2_ratio = 1.41; % Ratio of height to thickness of spring 2 (the bottom spring)
tau2 = 0.5; % Thickness of spring 2 (mm)
h01 = ht1_ratio * tau1; % Height of spring 1 (mm)
h02 = ht2_ratio * tau2; % Height of spring 2 (mm)
l0 = 4*(h01 + h02); % Natural length of the spring + spacers (mm)
l1 = 4*h02 + h01; % Distance between the base and the top of spring 2 (mm)

% Temporal parameters
dt = 1e-3; % Time step (s)
Tmax =  15; % Maximum time (s)
tspan = 0:dt:Tmax; % Duration of the simulation (s)

% Base excitation parameters
A = 0.5; % Peak to peak amplitude in mm
f_no = 4; % Frequency sweep index (1 corresponds to 4 Hz, 2 corresponds to 6 Hz)
f_excit = 4 + 2 * (f_no - 1); % Base excitation frequency in Hz

%% Solve for static equilibrium

% Obtain the displacement of each spring in static equilibrium
% x_eq(1) and x_eq(2) are the displacements of springs 1 and 2 respectively

x_eq = static_eqb(ht1_ratio, tau1, ht2_ratio, tau2, m1*g);


%% Set the initial state

x_base_init = 0; % Initial displacement of the base (mm)
U_init = 0;
x0 = [x_eq(1);0;x_base_init;U_init]; % Initial state vector (Static equilibrium state)

%% Solve the ODEs to obtain the displacements and velocities numerically 

% Using ode89
options=odeset('abstol',1e-9,'reltol',1e-9); % Tolerances
[t,x] = ode89(@(t,x)dyn_dspring_stack(t,x,m1,g,ht1_ratio,tau1,ht2_ratio,tau2,A,f_excit),tspan,x0,options);

% Calculate displacement of spring 2 (y)
y = zeros(length(t),1); % Placeholder array to store the displacement of spring 2
for i=1:length(t)
    fun = @(y)equality(y, x(i,1), x(i,3), ht1_ratio, tau1, ht2_ratio, tau2); % Based on equality of forces in both springs
    y(i,1) = fsolve(fun, 0); % Calculate spring 2 displacement  
end

% Calculate forces in springs 1 and 2
fs1 = zeros(length(t),1); % Placeholder array for forces in spring 1
fs2 = zeros(length(t),1); % Placeholder array for forces in spring 2
for i=1:length(t)
    % Calculate based on non-linear disc-spring model
    fs1(i,1) = disc_spring_force(x(i,1) - y(i,1),ht1_ratio,tau1); 
    fs2(i,1) = disc_spring_force(y(i,1) - x(i,3),ht2_ratio,tau2);  
end

% Obtain acceleration of the mass

dx_arr = zeros(length(t),4); % Placeholder array for the state evolution
% Obtain the state evolution vector for each time instant
for i=1:length(t)
    dx_arr(i,:) = dyn_dspring_stack(t(i),x(i,:).',m1,g,ht1_ratio,tau1,ht2_ratio,tau2,A,f_excit).';
end

acc = dx_arr(:,2); % Store the acceleration of the mass

%% Interface for the plots and animation

% Window position and xtick settings
set(gcf,'Position',[0,0,1500,1000],'defaultAxesTickLabelInterpreter','latex'); 

% Displacement plot
subplot(3,2,1)
plot(tspan,x(:,1),'LineWidth',2,'Color','Blue') % Displacement of the mass
hold on
plot(tspan,y(:,1),'LineWidth',2,'Color','Red')
hold on
plot(tspan,x(:,3),'LineWidth',2,'Color','Black') % Base displacement
title("Displacement profile","Interpreter","latex","FontSize",15)
legend("$x_{st}$","$x_1$","$x_{base}$","Interpreter","latex","FontSize",15)
%xlim([0 0.5])
ylim([-0.5 3.5])
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$x$ [mm]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

% Velocity plot
subplot(3,2,3)
plot(tspan,x(:,2),'LineWidth',2,'Color','Blue')
title("Velocity profile","Interpreter","latex","FontSize",15)
legend("$\dot{x}_{st}$","Interpreter","latex","FontSize",15)
%xlim([0 0.5])
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$\dot{x}$ [mm/s]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

% Acceleration plot
subplot(3,2,5)
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(t,acc,'LineWidth',2,'Color','Blue')
title("Acceleration profile","Interpreter","latex","FontSize",15)
legend("$\ddot{x}_{st}$","Interpreter","latex","FontSize",15)
%xlim([0 0.5])
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$\ddot{x}$ [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)


% Animation
for i=1:10:length(tspan)
    subplot(3,2,[2 6])
    plot([-17.25,-17.25],[x_base_init - x(i,3), 2*h02 - x(i,3)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([17.25,17.25],[x_base_init - x(i,3), 2*h02 - x(i,3)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([-17.25,-17.25],[l0 - x(i,1) - 2*h01, l0 - x(i,1)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([17.25,17.25],[l0 - x(i,1) - 2*h01, l0 - x(i,1)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([-11.2,-11.2],[l1 - y(i,1) - (h01 + h02), l1 - y(i,1)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([11.2,11.2],[l1 - y(i,1) - (h01 + h02), l1 - y(i,1)],'LineWidth',10,'Color','#5A5A5A')
    hold on

    plot([-17.25,-11.2],[2*h02 - x(i,3), l1 - y(i,1) - (h01 + h02)],'LineWidth',10,'Color','Blue')
    hold on
    plot([-11.2,-17.25],[l1 - y(i,1), l0 - x(i,1) - 2*h01],'LineWidth',10,'Color','Blue')
    hold on
    plot([17.25,11.2],[2*h02 - x(i,3), l1 - y(i,1) - (h01 + h02)],'LineWidth',10,'Color','Blue')
    hold on
    plot([11.2,17.25],[l1 - y(i,1), l0 - x(i,1) - 2*h01],'LineWidth',10,'Color','Blue')
    hold on


    rectangle('Position',[-21 l0-x(i,1) 42 1],'LineWidth',3,'FaceColor','#853A00')
    hold on
    rectangle('Position',[-23 x_base_init-x(i,3)-0.3 46 0.3],'LineWidth',3,'FaceColor','#808090')
    hold off
    axis([-35,35,-2,15])
    title("Motion of the Spring","Interpreter","latex","FontSize",15)
    xlabel('Lateral location [mm]',"Interpreter","latex","FontSize",15)
    ylabel('Longitudinal location [mm]',"Interpreter","latex","FontSize",15)
    set(gca, 'FontSize',15)
    pause (0.01)
end


%% Phase portraits

% Velocity vs Displacement
figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(x(:,1),x(:,2),'LineWidth',2,'Color','Blue')
title("Phase portrait (Velocity vs Displacement)","Interpreter","latex","FontSize",15)
xlabel('$x$ [mm]',"Interpreter","latex","FontSize",15)
ylabel('$\dot{x}$ [mm/s]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

% Acceleration vs Velocity
figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(x(:,2),acc,'LineWidth',2,'Color','Blue')
title("Phase portrait (Acceleration vs Velocity)","Interpreter","latex","FontSize",15)
xlabel('$\dot{x}$ [mm]',"Interpreter","latex","FontSize",15)
ylabel('$\ddot{x}$ [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

%% Spring Force(s) vs Time

% Compare forces in springs 1 and 2
figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(t,fs1,'LineWidth',2,'Color','Blue')
hold on
plot(t,fs2,'LineWidth',2,'Color','Red')
title("Force vs Time","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$F(t)$ [N]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

% Visualize the difference in forces present in springs 1 and 2
figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(t,fs1 - fs2,'LineWidth',2,'Color','Black')
title("Difference in Force vs Time","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$F(t)$ [N]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)



%% Force (F) vs Displacement of the mass (x) - 

%{
figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(x(:,1), fs1, 'LineWidth', 2, 'Color', 'Blue')
title( "Force vs Displacement of the mass","Interpreter","latex","FontSize",15)
xlabel('$x$ [mm]',"Interpreter","latex","FontSize",15)
ylabel('$F$ [N]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)
%}

figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
x_rel = x(:,1) - x(:,3); % Relative displacement between the mass and and the base
plot(x_rel, fs1, 'LineWidth', 2, 'Color', 'Blue')
title( "Force vs Displacement of the mass","Interpreter","latex","FontSize",15)
xlabel('$x$ [mm]',"Interpreter","latex","FontSize",15)
ylabel('$F$ [N]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

%% Stiffness (k) vs Time (t)

% Is this the best way of numerical differentiation ?
dP = diff(fs1);
d_delta = diff(x(:,1));

k_arr = dP./d_delta; % Local stiffness (k) = dP/d_delta

figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(tspan(1:end-1), k_arr, 'LineWidth', 2, 'Color', 'Blue')
title( "Stiffness vs Time","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$k$ [N/mm]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

%% Energy(En) vs Time(t)

% Useful as a check. Energy should be conserved !

KE_arr = zeros(length(t),1); % Array of Kinetic Energy (KE) Values
PE_arr = zeros(length(t),1); % Array of Potential Energy (PE) Values
TE_arr = zeros(length(t),1); % Array of Total Energy (TE) Values

for i=1:length(t)
    KE_arr(i) = 0.5 * m1 * (x(i,2)/1000)^2; % KE at time instant i
    PE_arr(i) = x(i,4); % PE at time instant i
    TE_arr(i) = KE_arr(i) + PE_arr(i); % TE at time instant i
end

figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(t, KE_arr, 'LineWidth', 2, 'Color', 'Blue')
hold on
plot(t, PE_arr, 'LineWidth', 2, 'Color', 'Red')
hold on
plot(t, TE_arr, 'LineWidth', 2, 'Color', 'Black')
legend("Kinetic","Potential","Total","Interpreter","latex","FontSize",15)
title( "Energy vs Time","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('Energy [J]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)


%% FFT Analysis of Acceleration

% Performing FFT analysis on the acceleration data in steady-state
acc_ss = acc(floor(end/2):end,1); % Acceleration values in steady-state
fft_analysis(tspan(floor(end/2):end),acc_ss.');


% Extra frequencies apart from the base excitation frequency visible possibly due to
% non-linear behaviour of the system

%% Calculate motion transmissibility (MT)

acc_base = -(2 * pi * f_excit)^2 * x(:,3);
acc_base_ptp = max(acc_base) - min(acc_base);

MT = (max(acc_ss) - min(acc_ss))/acc_base_ptp;
dB = 20*log10(MT);
disp("Motion Transmissibility = " + string(dB));

%% Try
%{
T = tiledlayout(2,2);
set(gcf,'Position',[0,0,800,500],'defaultAxesTickLabelInterpreter','latex'); 

nexttile
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(t,x(:,1),'LineWidth',2,'Color','Blue')
title("Displacement profile","Interpreter","latex","FontSize",15)
legend("$x_{st}$","Interpreter","latex","FontSize",15)
xlim([min(t), max(t)])
ylim([min(x(:,1))-0.1*abs(min(x(:,1))-mean(x(:,1))),max(x(:,1))+ 0.1*abs(max(x(:,1))-mean(x(:,1)))])
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$x$ [mm]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

nexttile
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(t,x(:,2),'LineWidth',2,'Color','Blue')
title("Velocity profile","Interpreter","latex","FontSize",15)
legend("$\dot{x}_{st}$","Interpreter","latex","FontSize",15)
xlim([min(t), max(t)])
ylim([min(x(:,2))-0.1*abs(min(x(:,2))),max(x(:,2))+ 0.1*abs(max(x(:,2)))])
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$\dot{x}$ [mm/s]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

nexttile
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(t,acc,'LineWidth',2,'Color','Blue')
title("Acceleration profile","Interpreter","latex","FontSize",15)
legend("$\ddot{x}_{st}$","Interpreter","latex","FontSize",15)
xlim([min(t), max(t)])
ylim([min(acc)-0.1*abs(min(acc)),max(acc)+ 0.1*abs(max(acc))])
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$\ddot{x}$ [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

% Velocity vs Displacement
nexttile;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(x(:,1),x(:,2),'LineWidth',2,'Color','Blue')
title("Phase portrait (Velocity vs Displacement)","Interpreter","latex","FontSize",15)
xlim([min(x(:,1))-0.1*abs(min(x(:,1))-mean(x(:,1))),max(x(:,1))+ 0.1*abs(max(x(:,1))-mean(x(:,1)))])
ylim([min(x(:,2))-0.1*abs(min(x(:,2))),max(x(:,2))+ 0.1*abs(max(x(:,2)))])
xlabel('$x$ [mm]',"Interpreter","latex","FontSize",15)
ylabel('$\dot{x}$ [mm/s]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

exportgraphics(T,'plots.pdf','Resolution',300)
%}

num_cyc = 6;
start_t = max(t) - (num_cyc/f_excit);
end_t = max(t);
min_y = min(min(acc),min(acc_base));
max_y = max(max(acc),max(acc_base));

T = tiledlayout(1,1);
set(gcf,'Position',[0,0,800,500],'defaultAxesTickLabelInterpreter','latex'); 

nexttile
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 

plot(t - start_t,acc,'LineWidth',2,'Color','Blue')
hold on
plot(t - start_t, acc_base,'LineWidth',2,'Color','Red')
title("Acceleration profile","Interpreter","latex","FontSize",15)
legend("$\ddot{x}_{st}$","$\ddot{x}_{base}$","Interpreter","latex","FontSize",15)
xlim([0, end_t - start_t])
ylim([min_y-0.1*abs(min_y),max_y+ 0.1*abs(max_y)])
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$\ddot{x}$ [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

plotpath = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Dynamic/Sim3', 'acc_' + string(f_excit) + '_sim.pdf');

%{
exportgraphics(T,plotpath,'Resolution',300)
%}

%% Superimpose the dynamic behaviour on the static force-deflection curve

% superimpose