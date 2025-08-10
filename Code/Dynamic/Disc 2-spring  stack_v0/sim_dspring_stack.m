%% Simulation of 2-spring stack - Mass system in MATLAB

% System Parameters
m1 = 10.7; % Mass (kg)
m2 = 1e-6; % Virtual mass (kg)
b1 = 0; % Damping coefficient for spring 1
b2 = 0; % Damping coefficient for spring 2
g = 9.81; % Acceleration due to gravity (m/s^2)
ht1_ratio = 1.6; % Ratio of height to thickness of spring 1
tau1 = 0.5; % Thickness of spring 1 (mm)
ht2_ratio = 1.6; % Ratio of height to thickness of spring 2
tau2 = 0.5; % Thickness of spring 2 (mm)
h01 = ht1_ratio * tau1; % Height of spring 1 (mm)
h02 = ht2_ratio * tau2; % Height of spring 2 (mm)
l0 = 4*(h01 + h02); % Natural length of the spring + spacers (mm)
l1 = 4*h02 + h01;

% Temporal parameters
dt = 0.001; % Time step (s)
Tmax =  20; % Maximum time (s)
tspan = 0:dt:Tmax; % Duration of the simulation (s)

x_base_init = 0; % Initial displacement of the base (mm)
x0 = [0;0;0;0;x_base_init]; % Initial state vector

%% Solve the ODEs to obtain the states numerically 

options=odeset('abstol',1e-9,'reltol',1e-9);
[t,x] = ode45(@(t,x)dyn_dspring_stack(t,x,m1,m2,g,b1,b2,ht1_ratio,tau1,ht2_ratio,tau2),tspan,x0,options);


%% Interface for the plots and animation

set(gcf,'Position',[0,0,1500,1000],'defaultAxesTickLabelInterpreter','latex'); 

% Displacement plot
subplot(2,2,1)
plot(tspan,x(:,1),'LineWidth',2,'Color','Blue')
hold on
plot(tspan,x(:,3),'LineWidth',2,'Color','Red')
title("Displacement profile","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$x$ [mm]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

% Velocity plot
subplot(2,2,3)
plot(tspan,x(:,2),'LineWidth',2,'Color','Blue')
hold on
plot(tspan,x(:,4),'LineWidth',2,'Color','Red')
title("Velocity profile","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$\dot{x}$ [mm/s]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

% Animation
for i=1:100:length(tspan)
    subplot(2,2,[2 4])
    plot([-17.25,-17.25],[x_base_init - x(i,5), 2*h02 - x(i,5)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([17.25,17.25],[x_base_init - x(i,5), 2*h02 - x(i,5)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([-17.25,-17.25],[l0 - x(i,1) - 2*h01, l0 - x(i,1)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([17.25,17.25],[l0 - x(i,1) - 2*h01, l0 - x(i,1)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([-11.2,-11.2],[l1 - x(i,3) - (h01 + h02), l1 - x(i,3)],'LineWidth',10,'Color','#5A5A5A')
    hold on
    plot([11.2,11.2],[l1 - x(i,3) - (h01 + h02), l1 - x(i,3)],'LineWidth',10,'Color','#5A5A5A')
    hold on

    plot([-17.25,-11.2],[2*h02 - x(i,5), l1 - x(i,3) - (h01 + h02)],'LineWidth',10,'Color','Blue')
    hold on
    plot([-11.2,-17.25],[l1 - x(i,3), l0 - x(i,1) - 2*h01],'LineWidth',10,'Color','Blue')
    hold on
    plot([17.25,11.2],[2*h02 - x(i,5), l1 - x(i,3) - (h01 + h02)],'LineWidth',10,'Color','Blue')
    hold on
    plot([11.2,17.25],[l1 - x(i,3), l0 - x(i,1) - 2*h01],'LineWidth',10,'Color','Blue')
    hold on


    rectangle('Position',[-21 l0-x(i,1) 42 1],'LineWidth',3,'FaceColor','#853A00')
    hold on
    rectangle('Position',[-23 x_base_init-x(i,5)-0.3 46 0.3],'LineWidth',3,'FaceColor','#808090')
    hold off
    axis([-35,35,-2,15])
    title("Motion of the Spring","Interpreter","latex","FontSize",15)
    xlabel('Lateral location [mm]',"Interpreter","latex","FontSize",15)
    ylabel('Longitudinal location [mm]',"Interpreter","latex","FontSize",15)
    set(gca, 'FontSize',15)
    pause (0.01)
end

%% Phase portrait

figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(x(:,1),x(:,2),'LineWidth',2,'Color','Blue')
title("Phase portrait","Interpreter","latex","FontSize",15)
xlabel('$x$ [mm]',"Interpreter","latex","FontSize",15)
ylabel('$\dot{x}$ [mm/s]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

fs1 = zeros(length(t),1);
fs2 = zeros(length(t),1);

%% Spring force vs Time

for i=1:length(t)
    fs1(i,1) = disc_spring_force(x(i,1) - x(i,3),ht1_ratio,tau1);
    fs2(i,1) = disc_spring_force(x(i,3) - x(i,5),ht2_ratio,tau2);
end

figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(t,fs1,'LineWidth',2,'Color','Blue')
hold on
plot(t,fs2,'LineWidth',2,'Color','Red')
title("Force vs time","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$F(t)$',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(t,fs1 - fs2,'LineWidth',2,'Color','Black')
title("Difference in Force vs time","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$F(t)$',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

%% Questions - 

% 1) If multiple solutions for x1 given x_stack (refer static case for
% notation), then how do we determine the actual dynamics 
% 2) Snap-through events and other non-linear behaviour (depends on ICs,
% etc)
