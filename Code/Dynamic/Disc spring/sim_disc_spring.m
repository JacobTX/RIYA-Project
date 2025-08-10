%% Simulation of Single Disc Spring - Mass system in MATLAB

% System Parameters
m = 10.7; % Mass
b = 10; % Damping coefficient
g = 9.81; % Acceleration due to gravity

ht_ratio = 1.4; % ratio of height to thickness
tau = 0.5; % thickness (mm)
h0 = ht_ratio * tau; % height (mm)
l0 = 5*h0; % Natural length of the spring + spacers
x_base_init = 0; % Initial displacement of the base

% Temporal parameters
dt = 0.001; % Time step
Tmax = 20; % Maximum time
tspan = 0:dt:Tmax; % Duration of the simulation

% Initial state vector
x0 = [0;0;x_base_init];

%% Solve the ODEs to obtain the states numerically 

options=odeset('abstol',1e-9,'reltol',1e-9);
[t,x] = ode45(@(t,x)dyn_disc_spring(t,x,m,g,b,ht_ratio,tau),tspan,x0,options);


%% Interface for the plots and animation

set(gcf,'Position',[0,0,1500,1000],'defaultAxesTickLabelInterpreter','latex'); 

% Displacement plot
subplot(2,2,1)
plot(tspan,x(:,1),'LineWidth',2,'Color','Blue')
title("Displacement profile","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$x$ [mm]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

% Velocity plot
subplot(2,2,3)
plot(tspan,x(:,2),'LineWidth',2,'Color','Red')
title("Velocity profile","Interpreter","latex","FontSize",15)
xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
ylabel('$\dot{x}$ [mm/s]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

% Animation
for i=1:10:length(tspan)
    subplot(2,2,[2 4])
    plot([-17.25,-17.25],[x_base_init - x(i,3), 2*h0 - x(i,3)],'LineWidth',10,'Color','	#5A5A5A')
    hold on
    plot([17.25,17.25],[x_base_init - x(i,3), 2*h0 - x(i,3)],'LineWidth',10,'Color','	#5A5A5A')
    hold on
    plot([-11.2,-11.2],[l0 - x(i,1) - 2*h0, l0 - x(i,1)],'LineWidth',10,'Color','	#5A5A5A')
    hold on
    plot([11.2,11.2],[l0 - x(i,1) - 2*h0, l0 - x(i,1)],'LineWidth',10,'Color','	#5A5A5A')
    hold on
    plot([-17.25,-11.2],[2*h0 - x(i,3), l0 - x(i,1) - 2*h0],'LineWidth',10,'Color','Blue')
    hold on
    plot([17.25,11.2],[2*h0 - x(i,3), l0 - x(i,1) - 2*h0],'LineWidth',10,'Color','Blue')
    hold on
    rectangle('Position',[-12.5 l0-x(i,1) 25 1],'LineWidth',3,'FaceColor','#853A00')
    hold on
    rectangle('Position',[-20 x_base_init-x(i,3)-0.3 40 0.3],'LineWidth',3,'FaceColor','#808090')
    hold off
    axis([-35,35,-2,6])
    title("Spring-Mass Animation","Interpreter","latex","FontSize",15)
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