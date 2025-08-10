%% Simulation of Simple Linear Spring - Mass system in MATLAB

% System Parameters
m1 = 10.7;
m2 = 1e-3;
b1 = 10;
b2 = 10;
g = 9.81;
l0 = 3.5;
x_base_init = 0;

% Temporal parameters
dt = 0.001; % Time step
Tmax = 100; % Maximum time
tspan = 0:dt:Tmax; % Duration of the simulation

% Initial state vector
x0 = [0;0;0;0;x_base_init];

%% Solve the ODEs to obtain the states numerically 

options=odeset('abstol',1e-9,'reltol',1e-9);
[t,x] = ode45(@(t,x)dyn_spring_stack(t,x,m1,m2,g,b1,b2),tspan,x0,options);


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
for i=1:100:length(tspan)
    subplot(2,2,[2 4])
    plot([0,0],[x_base_init - x(i,5),l0 - x(i,1)],'LineWidth',3,'Color','Blue')
    hold on
    rectangle('Position',[-0.3 l0-x(i,1) 0.6 0.6],'LineWidth',3,'FaceColor','#853A00')
    hold on
    rectangle('Position',[-1 x_base_init-x(i,5)-0.3 2 0.3],'LineWidth',3,'FaceColor','#808090')
    hold off
    axis([-2,2,-2,10])
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