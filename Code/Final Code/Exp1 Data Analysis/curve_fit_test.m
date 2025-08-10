%% Force-deflection model of the 2-spring stack based on curve-fitting

%% Clear and Close all windows
clear;
close all;

%% Read Data

data_forward = readtable('static_forward.xlsx','Sheet','F-s Measurement');
data_reverse = readtable('static_reverse.xlsx','Sheet','F-s Measurement');

x_f = data_forward.Laser_mm_;
F_f = data_forward.LoadCell1kN;

x_r = data_reverse.Laser_mm_;
F_r = data_reverse.LoadCell1kN;

x_f = x_f + 3.25;
x_r = x_r + 3.25;

%% Plot Experimental Data

figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(x_f, F_f, 'LineWidth', 1, 'Color', 'Cyan')
hold on
plot(x_r, F_r, 'LineWidth', 1, 'Color', '#FFA500')
hold on
xlim([0,2.5])
title( "Static Force-Deflection curve","Interpreter","latex","FontSize",15)
xlabel('$x_{st}$ [mm]',"Interpreter","latex","FontSize",15)
ylabel('Force',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

Colors = '#013220';

%% Key Points

x_key = [0, 0.241, 0.5, 1, 1.5, 2, 2.5];
F_key = [0, 41.8, 80.5, 111, 116, 123, 198];

x_key_2 = [0, 0.241, 0.5, 1, 1.5, 1.75, 2, 2.27, 2.5];
F_key_2 = [0, 37, 76, 106, 110, 112, 118, 148, 192];

%% Interpolate using Cubic splines

x_interp = 0:0.01:2.5;
F_interp = spline(x_key,F_key,x_interp);
F_interp_2 = spline(x_key_2, F_key_2, x_interp);

%% Plot force-deflection curve based on curve-fit

scatter(x_key, F_key, 50,'Red','filled');
hold on
plot(x_interp,F_interp,'LineWidth',2,'Color','Red');
hold on
scatter(x_key_2, F_key_2, 50,'Blue','filled');
hold on
plot(x_interp,F_interp_2,'LineWidth',2,'Color','Blue');
hold on

%% Test function

x_trial = 1.8; % Trial value of deflection in mm
F_trial = disc_spring_force_exp(x_trial); % Force exerted by the spring corresponding to the trial deflection value
scatter(x_trial, F_trial, 50,'black','filled');

