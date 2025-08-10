%% Code to simulate the dynamics of the system iteratively for different frequencies

% For conditions in experiment 1

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
f_excit_array = 4:2:26; % Array of base excitation frequency values (in Hz)

dx_arr = zeros(length(tspan),4); % Placeholder array for the state evolution

MT_array = zeros(1,length(f_excit_array)); % Place holder array for motion transmissibility values
Amp_array = zeros(1,length(f_excit_array)); % Placeholder array of Amplitude (pk-to-pk) values 


%% Settings - 

% VISCOUS DAMPING
% 1) m = 11.2 kg, c = 0.05 Ns/mm, h1/t = h2/t = 1.41
% 2) m = 11.2 kg, c = 0.1 Ns/mm, h1/t = h2/t = 1.41
% 3) m = 11.2 kg, c = 0.075 Ns/mm, h1/t = h2/t = 1.41
% 4) m = 11.2 kg, c = 0.5 Ns/mm, h1/t = h2/t = 1.41
% 5) m = 11.2 kg, c = 0.2 Ns/mm, h1/t = h2/t = 1.41
% 6) m = 10.9 kg, c = 0.05 Ns/mm, h1/t = h2/t = 1.41
% 7) m = 10.9 kg, c = 0.1 Ns/mm, h1/t = h2/t = 1.41
% 8) m = 10.9 kg, c = 0.2 Ns/mm, h1/t = h2/t = 1.41
% 9) m = 11.1 kg, c = 0.05 Ns/mm, h1/t = h2/t = 1.41
% 10) m = 11.1 kg, c = 0.1 Ns/mm, h1/t = h2/t = 1.41
% 11) m = 11.1 kg, c = 0.01 Ns/mm, h1/t = h2/t = 1.41
% 12) m = 11.2 kg, c = 0.01 Ns/mm, h1/t = h2/t = 1.41
% 13) m = 11.2 kg, c = 0.5 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05
% 14) m = 11.2 kg, c = 0.4 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05
% 15) m = 11.2 kg, c = 0.3 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05
% 16) m = 11.2 kg, c = 0.2 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05
% 17) m = 11.2 kg, c = 0.5 Ns/mm, h1/t = h2/t = 1.41, scale = 1.03
% 18) m = 11.2 kg, c = 0.5 Ns/mm, Experimental curve-fit model
% 19) m = 11.2 kg, c = 0.2 Ns/mm, Experimental curve-fit model
% 20) m = 11.2 kg, c = 0.35 Ns/mm, Experimental curve-fit model
% 21) m = 11.2 kg, c = 0.5 Ns/mm, h1/t = h2/t = 1.41, scale = 1.04
% 22) m = 11.2 kg, c = 0.1 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05
% 23) m = 11.2 kg, c = 0.35 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05
% 24) m = 11.2 kg, c = 0.45 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05

% MIXED DAMPING
% 1) m = 11.2 kg, c = 0.4 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.1
% 2) m = 11.2 kg, c = 0.4 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.5
% 3) m = 11.2 kg, c = 0.05 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.05, Fn = mg
% 4) m = 11.2 kg, c = 0.2 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.05, Fn = mg
% 5) m = 11.2 kg, c = 0.2 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.03, Fn = mg
% 6) m = 11.2 kg, c = 0.1 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.03, Fn = mg

% COULOMB FRICTION
% 1) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.1, Fn = 1
% 2) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.5, Fn = 1
% 3) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.6, Fn = 1
% 4) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.7, Fn = 1
% 5) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 1.0, Fn = 1
% 6) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.1, Fn = mg
% 7) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.05, Fn = mg
% 8) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.02, Fn = mg
% 9) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.07, Fn = mg
% 10) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.04, Fn = mg
% 11) m = 11.2 kg, c = 0.0001 Ns/mm, h1/t = h2/t = 1.41, scale = 1.05, mu = 0.03, Fn = mg

%% Simulate the dynamics of the system for different base excitation freqeuncies

% f_no is the frequency sweep index (1 corresponds to 4 Hz, 2 corresponds to 6 Hz)
for i=1:12

    f_excit = f_excit_array(i); % Base excitation frequency in Hz
 
    % Solve for static equilibrium
    x_eq = static_eqb(ht1_ratio, tau1, ht2_ratio, tau2, m1*g);
    
    % Set initial state
    x_base_init = 0; % Initial displacement of the base (mm)
    U_init = 0;
    x0 = [x_eq(1);0;x_base_init;U_init]; % Initial state vector (Static equilibrium state)
    
    % Solve the ODEs to obtain the displacements and velocities numerically 
    options=odeset('abstol',1e-9,'reltol',1e-9); % Tolerances
    [t,x] = ode89(@(t,x)dyn_dspring_stack(t,x,m1,g,ht1_ratio,tau1,ht2_ratio,tau2,A,f_excit),tspan,x0,options);


    % Obtain the state evolution vector for each time instant
    for j=1:length(t)
        dx_arr(j,:) = dyn_dspring_stack(t(j),x(j,:).',m1,g,ht1_ratio,tau1,ht2_ratio,tau2,A,f_excit).';
    end

    acc = dx_arr(:,2); % Acceleration of the mass
    acc_base = -(2 * pi * f_excit)^2 * x(:,3); % Acceleration of the base

    acc_ss = acc(floor(end/2):end,1); % Mass acceleration values in steady-state

    % Caclulate and store Motion Transmissibility Values
    acc_base_ptp = max(acc_base) - min(acc_base);
    acc_ss_ptp = max(acc_ss) - min(acc_ss);

    MT = acc_ss_ptp/acc_base_ptp;
    dB = 20*log10(MT);
    MT_array(i) = dB;

    Amp_array(i) = acc_ss_ptp;

    % Export state and acceleration data as csv file
    export_array = [t, x, acc, acc_base];
    writematrix(export_array,'/Users/jacobsony/Desktop/try3/a_sim_'+ string(f_excit)+'.csv')

    % Plot acceleration profiles and save them 
    num_cyc = 6;
    start_t = max(t) - (num_cyc/f_excit);
    end_t = max(t);
    min_y = min(min(acc),min(acc_base));
    max_y = max(max(acc),max(acc_base));

    T1 = tiledlayout(1,1);
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

    plotpath = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Dynamic/Sim6', 'acc_' + string(f_excit) + '_sim.pdf');
    %exportgraphics(T1,plotpath,'Resolution',300)

end



%% Plot Motion Transmissibility Curve

% Makima Interpolation of the Motion Transmissibility Curves
f_excit_array_interp = 4:0.1:26;
MT_array_interp = makima(f_excit_array, MT_array, f_excit_array_interp);

T2 = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_excit_array_interp, MT_array_interp, 'LineWidth', 2, 'Color', 'Blue')
hold on
scatter(f_excit_array, MT_array, 50,'Blue','filled')
title( "Motion Transmissibility","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Motion Transmissibility [dB]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

plotpath2 = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Dynamic/Sim6/MT_curve.pdf');
exportgraphics(T2,plotpath2,'Resolution',300)

%% Plot Peak-to-Peak Amplitude Curve

% Makima Interpolation of the Peak-to-Peak Amplitude Curve
Amp_array_interp = makima(f_excit_array, Amp_array, f_excit_array_interp);

T3 = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_excit_array_interp, Amp_array_interp, 'LineWidth', 2, 'Color', 'Blue') 
hold on
scatter(f_excit_array, Amp_array, 50,'Blue','filled')
title( "Peak-to-Peak Amplitude","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Peak-to-Peak Amplitude [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

plotpath3 = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Dynamic/Sim6/Amp_curve.pdf');
exportgraphics(T3,plotpath3,'Resolution',300)

%% Close All windows

close all;

%% Observations 

% Increase c, Resonance frequency shifts to the left
% Increase c, Steady state amplitude increases (woah !) (MT at high freq
% low but MT at low freq is higher !)
% Increase effective k, Resonance frequency shifts to the right
% Intuition, need to increase c and increase effective k to get closer match with experimental results
% As of now, best match when m = 11.2 kg, c = 0.35 - 0.5 Ns/mm, h1/t = h2/t = 1.41 and scale = 1.05
% Experimental results - Higher frequency components at higher base excitation frequencies
% Effect of Coulomb friction relative to viscous damping on the acceleration decreases at higher base excitation frequencies