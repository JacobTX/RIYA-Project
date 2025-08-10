%% Dynamic Analysis of Experiment 2 Data

%% Clear and Close All Windows

clear;
close all;

%% Iteratively read all the files, plot the accelerometer data and aclculate Motion Transmissibility

f_arr = 4:2:20; % Array of frequency values
MT_arr = zeros(1,length(f_arr)); % Placeholder array of motion transmissibility values
Amp_arr = zeros(1,length(f_arr)); % Placeholder array of Amplitude (pk-to-pk) values 

for i = 4:4
    % Extract the Data
    filename = 'set' + string(i) + '.xlsx'; % File Name corresponding to the File Number (Use appropriate file path)
    data = readtable(filename,'Sheet','Vibration');
    f = f_arr(i); % Frequency of base excitation corresponding to the file number

    a_base = data.pcb4; % Accelerometer data corresponding to the base (in g)
    a_mass = data.pcb2; % Accelerometer data corresponding to the mass (in g)

    len = min(length(a_base),length(a_mass));

    a_base = 9.81 * 1000 * a_base(1:len); % in mm/s^2
    a_mass = 9.81 * 1000 * a_mass(1:len); % in mm/s^2
    t = 0:(1/2000):((len-1)/2000);

    num_cyc = 6; % Number of cycles to be displayed in the plot
    start_t = max(t) - num_cyc * (1/f);
    end_t = max(t);

    % Calculate Motion Transmissibility using acceleration data in Time domain
    a_base_ptp = max(a_base) - min(a_base);
    a_mass_ptp = max(a_mass) - min(a_mass);

    MT_t = 20*log10(a_mass_ptp/a_base_ptp);

    MT_arr(i) = MT_t;
    Amp_arr(i) = a_mass_ptp;
    
    % Plot accelerometer data
    F1 = tiledlayout(1,1);
    set(gcf,'Position',[0,0,800,500],'defaultAxesTickLabelInterpreter','latex'); 
    
    nexttile
    set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
    
    plot(t - start_t, a_mass, 'LineWidth', 2, 'Color', 'Blue')
    hold on
    plot(t - start_t, a_base, 'LineWidth', 2, 'Color', 'Red')
    hold on
    legend("$\ddot{x}_{st}$","$\ddot{x}_{base}$","Interpreter","latex","FontSize",15)
    title( "Accelerometer Data vs Time","Interpreter","latex","FontSize",15)
    xlim([0, end_t - start_t])
    ylim([-1000,800])
    xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
    ylabel('Acceleration [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
    set(gca, 'FontSize',15)
    
    plotpath = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Data Analysis/Exp', 'acc_' + string(f) + '_exp.pdf');
    %exportgraphics(F1,plotpath','Resolution',300)
end

%% Plot Motion Transmissibility curve

% Makima Interpolation of the Motion Transmissibility Curves
f_arr_interp = 4:0.1:20;
MT_arr(1:2) = [0 0.5]; % Correction since data is not pre-processed
MT_arr_interp = makima(f_arr, MT_arr, f_arr_interp);

% Based on Simulation Data


% For Experimental Data 2
f_arr_2 = 4:2:26;
f_arr_2_interp = 4:0.1:26;

MT_arr_2 = [0.2429    0.5523    0.9968    1.5901    2.3709    3.3177  4.4612    5.8065    7.2423    8.2602    8.1276    6.8076];
MT_arr_3 = [0.2658    0.6045    1.0911    1.7445    2.5986    3.6282    4.8668    6.2655    7.5713    8.0510    7.2203    5.5559];
MT_arr_4 = [0.2052    0.3878    0.5973    1.0089    2.5543    4.6137    5.9980    7.8093   11.0573   18.1993   14.8238   11.1735];
MT_arr_5 = [0.2159    0.4270    0.6782    1.1479    2.4665    3.9255    5.1184    6.6031    8.7362   11.5762   11.7702    9.2366];
%MT_arr_6 = [0.1387    0.2028    0.2520    0.3154    0.4654    1.0386    4.4314    7.2115   10.2921   18.2813   14.8878   11.2451];


MT_arr_2_interp = makima(f_arr_2,MT_arr_2,f_arr_2_interp);
MT_arr_3_interp = makima(f_arr_2,MT_arr_3,f_arr_2_interp);
MT_arr_4_interp = makima(f_arr_2,MT_arr_4,f_arr_2_interp);
MT_arr_5_interp = makima(f_arr_2,MT_arr_5,f_arr_2_interp);

F2 = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_arr_interp, MT_arr_interp, 'LineWidth', 2, 'Color', 'Blue')  
hold on
plot(f_arr_2_interp, MT_arr_2_interp, 'LineWidth', 2, 'Color', '#A020F0') 
hold on
plot(f_arr_2_interp, MT_arr_4_interp, 'LineWidth', 2, 'Color', 'Red') 
hold on
plot(f_arr_2_interp, MT_arr_5_interp, 'LineWidth', 2, 'Color', '#037D50') 
hold on
plot(f_arr_interp, MT_arr_interp + 3, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
hold on
plot(f_arr_interp, MT_arr_interp - 3, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
hold on
scatter(f_arr, MT_arr, 50,'Blue','filled')
hold on
scatter(f_arr_2, MT_arr_4, 50,'Red','filled')
title( "Motion Transmissibility","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Motion Transmissibility [dB]',"Interpreter","latex","FontSize",15)
legend("Experimental","$c = 0.5, \mu = 0$","$c = 10^{-4}, \mu = 0.04$","$c = 0.2, \mu = 0.03$","Interpreter","latex","FontSize",15,'Location','Best')
set(gca, 'FontSize',15)

plotpath2 = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Data Analysis/Exp/MT_curve.pdf');
exportgraphics(F2,plotpath2,'Resolution',300)

%% Plot Peak-to-Peak Amplitude Curve

% Makima Interpolation of the Peak-to-Peak Amplitude Curve
Amp_arr_interp = makima(f_arr, Amp_arr, f_arr_interp);

F3 = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_arr_interp, Amp_arr_interp, 'LineWidth', 2, 'Color', 'Blue') 
hold on
scatter(f_arr, Amp_arr, 50,'Blue','filled')
title( "Peak-to-Peak Amplitude","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Peak-to-Peak Amplitude [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

plotpath3 = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Data Analysis/Exp/Amp_curve.pdf');
%exportgraphics(F3,plotpath3,'Resolution',300)


%% Close Windows

close all;

%% FFT Analysis

% FFT for input and output as well in experimental data
% MT from FFT based on dominant peak of output vs input

% Remove DC component from accelerometer data
%a_b = a_base - mean(a_base);
%a_m = a_mass - mean(a_mass);

% FFT Spectrum of Accelerometer data corresponding to the mass
%fft_analysis(t,a_b.');

% FFT Spectrum of Accelerometer data corresponding to the mass
%fft_analysis(t,a_m.');

% Calculation using acceleration data in Frequency domain using the
% dominant peak

%% Sources of error

% Harmonic distortion due to loading effect
% The input displacement is not sinousoidal, a more accurate model would be
% to consider the dynamics of the base coupled with that of the mass

% Figure out Damping model that works !
% Non-linear damping = k based on simulation / experimental data

% Use experimental MT to find c 
% Coulomb friction - read paper, incorporate (increase coulomb, keep small
% amount of viscous)
% Simulation for new experimental data

% Beyond 20 Hz, system behaves more like 2 DoF system so theory breaks down
% possibly !

% computational issues, integration method, equations of motion,  
