%% Dynamic Analysis of Experiment 1 Data


%% Clear and Close All Windows

clear;
close all;

%% Iteratively read all the files, plot the accelerometer data and calculate Motion Transmissibility

f_arr = 4:2:26; % Array of frequency values
MT_arr = zeros(1,length(f_arr)); % Placeholder array of motion transmissibility values
Amp_arr = zeros(1,length(f_arr)); % Placeholder array of Amplitude (pk-to-pk) values 

for i = 1:12
    % Extract the Data
    filename = 'set' + string(i) + '.xlsx'; % File Name corresponding to the File Number (Use appropriate file path)
    data = readtable(filename,'Sheet','Vibration');
    f = f_arr(i); % Frequency of base excitation corresponding to the file number

    a_base = data.pcb4; % Accelerometer data corresponding to the base (in g)
    a_mass = data.pcb5; % Accelerometer data corresponding to the mass (in g)

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
    xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
    ylabel('Acceleration [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
    set(gca, 'FontSize',15)
    
    plotpath = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Data Analysis/Exp', 'acc_' + string(f) + '_exp.pdf');
    %exportgraphics(F1,plotpath','Resolution',300)
end

%% Interpolation of Experimental Data

% Makima Interpolation of the Motion Transmissibility Curves
f_arr_interp = 4:0.1:26;
MT_arr_interp = makima(f_arr, MT_arr, f_arr_interp);

%% Interpolation of Simulation Data

% Simulation data is obtained after first running the simulation and
% copying the MT array values to the file

% Viscous Damping

% c = 0.35 Ns/mm
MT_arr_v1 = [2.5728    5.4971    3.5577   -0.0768   -3.0435   -5.4081 -7.3216   -8.9034  -10.2224  -11.3866  -12.3893  -13.2785]; 
MT_arr_interp_v1 = makima(f_arr,MT_arr_v1,f_arr_interp);

% c = 0.4 Ns/mm
MT_arr_v2 = [2.4142    4.7777    3.2931    0.0625   -2.6742   -4.8793 -6.6747   -8.1628   -9.4143  -10.5142  -11.4707  -12.3224]; 
MT_arr_interp_v2 = makima(f_arr,MT_arr_v2,f_arr_interp);

% c = 0.45 Ns/mm
MT_arr_v3 = [2.2747    4.1999    3.0056    0.1596   -2.3464   -4.3957   -6.0761   -7.4813   -8.6617   -9.7188  -10.6371  -11.4591];
MT_arr_interp_v3 = makima(f_arr,MT_arr_v3,f_arr_interp);

% c = 0.5 Ns/mm
MT_arr_v4 = [2.1478    3.7218    2.7133    0.2252   -2.0569   -3.9565   -5.5311   -6.8584   -7.9907   -8.9962   -9.8815  -10.6791];
MT_arr_interp_v4 = makima(f_arr,MT_arr_v4,f_arr_interp);

% Coulomb Friction (Friction force magnitude = mu * m * g)

% c !=0 for some numerical stability

% c = 0.0001 Ns/mm, mu = 0.02 
MT_arr_c1 = [0.7740   15.5273   17.7507   19.4434   26.3035   -0.9123   -2.6362   -9.6428  -14.3745  -16.1213  -17.7086  -19.1547];
MT_arr_interp_c1 = makima(f_arr,MT_arr_c1,f_arr_interp);

% c = 0.0001 Ns/mm, mu = 0.03
MT_arr_c2 = [0.2727    3.6364    4.9633    2.2697    4.8226    5.0399   -4.7386   -9.0554  -13.3085  -15.0610  -16.6368  -18.0739];
MT_arr_interp_c2 = makima(f_arr,MT_arr_c2,f_arr_interp);

% c = 0.0001 Ns/mm, mu = 0.04 
MT_arr_c3 = [0.1385    0.3178    4.8717    1.9242    2.1246   -0.2822   -4.3041   -8.6037  -12.4413  -14.1179  -15.6923  -17.1141];
MT_arr_interp_c3 = makima(f_arr,MT_arr_c3,f_arr_interp);

% c = 0.0001 Ns/mm, mu = 0.05 
MT_arr_c4 = [0.0843    0.1190    2.8102    1.4215    2.0661   -0.5236   -3.9823   -8.3238  -11.5642  -13.2708  -14.8258  -16.2513];
MT_arr_interp_c4 = makima(f_arr,MT_arr_c4,f_arr_interp);

% c = 0.0001 Ns/mm, mu = 0.07 
MT_arr_c5 = [0.0416    0.0428   -0.0079    0.7822   -0.3722   -1.0583   -3.5662   -8.2420  -10.1498  -11.7974  -13.3360  -14.7507];
MT_arr_interp_c5 = makima(f_arr,MT_arr_c5,f_arr_interp);

% c = 0.0001 Ns/mm, mu = 0.1 
MT_arr_c6 = [0.0196    0.0177    0.0114   -0.2218   -0.9104   -2.4473   -4.5340   -6.5219   -8.3673   -9.9882  -11.5082  -12.9038];
MT_arr_interp_c6 = makima(f_arr,MT_arr_c6,f_arr_interp);


% Mixed Friction

% c = 0.05 Ns/mm, mu = 0.05
MT_arr_m1 = [0.0825    0.1141    1.3460    1.0191    0.6682   -1.1215  -4.0599   -9.3067  -11.1550  -12.8304  -14.3400  -15.7188];
MT_arr_interp_m1 = makima(f_arr,MT_arr_m1,f_arr_interp);

% c = 0.2 Ns/mm, mu = 0.05
MT_arr_m2 = [0.0775    0.1014    0.3041    0.5900   -1.5508   -3.8484   -5.9173   -7.7450   -9.3465  -10.7974  -12.0844  -13.2467];
MT_arr_interp_m2 = makima(f_arr,MT_arr_m2,f_arr_interp);

% c = 0.2 Ns/mm, mu = 0.03
MT_arr_m3 = [0.2317    0.9744    3.4067    0.6077   -1.1367   -2.3326 -7.1396   -9.0102  -10.6233  -12.0670  -13.3371  -14.4757];
MT_arr_interp_m3 = makima(f_arr,MT_arr_m3,f_arr_interp);

% c = 0.1 Ns/mm, mu = 0.03
MT_arr_m4 = [0.2508    1.5839    4.5055    1.5946    1.5096   -1.2166   -4.5525  -10.3273  -12.1204  -13.7411  -15.1791  -16.4804];
MT_arr_interp_m4 = makima(f_arr,MT_arr_m4,f_arr_interp);

% c = 0.4 Ns/mm, mu = 0.001
MT_arr_m5 = [2.4342    4.6455    3.2788    0.0916   -2.6265   -4.8256 -6.6193   -8.1081   -9.3610  -10.4638  -11.4228  -12.2767];
MT_arr_interp_m5 = makima(f_arr,MT_arr_m5,f_arr_interp);

%% Plot Motion Transmissibility Curve for different Damping coefficients

F2 = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_arr_interp, MT_arr_interp, 'LineWidth', 2, 'Color', 'Blue')  
hold on
plot(f_arr_interp, MT_arr_interp_v1, 'LineWidth', 2, 'Color', '#FFA500') 
hold on
plot(f_arr_interp, MT_arr_interp_v2, 'LineWidth', 2, 'Color', '#037D50') 
hold on
plot(f_arr_interp, MT_arr_interp_v3, 'LineWidth', 2, 'Color', '#A020F0') 
hold on
plot(f_arr_interp, MT_arr_interp_v4, 'LineWidth', 2, 'Color', 'Red') 
hold on
plot(f_arr_interp, MT_arr_interp + 5, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
hold on
plot(f_arr_interp, MT_arr_interp - 5, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
hold on
scatter(f_arr, MT_arr, 50,'Blue','filled')
hold on
scatter(f_arr, MT_arr_v4, 50,'Red','filled')
title( "Motion Transmissibility","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Motion Transmissibility [dB]',"Interpreter","latex","FontSize",15)
legend("Experimental","$c = 0.35$","$c = 0.40$","$c = 0.45$","$c = 0.50$","Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

plotpath2 = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Dynamic/MT_curve_visc.pdf');
exportgraphics(F2,plotpath2,'Resolution',300)

%% Plot Motion Transmissibility Curve for mu values of Coulomb friction

F3 = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_arr_interp, MT_arr_interp, 'LineWidth', 2, 'Color', 'Blue')  
hold on
plot(f_arr_interp, MT_arr_interp_c2, 'LineWidth', 2, 'Color', '#037D50') 
hold on
plot(f_arr_interp, MT_arr_interp_c3, 'LineWidth', 2, 'Color', 'Red') 
hold on
plot(f_arr_interp, MT_arr_interp_c4, 'LineWidth', 2, 'Color', '#A020F0') 
hold on
plot(f_arr_interp, MT_arr_interp_c5, 'LineWidth', 2, 'Color', '#FFA500') 
hold on
plot(f_arr_interp, MT_arr_interp_c6, 'LineWidth', 2, 'Color', 'Magenta') 
hold on
plot(f_arr_interp, MT_arr_interp + 5, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
hold on
plot(f_arr_interp, MT_arr_interp - 5, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
hold on
scatter(f_arr, MT_arr, 50,'Blue','filled')
hold on
scatter(f_arr, MT_arr_c3, 50,'Red','filled')
title( "Motion Transmissibility","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Motion Transmissibility [dB]',"Interpreter","latex","FontSize",15)
legend("Experimental","$\mu = 0.03$","$\mu = 0.04$","$\mu = 0.05$","$\mu = 0.07$","$\mu = 0.1$","Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

plotpath3 = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Dynamic/MT_curve_coul.pdf');
exportgraphics(F3,plotpath3,'Resolution',300)


%% Plot Motion Transmissibility Curve for Mixed conditions

F4 = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_arr_interp, MT_arr_interp, 'LineWidth', 2, 'Color', 'Blue')  
hold on
plot(f_arr_interp, MT_arr_interp_m1, 'LineWidth', 2, 'Color', '#037D50') 
hold on
plot(f_arr_interp, MT_arr_interp_m2, 'LineWidth', 2, 'Color', '#A020F0') 
hold on
plot(f_arr_interp, MT_arr_interp_m3, 'LineWidth', 2, 'Color', 'Red') 
hold on
plot(f_arr_interp, MT_arr_interp_m4, 'LineWidth', 2, 'Color', '#FFA500') 
hold on
plot(f_arr_interp, MT_arr_interp + 5, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
hold on
plot(f_arr_interp, MT_arr_interp - 5, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
hold on
scatter(f_arr, MT_arr, 50,'Blue','filled')
hold on
scatter(f_arr, MT_arr_m3, 50,'Red','filled')
title( "Motion Transmissibility","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Motion Transmissibility [dB]',"Interpreter","latex","FontSize",15)
legend("Experimental","$c = 0.05, \mu = 0.05$","$c = 0.2, \mu = 0.05$","$c = 0.2, \mu = 0.03$","$c = 0.1, \mu = 0.03$","Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

plotpath4 = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Dynamic/MT_curve_mix.pdf');
exportgraphics(F4,plotpath4,'Resolution',300)

%% Plot Peak-to-Peak Amplitude Curve

% Makima Interpolation of the Peak-to-Peak Amplitude Curve
Amp_arr_interp = makima(f_arr, Amp_arr, f_arr_interp);

F5 = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_arr_interp, Amp_arr_interp, 'LineWidth', 2, 'Color', 'Blue') 
hold on
scatter(f_arr, Amp_arr, 50,'Blue','filled')
title( "Peak-to-Peak Amplitude","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Peak-to-Peak Amplitude [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

plotpath5 = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Data Analysis/Exp/Amp_curve.pdf');
%exportgraphics(F5,plotpath5,'Resolution',300)

%% Compare different damping models

F6 = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_arr_interp, MT_arr_interp, 'LineWidth', 2, 'Color', 'Blue')
hold on
plot(f_arr_interp, MT_arr_interp_v4, 'LineWidth', 2, 'Color', 'Red') 
hold on
plot(f_arr_interp, MT_arr_interp_c3, 'LineWidth', 2, 'Color', '#037D50') 
hold on
plot(f_arr_interp, MT_arr_interp_m3, 'LineWidth', 2, 'Color', '#A020F0') 
hold on
plot(f_arr_interp, MT_arr_interp + 5, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
hold on
plot(f_arr_interp, MT_arr_interp - 5, 'LineWidth', 2, 'Color', 'Blue', 'LineStyle',':')
scatter(f_arr, MT_arr, 50,'Blue','filled')
title( "Motion Transmissibility","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Motion Transmissibility [dB]',"Interpreter","latex","FontSize",15)
legend("Experimental","$c = 0.5, \mu = 0$","$c = 10^{-4}, \mu = 0.04$","$c = 0.2, \mu = 0.03$","Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)

plotpath6 = fullfile('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Dynamic/MT_curve_comparison.pdf');
exportgraphics(F6,plotpath6,'Resolution',300)

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