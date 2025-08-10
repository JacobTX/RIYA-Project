%% Motion Transmissibility Analysis in Simulation

%% Simulation Parameters
% m = 10.7 kg (for Non-QZS), 11.2 kg (for QZS)
% h1/tau = h2/tau = 1.41, E = 200 GPa
% c = 0.05 Ns/mm, mu = 0
% Peak-to-Peak amplitude = 0.1 mm

% For comparison of QZS vs Non-QZS motion transmissibility

% Non-QZS Regime
f_NQZS = [5,7,8,9,11,13];
MT_NQZS = [4.91, 16.65, 11.37, 6.23, -1.05, -5.72];

% QZS Regime
f_QZS = [1,2,3,5,7];
MT_QZS = [5.42, 17.48, -11.06, -21.49, -27.41];

% Interpolation (Modified Akima Interpolation)

f_NQZS_interp = 5:0.1:13;
MT_NQZS_interp = makima(f_NQZS, MT_NQZS, f_NQZS_interp);

f_QZS_interp = 1:0.1:7;
MT_QZS_interp = makima(f_QZS, MT_QZS, f_QZS_interp);

% Plot
figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(f_NQZS_interp, MT_NQZS_interp, 'LineWidth', 2, 'Color', 'Blue') % Non-QZS Transmissibility Curve
hold on
plot(f_QZS_interp, MT_QZS_interp, 'LineWidth', 2, 'Color', 'Red') % QZS Transmissibility Curve
hold on
scatter(f_QZS, MT_QZS, 50,'Red','filled')
hold on
scatter(f_NQZS, MT_NQZS, 50,'Blue','filled')

legend("Non-QZS Regime","QZS Regime", "Interpreter","latex","FontSize",15)
title( "Motion Transmissibility Curves","Interpreter","latex","FontSize",15)
xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
ylabel('Motion Transmissibility [dB]',"Interpreter","latex","FontSize",15)
set(gca, 'FontSize',15)
