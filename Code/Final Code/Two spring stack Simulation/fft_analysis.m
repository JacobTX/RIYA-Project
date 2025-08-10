%% Function for FFT analysis of a signal with calibration

% References - https://www.mathworks.com/help/matlab/ref/fft.html

function [] = fft_analysis(t,a)

    dt = t(2) - t(1); % Time step (s)
    Fs = 1/dt; % Sampling frequency (Hz)
    LT = length(t); % Number of elements in the time vector

    % Set the number of FFT points (N)
    N = 2^(nextpow2(LT)+1);
    %N = LT;

    df = Fs/N; % Frequency resolution (Hz)

    disp("DSP Parameters")
    disp(["Time resolution", dt])
    disp(["Time record", max(t)])
    disp(["Sampling frequency", Fs])
    disp(["Frequency resolution", df])
    disp(["Number of points", N])

    a = a.*hanning(LT)'; % Windowing to prevent spectral leakage
    a = [a zeros(1, N - LT)]; % Zero Padding to increase frequency resolution
    S = fft(a);
    
    % Calibration
    % Scaling (works for perfect signals with integer number of periods)
    S1 = S(1:floor(N/2));
    f = df*(0:(N/2)-1); % Frequency vector
    S2 = abs(S1)/(N/2);
    k = 6.6; % For windowed signals ??
    S2 = k * S2;
    
    figure;
    set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
    plot(f,S2,'LineWidth',3,'Color','Blue');
    title("FFT Spectrum","Interpreter","latex","FontSize",15)
    xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
    ylabel('$a$ [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
    %ylabel('$a$ [$g$]',"Interpreter","latex","FontSize",15)
    %xlim([0,max(f)])
    xlim([0,100])
    set(gca, 'FontSize',15)

    figure;
    set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
    plot(f,20*log10(S2),'LineWidth',3,'Color','Blue');
    title("FFT Spectrum","Interpreter","latex","FontSize",15)
    xlabel('$f$ [Hz]',"Interpreter","latex","FontSize",15)
    ylabel('dB',"Interpreter","latex","FontSize",15)
    %xlim([0,max(f)])
    xlim([0,100])
    set(gca, 'FontSize',15)

end