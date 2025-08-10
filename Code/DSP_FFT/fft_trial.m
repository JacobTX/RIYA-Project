%% FFT Trial

%% Case 1

fs = 500;
duration = 3;
N = fs * duration;
t = 0:1/fs:duration-1/fs;

a1 = 3; f1 = 30; phi1 = 0.6;
a2 = 2; f2 = 45; phi2 = -0.8;
a3 = 1; f3 = 70; phi3 = 2;

s1 = a1*cos(2*pi*f1*t + phi1);
s2 = a2*cos(2*pi*f2*t + phi2);
s3 = a3*cos(2*pi*f3*t + phi3);

s = s1 + s2 + s3 + 0.5*(randn(size(t)));
%s = s1 + s2 + s3;

fft_analysis(t,s)
basicFFT(t,s)

%% Case 2

fs2 = 1024;
duration2 = 30;
N2 = fs2 * duration2;
t2 = 0:1/fs2:duration2-1/fs2;
f2 = 3*sin(25*2*pi*t2) + 0.4*sin(50*2*pi*t2 - 49*pi/180) + 5*(randn(size(t2)));
%f2 = 3*sin(25*2*pi*t2) + 0.4*sin(50*2*pi*t2 - 49*pi/180);

fft_analysis(t2,f2)
basicFFT(t2,f2)

%% Case 3

fs3 = 1024; % Sampling frequency
duration3 = 30; % Duration of signal
N3 = fs3 * duration3; % Number of samples
t3 = 0:1/fs3:duration3-1/fs3; % Time vector
f3 = 0.1*sqrt(2)*sin(250*2*pi*t3); % Vector respresenting 250 Hz sinusoidal signal with 0.1g RMS

fft_analysis(t3,f3)
basicFFT(t3,f3);

%% Case 4 - Rectified Sinusoid

% Reference - https://www.wolframalpha.com/input?i=fourier+series+abs%28sin%28x%29%29


fs4 = 2048; % Sampling frequency
duration4 = 32; % Duration of signal
N4 = fs4 * duration4; % Number of samples
t4 = 0:1/fs4:duration4-1/fs4; % Time vector
f4 = abs(sin(250*2*pi*t4)); % Vector respresenting 250 Hz rectified sinusoidal signal (has multiple frequency components)

fft_analysis(t4,f4)
basicFFT(t4,f4);

%% Problems - 

% Calibration when non-integer periods + effect of windowing ?
% Why dont we set N = total number of samples in the time window

% Intersting mappings (what conditions (frequency and amplitude) are good
% for experiments)
% Linearized system => Will give approximate 