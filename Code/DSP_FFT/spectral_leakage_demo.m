% Spectral Leakage Demonstration Code:
% % Change nT to an integer to eliminate spectral leakage. This would % %
% % suggest that the entire window could be periodic. A non-integer   % %
% % nT precludes that our entire dataset can be periodic, so we get   % %
% % spectral leakage.                                                 % %

nT = 6; % number of periods to simulate
fe = 5; %(Hz) excitation frequency
N = 2^11; % number of points in dataset
phi = 2*pi*rand; %(rad) phase shift
fs = (fe*N)/nT; %(Hz) sample freq
df = fs/N; %(Hz) frequency resolution, using 
t = (0:N-1)/fs; %(s) define sampled time vector
x = sin(2*pi*fe*t - phi); %(mm) sampled signal using dt*df*N = 1
X = fft(x,N)/N; %(mm) Fourier Transform in amplitude units
X = abs([X(1) 2*X(2:N/2)]); %(mm) make it a single-sided amplitude spectrum 
f = (0:N/2-1)*df; %(Hz) frequency vector
figure(2010)
subplot(2,1,1)
  plot(t,x,'r-','linewidth',2)
  xlabel('t (s)');
  ylabel('x (mm)');
  title('time domain)')
subplot(2,1,2)
  semilogy(f,X,'r-','linewidth',2)
  xlabel('f (Hz)')
  ylabel('X (mm)')
  xlim([0 20]) % zoom in on relevant frequency range
  title('frequency domain')
