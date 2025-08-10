%% FFT Analysis

% Here's a basic FFT function. It calculates as many averages as it can 
% without overlap using a specified number of FFT points, N and plots both 
% a single-sided magnitude spectrum and a level (in dB with a user specified level, 1 g by default).

% Suggestion – assign physical system units to your data of choice; check to make sure that consistent units are used for calculations.
% Check calculations for a known periodic function.

function [] = basicFFT(t,a,aref)

% DSP Parameters

  dt = mean(t(2:end)-t(1:end-1)); %(s) sample period

  NFFT = 1024; % fft points (powers of 2 preferred)
  %NFFT = length(t); % Can we set NFFT = length(t) for increased resolution ?

  df = 1/(NFFT*dt); %(Hz) freq resolution
  fv = (0:NFFT/2-1)*df; %(Hz) frequency vector for single-sided spectra
  Navg = floor(length(a)/NFFT); % # of no-overlap windows for averaging
  A = zeros(Navg,NFFT/2);

  % Reason for this ??
  for j = 1:Navg
    Aj = fft(a((j-1)*NFFT + (1:NFFT)),NFFT)/NFFT; %(g) double-sided amplitude spectrum
    A(j,:) = abs([Aj(1) 2*Aj(2:NFFT/2)]); %g) single-sided amplitude spectrum (magnitude)
  end

  Amag = mean(A,1);

  figure;
  semilogy(fv, Amag,'r-')
  xlabel('f (Hz)')
  ylabel('a (g)')
  xlim([fv(1) fv(end)])

  % Level calculations
  if nargin < 3
    aref = 1; %(g) default dB reference level
  end

  Alvl = 20*log10(Amag/aref);


  figure;
  plot(fv,Alvl,'b-','linewidth',2);
  xlabel('f (Hz)')
  ylabel(['A (dB, ref ', num2str(aref),' g)'])
  xlim([min(fv) max(fv)])

end

% Sample Data. Add random noise to the periodic function.
%t = 0:1/1024:30; %(S) time vector
%a = 3*sin(25*2*pi*t) + 0.4*sin(50*2*pi*t - 49*pi/180) + 5*(randn(size(t))); %(g) accel data