% Simultaneous Reconstruction of Quantized and Downsampled Audio Signals
% main script CVA1
% Vojtěch Kovanda
% Brno University of Technology, 2026


% using LTFAT toolbox
ltfatstart

%% input signal
 audiofile = 'dataset/1.wav';
[x, param.fs] = audioread(audiofile);

% normalization
maxval = max(abs(x));
temp = mod(length(x), 2);
x = x(1:end-temp);
x = x/maxval;

% signal length
param.L = length(x);

%% generate observations y

% setting conversion parameters
param.w = 6;           % bit depth (bps)
param.k = 2;           % downsampling factor

% quantization
y2 = quant(x, param.w);

% downsampling
y2 = y2(1:param.k:end);
y2 = y2(1:floor(param.L/param.k));

%% settings for proposed algorithm (CVA)

% frame settings
param.winlen = 2048*4;            % window length
param.wtype = 'hann';           % window type
param.a = param.winlen/4;       % window shift
param.M = 2*param.winlen;       % number of frequency channels

% frame construction
param.F = frametight(frame('dgtreal', {param.wtype, param.winlen}, param.a, param.M));
param.F = frameaccel(param.F, param.L);

% algorithm parameters
param.lam = 0.0001;
param.tau = 1;
param.sig = 1/2;
param.rho = 1;

% maximal number of iteration
param.maxit = 100;

%% calling optimization algorithm

[xhat, SDR_t] = CVA1(y2, param, x);

%% evaluation

SDR = 20*log10(norm(x,2)./norm(x-xhat, 2));


%fprintf('SDR of the reconstructed signal is %4.3f dB.\n', SDR);
%fprintf('ODG of the reconstructed signal is %4.3f.\n', ODG);

fprintf('SDR of the reconstructed signal is %4.3f dB.\n', SDR);

% plot results

%figure;
plot(SDR_t);
ylabel('SDR (dB)');
xlabel('number of iteration');