% Simultaneous Reconstruction of Quantized and Downsampled Audio Signals
% main script CVA2
% Vojtěch Kovanda
% Brno University of Technology, 2026


addpath('phase_correction');

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


%% settings for proposed algorithm (CVA-2)

w = 2048*4;
a = w/4;
M = w*2;

[win, ~] = generalizedCosWin(w, 'hann');
tight_win = calcCanonicalTightWindow(win, a);
tight_win = tight_win/norm(tight_win)*sqrt(a/w);
diff_win = numericalDiffWin(tight_win);
    
zeroPhaseFlag = true;
rotateFlag = true;

insig = x;

[sigIdx, sumIdx, sumArray, ifftArray, rotIdx] = precomputationForFDGT(length(insig), w, a, M);

% DGTs (original, its adjoint, and one with the differentiated window)
G = @(x) FDGT(x, tight_win, sigIdx, M, rotIdx, zeroPhaseFlag);
G_adj = @(u) invFDGT(u, tight_win, sumIdx, sumArray, ifftArray, rotIdx, zeroPhaseFlag)*w;
G_diff = @(x) FDGT(x, diff_win, sigIdx, M, rotIdx, zeroPhaseFlag);

% function to calculate the instantaneous frequency of the input signal
omega = @(x) calcInstFreq(G(x), G_diff(x), M, w, rotateFlag);

% operator to correct phase rotation and its adjoint
R = @(z, omega) instPhaseCorrection(z, omega, a, M);
R_adj = @(z, omega) invInstPhaseCorrection(z, omega, a, M);

% time-directional difference
D = @(z) z(:,1:end-1) - z(:,2:end);
D_adj = @(z) [z(:,1), (z(:,2:end) - z(:,1:end-1)), -z(:,end)];

% iPC-DGT
hatG = @(x, omega) D(R(G(x), omega));
hatG_adj = @(u, omega) G_adj(R_adj(D_adj(u), omega));

y_L = interp(y2, 2);
y_L = y_L(1:param.L)+0.0001*randn(param.L, 1);

    omega_y = omega(y_L);
    param.L1 = @(x) hatG(x, omega_y);
    param.L1_adj = @(u) hatG_adj(u, omega_y);
    param.x0 = y_L;


% algorithm parameters
param.lam = 0.0001;
param.rho = 0.8;
param.tau = 1;
param.sig = 1/2;

% maximal number of iteration
param.maxit = 100;

%% calling optimization algorithm



[xhat, SDR_t] = CVA2(y2, param, x);


%% evaluation

SDR = 20*log10(norm(x,2)./norm(x-xhat, 2));


fprintf('SDR of the reconstructed signal is %4.3f dB.\n', SDR);

% plot results

figure;
plot(SDR_t);
ylabel('SDR (dB)');
xlabel('number of iteration');