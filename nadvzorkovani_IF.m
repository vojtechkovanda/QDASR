% nadvzorkovani

%function [SDR2, ODG2] = nadvzorkovani_IF(w, audiofile)

%          x
%       ___|
%      |         
%      B    
%     D_k     
%    Q_fine   
%     y_1
%
% B is anti-aliasing filter, D_k is downsampling, Q_fine is fine
% quantization
%
% the PEMO-Q audioqual is now inactive
%
% Vojtěch Kovanda
% Brno University of Technology, 2024


% using LTFAT toolbox
%ltfatstart


%% input signal
 % 
  audiofile = '2.wav';
 % % 
  w = 7;

[x, param.fs] = audioread(audiofile);

% signal length


% normalization
maxval = max(abs(x));
temp = mod(length(x), 2);
x = x(1:end-temp);
x = x/maxval;

param.L = length(x);

%% generate observations y

% setting conversion parameters
param.w = w;           % bit depth (bps) of Q_fine
param.k = 2;           % downsampling factor

% load impulse response of B for downsampling factor k = 4 and sampling
% frequency f_s = 48kHz or 44.1kHz
% load("filter_coeffs_44_1.mat");
% load("filter_coeffs.mat");
% param.B = Num;
% param.Bt = flip(param.B);


% filtering (using convolution)
%y = conv(x, param.B);


% signal length after filtering
%param.L1 = length(y);

% quantization
%y = quant(y, param.w);

y2 = quant(x, param.w);

% downsampling
%y = y(1:param.k:end);
%y = y(1:floor(param.L1/param.k));

y2 = y2(1:param.k:end);
y2 = y2(1:floor(param.L/param.k));

%xq = quant(x, param.w/param.k);

%% settings for proposed algorithm (CVA)


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
param.lam = 0.0001; %[0.0012 0.0012 0.0012 0.0012 0.0012 0.0012 0.0012 0.0012 0.0001 0.00005 0.00002 0.00001 0.000005 0.000001 0.0000005 0.0000001]; % different clipping thresholds for different bit depths of y_2
param.rho = 0.8;
param.tau = 1;
param.sig = 1/2;

% maximal number of iteration
param.maxit = 200;

param.outloop = 1;

%% calling optimization algorithm

%[xhat, SDR_t] = cv_alg_if(y, param, x);

for i = 1:param.outloop

[xhat2, SDR_t2] = without_B_IF(y2, param, x);

    omega_y = omega(xhat2);
    param.L1 = @(x) hatG(x, omega_y);
    param.L1_adj = @(u) hatG_adj(u, omega_y);
    param.x0 = xhat2;

end

%% evaluation

% SDR of reconstructed signal, SDR(xhat, x)
%[SDR, bestit] = max(SDR_t);



SDR2 = 20*log10(norm(x,2)./norm(x-xhat2, 2));

% ODG of reconstructed signal, ODG(x, x)
%[~, ~, ODG] = audioqual(x, xhat, param.fs); 

[~, ~, ODG2] = audioqual(x, xhat2, param.fs); 

% fprintf('SDR of the reconstructed signal is %4.3f dB.\n', SDR);
% fprintf('ODG of the reconstructed signal is %4.3f.\n', ODG);

fprintf('SDR2 of the reconstructed signal is %4.3f dB.\n', SDR2);
fprintf('ODG2 of the reconstructed signal is %4.3f.\n', ODG2);

% plot results

% figure;
% plot(SDR_t);
% ylabel('SDR (dB)');
% xlabel('number of iteration');
% 
% y_zero = downsample(x, 2);
% y_zero = upsample(y_zero, 2);
% y_zero = quant(y_zero, param.w);
% 
% y_zero_anti = conv(x, param.B);
% y_zero_anti = quant(y_zero_anti, 8);


%figure;
% plot(SDR_t2);
% ylabel('SDR (dB)');
% xlabel('number of iteration');



y_L = interp(y2, 2);
y_L = y_L(1:param.L);%+0.0001*randn(param.L, 1);

% figure;
% subplot(2, 3, 1);
% sgram(x, param.fs);
% clim([-100 10]);
% subplot(2, 3, 2);
% sgram(quant(x, param.w), param.fs);
% clim([-100 10]);
% subplot(2, 3, 3);
% sgram(y_zero, param.fs);
% clim([-100 10]);
% subplot(2, 3, 4);
% sgram(y_zero_anti, param.fs);
% clim([-100 10]);
% subplot(2, 3, 5);
% sgram(y_L, param.fs);
% clim([-100 10]);
% subplot(2, 3, 6);
% sgram(xhat2, param.fs);
% clim([-100 10]);

[~, ~, ODG_int] = audioqual(x, y_L, param.fs); 
SDR_int = 20*log10(norm(x,2)./norm(x-y_L, 2));

fprintf('SDR_int is %4.3f dB.\n', SDR_int);
fprintf('ODG_int is %4.3f.\n', ODG_int);