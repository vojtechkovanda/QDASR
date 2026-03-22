function [x, SDR] = CVA1(y1, param, in)
% CVA1 is the Condat-Vu algorithm
%
% Vojtěch Kovanda
% Brno University of Technology, 2026

lam = param.lam; % setting threshold for clipping

% definition of clip function (result of the Fenchel-Rockafellar conjugate of soft thresholding)
clip = @(x) (sign(x).*min(abs(x), lam));

%% initial values
i = 0; % number of iteration

x = interp(y1, 2)+0.0001*randn(param.L, 1);
u1 = zeros(param.L, 1);
u1 = frana(param.F, u1);
u1 = zeros(size(u1));
u2 = zeros(floor(param.L/param.k), 1);
SDR = zeros(param.maxit-1, 1);

%% algorithm
while i < param.maxit

    i = i + 1;
    waitbar(i/param.maxit);
     
     
     U1 = frsyn(param.F, u1);
     U1 = U1(1:param.L);
     U2 = zeros(param.L,1);
     for n = 1:param.k:param.k*length(u2)
         U2(n) = u2((n+param.k-1)/param.k);
     end


     x_tild = x - param.tau * (U2+U1);
     x = param.rho * x_tild + (1 - param.rho) * x;

     bL = 2*x_tild-x;

     p1 = u1 + param.sig * frana(param.F, bL);
     u1_tild = clip(p1);
     u1 = param.rho * u1_tild + (1 - param.rho) * u1;

     aL = bL(1:param.k:end);
     aL = aL(1:floor(param.L/param.k));

     p2 = u2 + param.sig * aL;
     u2_tild = p2 - param.sig * projection(p2/param.sig, y1, param.w);
     u2 = param.rho * u2_tild + (1 - param.rho) * u2;
     
    % SDR through iterations
     SDR(i) = 20*log10(norm(in,2)./norm(in-x, 2));

end


