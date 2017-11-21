%   Experiment of Table 4.2.
%   Comparison of the approximation quality of Jackson-Cheby and Allen-Zhu.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
N = 500;
k = 180;

sensorparams.connected = 1;
G = gsp_sensor(N, sensorparams);
G = gsp_compute_fourier_basis(G);

C = [50, 100, 200, 400]; nb_c = numel(C);
KAP = [1e-2, 1e-3, 1e-4, 1e-5];

for ii = 1:nb_c
    c = C(ii);
    fprintf('Computations with polynomial order c=%d\n', c);
    exact_filt = [ones(k, 1); zeros(N-k, 1)];
    
    jch_approx = approx_filter_jch(G, (G.e(k) + G.e(k+1))/2, c);
    azparam.kap = KAP(ii);
    az_approx = approx_filter_allenzhu(G, (G.e(k) + G.e(k+1))/2, c, azparam);
    
    fprintf('Jackson-Cheby approx: %f\n', norm(jch_approx(G.e) - exact_filt));
    fprintf('Allen-Zhu approx: %f\n', norm(az_approx(G.e) - exact_filt));
end