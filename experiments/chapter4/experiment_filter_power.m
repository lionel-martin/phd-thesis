%EXPERIMENT_FILTER_POWER Script for the experiment on filter power (Figure 4.2)
%
% Added the comparison with Allen-Zhu approximation in my thesis.
% Modified by Lionel Martin
% Original work Johan Paratte and Lionel Martin
%
% Copyright (C) 2016 Johan Paratte, Lionel Martin.
% This file is part of the Reproducible Research code base for the Fast
% Eigenspace Approximation using Random Signals (FEARS) method.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% If you use this code please kindly cite
%     Paratte, Johan, and Lionel Martin. 
%     "Fast Eigenspace Approximation using Random Signals."
%     arXiv preprint arXiv:1611.00938 (2016).
% https://arxiv.org/abs/1611.00938

%% Study the effect of the polynomial power on the approximation

N = 500;
k = 150;

G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);

lk = (G.e(k) + G.e(k+1)) / 2;
h = @(x) x < lk; % ideal low-pass

%% Cheby vs Jackson-Cheby
m = 250;

figure; hold on;
paramplot.plot_eigenvalues = 1;
paramplot.cla = 0;
gsp_plot_filter(G, h, paramplot);


paramplot.plot_eigenvalues = 0; %avoid redrawing eigenvalues

cheby_approx = gsp_approx_filter(G, h, m);
gsp_plot_filter(G, cheby_approx, paramplot);

jch_approx = approx_filter_jch(G, lk, m);
gsp_plot_filter(G, jch_approx, paramplot);

az_param = struct('order', m, 'kap', 1e-4);
az_approx = approx_filter_allenzhu(G, lk, m, az_param);
gsp_plot_filter(G, az_approx, paramplot);

hold off;


%% Jackson-Cheby increasing m
ms = [50, 100, 200, 400];

figure; hold on;
paramplot.plot_eigenvalues = 1;
paramplot.cla = 0;

gsp_plot_filter(G, h, paramplot);
paramplot.plot_eigenvalues = 0;

for m = ms
    jch_approx = approx_filter_jch(G, lk, m);
    gsp_plot_filter(G, jch_approx, paramplot);

end

%% Allen-Zhu increasing m
ms = [50, 100, 200, 400];
kaps = [1e-2, 1e-3, 1e-3, 1e-4];

figure; hold on;
paramplot.plot_eigenvalues = 1;
paramplot.cla = 0;

gsp_plot_filter(G, h, paramplot);
paramplot.plot_eigenvalues = 0;

for ii = 1:numel(ms)
    m = ms(ii);
    az_param = struct('order', m, 'kap', kaps(ii));
    az_approx = approx_filter_allenzhu(G, lk, m, az_param);
    gsp_plot_filter(G, az_approx, paramplot);
end