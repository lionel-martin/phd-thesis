%   Experiment of Figure 4.1.
%   Quality of the eigenvectors recovery using FEARS with the ideal lowpass.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin

%% Top row (sensor)
N = 10000;
k = 2;

for gid = 1:2
    switch gid
        case 1
            paramsensor.connected = 1;
            G = gsp_sensor(N, paramsensor);
        case 2
            gparams = struct('deg', 100, 'eps_factor', 0.25);
            G = generate_sbm(N, k, gparams);
    end
    [Uk, ek] = compute_eigen(G, k);
    figure; gsp_plot_signal(G, Uk(:, 2));

    fearsparam.filter = 'exact';
    [Bk, Psi] = gsp_eigenspace_estimation(G, k, fearsparam);
    figure; gsp_plot_signal(G, Bk(:, 2));
    figure; gsp_plot_signal(G, Psi(:, 2));
end

for coordid = 1:3
    switch coordid
        case 1
            G.coords = Uk;
        case 2
            G.coords = Bk;
        case 3
            G.coords = Psi;
    end
    figure;
    gsp_plot_signal(G, G.info.node_com);
end