%   Experiment of Figure 5.2.
%   Replacing the assignment of k-medoids with that of the diffusion of
%   Kronecker positioned at the centers of k-medoids.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin
N = 500;
k = 3;

G = gsp_sensor(N);
G = gsp_create_laplacian(G, 'normalized');
[Uk, ek] = compute_eigen(G, k);

[assignment, ~, ~, ~, centers] = kmedoids(Uk, k, 'Replicates', 50);

Uc = sum(Uk(centers, :).^2, 2);
Uc = repmat(Uc', G.N, 1);

cent_sigs = zeros(G.N, k);
cent_sigs(sub2ind([G.N, k], centers', 1:k)) = 1;
dist_apx = Uk*Uk'*cent_sigs; % exact filtering with the ideal low-pass filter
[~, ass_apx] = max(dist_apx, [], 2);

dist_ex = dist_apx - Uc/2;

%% Plotting

[ord_class, ord_nod_idx] = sort(assignment);

figure(1); hold off; hold on; set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 500 400]);
xlabel('node index'); ylabel('objective function');
for ii = 1:k
    scatter(1:N, dist_ex(ord_nod_idx, ii), 'filled');
end

figure(2); hold off; hold on; set(gcf, 'Color', 'w');
set(gcf, 'Position', [100 100 500 400]);
xlabel('node index'); ylabel('objective function');
for ii = 1:k
    scatter(1:N, dist_apx(ord_nod_idx, ii), 'filled');
end

plotparam.vertex_highlight = centers;
figure(3); gsp_plot_signal(G, assignment, plotparam);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [100 100 500 400]); colorbar('off');
figure(4); gsp_plot_signal(G, ass_apx, plotparam);
set(gcf, 'Color', 'w'); set(gcf, 'Position', [100 100 500 400]); colorbar('off');