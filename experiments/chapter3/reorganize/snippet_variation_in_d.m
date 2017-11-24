col = cbrSelector('name', 'Paired', 'length', 12);
col = col{1};

load('/mnt/data/thesis_data/simulations_varying_k_thesis/measures_total.mat');
load('/mnt/data/thesis_data/simulations_varying_k_thesis/nocomp_meas_var_k.mat');
load('/mnt/data/thesis_data/simulations_varying_k_thesis/SC_measures_var_k.mat');
load('/mnt/data/thesis_data/simulations/D_small.mat');
D = logspace(log10(25), log10(200), 10);
D = [D_small, D]; nb_d = numel(D);
K = round(logspace(log10(10), log10(100), 6)); nb_k = numel(K);


for fi = 1:1
    figure(100); hold off;
    data = squeeze(const_comp_meas_var_k(fi, :, :, 1:2));
    SC_data = SC_measures_var_k(fi, :, 2);
    % data= squeeze(measures(:, :, 3, 2, 1, 1, 1:2));
    c1 = [.75, 0, .25]; c2 = [.25, 0, .75];
    %subplot(1,2, 1);

%     for di=1:nb_d
%         alph = (nb_d - di)/(nb_d-1);
%         ci = alph*c1 + (1-alph)*c2;
%         scatter((data(di, :, 2)-SC_data)./SC_data, data(di, :, 1), 100, col(mod(di-1, 12)+1, :)/255, 'filled');
%         hold on;
%     end
    sub_data = (data(:, :, 2) - repmat(SC_data, nb_d, 1)) ./ repmat(SC_data, nb_d, 1);
    med = median(sub_data, 2);
    low = prctile(sub_data, 16, 2);
    hig = prctile(sub_data, 84, 2);
    errorbar(round(D), med, low, hig);

    % title('P = 0.25');
    legend(cellstr(num2str(round(D)', 'D=%-d')));
    % title('n = 70 %');
    title(sprintf('k = %d %', K(fi)));
    % xlim([4, 8]);
    % ylim([2.5, 4.3]);
end
%%
min_perc = 5;
max_perc = 95;

figure;
data = squeeze(measures_var_k(:, :, :, 1, 1, 2));

hold on;
for di=1:10
    scatter(vec(data(:, di, :)), vec(repmat(D_full(di), 4, 1, 100)), 50, col(di, :)/255, 'filled');
end

for ki=1:4
    data_sub_k = permute(squeeze(data(ki, :, :)), [2, 1]);
    med_err = median(data_sub_k, 1); low_err = prctile(data_sub_k, min_perc, 1); hig_err = prctile(data_sub_k, max_perc, 1);
    errorbarxy(med_err, D_full(1:10)+cos(pi/3*(ki-1)), med_err - low_err, hig_err - med_err, zeros(1, 10), zeros(1, 10), {'ko-', colors{ki}, colors{ki}});
end

hold off;
% title('P = 0.25');
legend(cellstr(num2str(round(D)', 'D=%-d')));
% title('n = 70 %');
% xlim([4, 8]);
% ylim([2.5, 4.3]);
%%

data = squeeze(measures(:, :, 3, 1, 1:2));
%subplot(1,2, 2);
hold on;
for di=1:nb_d
    alph = (nb_d - di)/(nb_d-1);
    ci = [0,0,1];%alph*c1 + (1-alpha)*c2;
    scatter(data(di, :, 2), data(di, :, 1), 100, ci, 'MarkerFaceColor', ci);
    alpha(0.1);
end
% xlim([4, 8]);
% ylim([2.5, 4.3]);
title('P = 0');
hold off;

%%

figure;
data = squeeze(measures(:, :, 1, 2, 4:5)) - squeeze(measures(:, :, 1, 1, 4:5));
%hold on;
for di=1:nb_d
    alph = (nb_d - di)/(nb_d-1);
    ci = alph*c1 + (1-alph)*c2;
    scatter(data(di, :, 2), data(di, :, 1), 100, ci, 'MarkerFaceColor', ci);
    drawnow;
    xlim([-1/di, 1/di]); ylim([-2, 2]);
    pause(2);
end

title('P = 0.25 - P = 0');
hold off;

%% 

data = squeeze(measures(:, :, 3, 1, 4:5));
data_x = data(:, :, 2);
data_y = data(:, :, 1);
idx_rem = find(data_x > 10);
data_x(idx_rem) = [];
data_y(idx_rem) = [];