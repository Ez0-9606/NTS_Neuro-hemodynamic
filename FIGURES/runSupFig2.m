clear;
close all;
clc;

addpath('..\FUNCTIONS\');
params = load('..\PARAMETER_DATA\params.mat');

%%
figure(1); hold on;
figure(2); hold on;
color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];...
    [0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840];...
    [0.5 0 0];[0 0.5 0];[0 0 0.5]];
for i = 1:numel(params.rat)
    ratNumb = params.rat{i};
    data = load(['..\POSTPROCESSED_MAT_CODE\DATA\', ratNumb, '_matlab_save.mat']);
    h1 = data.h1;
    h1_len = h1(:,2) - h1(:,1);
%     h1_len = sort(h1_len/max(h1_len), 'descend');
%     boolIdx = h1_len > prctile(h1_len, 50);
%     h1_len = h1_len(boolIdx);

    h2 = data.h2;
    h2_len = h2(:,2) - h2(:,1);
%     h2_len = sort(h2_len/max(h2_len), "descend");

    snr1(i) = (max(h1_len)-mean(h1_len))/std(h1_len);
    snr2(i) = (max(h2_len)-mean(h2_len))/std(h2_len);
    
    figure(1);
    %plot((1:numel(h1_len)), h1_len, 'Color', color(i,:), 'LineWidth',1);
    scatter(ones(size(h1_len))*i+0.2*(rand(size(h1_len))-0.5), h1_len, 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', color(i,:));
%     yline(mean(h1_len) + std(h1_len));
%     figure(10);
%     scatter(1+rand(1)-0.5, numel(h1_len), 50,'MarkerEdgeColor', 'none', 'MarkerFaceColor', color(i,:));

    figure(2);
    scatter(ones(size(h2_len))*i+0.2*(rand(size(h2_len))-0.5), h2_len, 50, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', color(i,:));
%     scatter(max(h1_len), max(h2_len), 50, color(i,:), 'filled');

end
figure(1);xlim([0 11]);
figure(2);xlim([0 11]);
figure;
ax = axes();
hold(ax);
boxchart(ones(size(snr1)), snr1,...
    'BoxFaceColor', 'b', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(snr2)), snr2,...
    'BoxFaceColor', 'g', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [snr1; snr2], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1)+0.2*(rand(10,1)-0.5), snr1, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(10, 1)+0.2*(rand(10,1)-0.5), snr2, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');

%%
figure(4); hold on;
    x = (-5:0.01:5);
    for j = 1:30
        plot(x, x-4+0.25*j,'k');
    end

figure(5); hold on;
for i = 1:numel(params.rat)
    ratNumb = params.rat{i};
    data = load(['..\POSTPROCESSED_MAT_CODE\DATA\', ratNumb, '_matlab_save.mat']);

    dt = params.DT(i);
    tr = params.TRIAL(i);
    bDelay = params.bpDELAY(i);
    nDelay = params.ntsDELAY(i);
    pre = params.PRE(i);

    [t, bp, fr, iso1, iso2] = dataPrepare(data, dt, tr, nDelay, bDelay, pre);

    figure(4);
    scatter(data.th(7:end-7), data.corr_integral(7:end-7), 50, 'filled');
    xlim([-0.25 2.5]);
    ylim([-1 2]);
    
    fitting = fitlm(data.th(7:end-7), data.corr_integral(7:end-7));
    slope(i) = table2array(fitting.Coefficients(2,1));
    %sem = table2array(fitting.Coefficients(2,2));

end
f = figure(5);f.Position(3) = 280; hold on;
bar(1, mean(slope), 'BarWidth', 0.5, 'FaceColor', 0.5*[1 1 1]);
errorbar(1, mean(slope), 0, std(slope), 'LineWidth', 2, 'Color', 'k');
scatter(ones(size(slope))+0.3*(rand(size(slope))-0.5), slope, 50, 'w', 'filled', 'MarkerEdgeColor','k');
xlim([0.5 1.5]);
ylim([0 1.5]);

%%
figure(6); hold on;
exp_var_ratio = zeros(6, numel(params.rat));
cum_exp_var_r = zeros(1, numel(params.rat));
for i = 1:numel(params.rat)
    ratNumb = params.rat{i};
    data = load(['..\POSTPROCESSED_MAT_CODE\DATA\', ratNumb, '_matlab_save.mat']);
    exp_var = var(data.fr_isomap, [], 1);
    exp_var_ratio(:,i) = exp_var./sum(exp_var);
    scatter((1:6)+0.3*(rand(1,6)-0.5), exp_var_ratio(:,i), 100,...
        'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor','none');

    cum_exp_var_r(i) = sum(exp_var_ratio(1:2,i));
end
c = 'r';
errorbar((1:6), mean(exp_var_ratio,2), std(exp_var_ratio,[],2), std(exp_var_ratio,[],2),...
    '-o', 'LineWidth', 3, 'Color', 'k', 'MarkerSize', 5);
% fig.CurrentAxes.YAxis(1).Color = c;
xlim([0.5 6.5]);
ylim([0 0.8]);

f = figure(7);f.Position(3) = 280; hold on;
bar(1, mean(cum_exp_var_r), 'BarWidth', 0.5, 'FaceColor', [189 117 230]/255, 'FaceAlpha',0.5);
errorbar(1, mean(cum_exp_var_r), 0, std(cum_exp_var_r), 'LineWidth', 2, 'Color', 'k');
scatter(ones(size(cum_exp_var_r))+0.3*(rand(size(cum_exp_var_r))-0.5), cum_exp_var_r, 50, 'w', 'filled', 'MarkerEdgeColor','k');
xlim([0.5 1.5]);