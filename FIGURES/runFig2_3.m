clear;
close all;
clc;

addpath('..\FUNCTIONS\');

data = load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\expDataAnalysis.mat');
params = load('..\PARAMETER_DATA\params.mat');

figInfo = figure; hold on;
figVar = figure; hold on;
figTuningPop = figure;
figTuningBP = figure;

color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];...
    [0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840];...
    [0.5 0 0];[0 0.5 0];[0 0 0.5]];


N = numel(params.rat);
for i = 1:N
    ratNumb = params.rat{i};
    
    % population state phase (range: [-1 1])
%     [pc1, pc2] = alignLatentSlide(data.iso2{i}, data.iso1{i},2);
%     popStatPhase = mean(atan2(pc2, pc1),2)/pi;
    popStatPhase = mean(data.unalgnRotPh{i},2)/pi;
    popStatPhase = round(popStatPhase*12);
    
    % blood pressure (range: [-1 1])
    bp = mean(data.bp{i},2);
    bp = round(bp*12);

    t = data.t{i};

%     for j = 1:numel(data.fr{i})
%         fr(:,j) = mean(data.fr{i}(j},2);
%     end
    fr = data.fr{i}';

    [~, tunPop] = get_tuning_curve(popStatPhase, fr);
    normTunPop = norm_tuning(tunPop);
    % Sort unsupervised tunning curve
    [~, argmax_tuning] = max(normTunPop, [], 1);
    [~, srt_idx_Pop] = sort(argmax_tuning);
    [~, tunBP] = get_tuning_curve(bp, fr);
    normTunBP = norm_tuning(tunBP);
    % Sort unsupervised tunning curve
    [~, argmax_tuning] = max(normTunBP, [], 1);
    [~, srt_idx_BP] = sort(argmax_tuning);
    if i == 1
        figure(figTuningBP);
        imagesc(normTunBP(:,srt_idx_BP)');
        colormap('parula');
        figure(figTuningPop);
        imagesc(normTunPop(:,srt_idx_Pop)');
        colormap('parula');
    end

    [~, infoPop] = get_information(popStatPhase, fr);
    [~, infoBP] = get_information(bp, fr);
    figure(figInfo);
    scatter(infoBP, infoPop, 50, color(i,:), 'filled');
    avgInfoPop(i) = mean(infoPop);
    avgInfoBP(i) = mean(infoBP);

    fovBP = exp_var_tuning(tunBP(bp-min(bp)+1,:), fr);
    fovPop = exp_var_tuning(tunPop(popStatPhase-min(popStatPhase)+1,:), fr);
    figure(figVar);
    scatter(fovBP, fovPop, 50, color(i,:), 'filled');
    avgVarPop(i) = mean(fovPop);
    avgVarBP(i) = mean(fovBP);

%     figure(fig);
%     for j = 1:size(fr,2)
%         tmp = corrcoef(bp, fr(:,j));
%         test(j) = tmp(1,2)^2;
%         plot(bp, fr(:,j));
%     end
%     Test(i) = mean(test);
%     clearvars tmp test;
end



figure(figInfo);
plot((0:6), (0:6), 'k', 'LineStyle','--');
figure(figVar);
plot((0:1), (0:1), 'k', 'LineStyle','--');

figure;
ax = axes();
hold(ax);
boxchart(ones(size(avgInfoBP)), avgInfoBP,...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(avgInfoPop)), avgInfoPop,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
plot(repmat([1;2],1,N), [avgInfoBP;avgInfoPop], 'k');
xlim([0.5 2.5]);ylim([0 1.2]);
scatter(ones(N,1), avgInfoBP, 50, 'white', 'filled', 'MarkerEdgeColor','k');
scatter(2*ones(N,1), avgInfoPop, 50, 'white', 'filled', 'MarkerEdgeColor','k');
disp('Information');
ranksum(avgInfoBP, avgInfoPop)

figure;
ax = axes();
hold(ax);
boxchart(ones(size(avgVarBP)), avgVarBP,...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(avgVarPop)), avgVarPop,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
plot(repmat([1;2],1,N), [avgVarBP;avgVarPop], 'k');
xlim([0.5 2.5]);ylim([0 1]);
scatter(ones(N,1), avgVarBP, 50, 'white', 'filled', 'MarkerEdgeColor','k');
scatter(2*ones(N,1), avgVarPop, 50, 'white', 'filled', 'MarkerEdgeColor','k');
disp('Explained Var. by Tuning curve');
ranksum(avgVarBP, avgVarPop)