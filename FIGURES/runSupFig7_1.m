clear;
close all;
clc;

addpath('..\FUNCTIONS\');

params = load('..\PARAMETER_DATA\params.mat');

orD = load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\expDataAnalysis.mat');

iso1 = orD.algn1{1};
iso2 = orD.algn2{1};
bp = orD.bp{1};
tst = orD.fr{1};
disp(orD.accuComm(1));
corr = corrcoef(mean(iso2,2), mean(bp,2));
disp(corr(1,2)^2);

accuAllCell = zeros(1,10);
for i = 1:numel(params.rat)
    fr = orD.fr{i};
    t = orD.t{i};
    tL = numel(t);
    [~, ~, allCellCorr, ~] = linRegress(fr', mean(bp,2));
    bpest(i) = allCellCorr^2;
    [reg, y, allCellCorr, ~] = linRegress(fr', mean(iso1,2));
    iso1est(i) = allCellCorr^2;
    [reg, y, allCellCorr, ~] = linRegress(fr', mean(iso2,2));
    iso2est(i) = allCellCorr^2;
    
end
figure;


accuLatent = orD.accuComm;
figure;
ax = axes();
hold(ax);
boxchart(ones(size(accuLatent)), accuLatent,...
    'BoxFaceColor', 'b', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(bpest)), bpest,...
    'BoxFaceColor', 'g', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
plot(repmat([1;2],1,10), [accuLatent;bpest], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1)+0.2*(rand(10,1)-0.5), accuLatent, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(10, 1)+0.2*(rand(10,1)-0.5), bpest, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
ylim([0 1]);
xlim([0.5 2.5]);
disp(ranksum(accuLatent, bpest));

%%
figure;
ax = axes();
hold(ax);
boxchart(ones(size(iso1est)), iso1est,...
    'BoxFaceColor', 'c', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(iso2est)), iso2est,...
    'BoxFaceColor', 'm', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
plot(repmat([1;2],1,10), [iso1est;iso2est], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1)+0.2*(rand(10,1)-0.5), iso1est, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(10, 1)+0.2*(rand(10,1)-0.5), iso2est, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
ylim([0 1]);
xlim([0.5 2.5]);
disp(ranksum(iso1est, iso2est));

%%
tst = load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\hBINDTestKFoldFeedforwardRepresentative1.mat');
