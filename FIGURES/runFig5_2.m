clc;
clear;
close all;

addpath('..\FUNCTIONS\');
rat = {...
    '19', '21', '22',...
    'B2', 'B3', 'B4', 'B9_a1',...
    'B10_a1', 'B14_a1', 'B15_a1'};

load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\hBINDTestKFold.mat');
load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\hBINDTestKFoldFeedforwardStat.mat');

c1 = [42 183 202]/255;
c2 = [0.9290 0.6940 0.1250];
c3 = [254 74 73]/255;%[0.7 0.7 0.7];

%% BP prediction accuracy
figure;
ax = axes();
hold(ax);
boxchart(ones(size(corND)), corND,...
    'BoxFaceColor', c2, 'BoxFaceAlpha',1, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(corND)), corhBIND,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',1, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [corND, corhBIND]', 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1), corND, 100,...
    'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor','none');
scatter(2*ones(10, 1), corhBIND, 100,...
    'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor','none');
ylim([0.6 1]);
xlim([0.5 2.5]);
yline(0.9, 'LineWidth',2,'LineStyle','--', 'Color','k');

figure; hold on;
scatter(corND, corhBIND, 50, c1, 'filled');
plot((0.6:0.1:1), (0.6:0.1:1), 'k', 'LineWidth', 1.5);
disp(['ND vs. hBIND: BP prediction accuracy: p = ', num2str(signrank(corND, corhBIND))]);

%% 
dE = 0.025;
edges = (-1.5:dE:1.5);
range = (-20:0.5:80);
L = numel(range);
e_ND = 10*ones(L, numel(rat));
e_hBIND = 10*ones(L, numel(rat));
eP_ND = 10*ones(L, numel(rat));
eP_hBIND = 10*ones(L, numel(rat));
eR_ND = 10*ones(L, numel(rat));
eR_hBIND = 10*ones(L, numel(rat));
for i = 1:numel(rat)
    tTmp = tt{i,1};
    
    idx = zeros(size(range));
    for k = 1:numel(range)
        [~, idx(k)] = min(abs(tTmp-range(k)));
    end
    
    t = tTmp(idx);
    tr = numel(errorND(i,:));
    for k = 1:numel(errorND(i,:))
        if isempty(errorND{i,k})
            tr = tr-1;
        end
    end
    errTmpND = zeros(L, tr);
    errTmphBIND = zeros(L, tr);
    errPhaseTmpND = zeros(L, tr);
    errPhaseTmphBIND = zeros(L, tr);
    errRadiTmpND = zeros(L, tr);
    errRadiTmphBIND = zeros(L, tr);
    for j = 1:tr
        errTmpND(:,j) = errorND{i,j}(idx);
        errTmphBIND(:,j) = errorhBIND{i,j}(idx);
        errPhaseTmpND(:,j) = errorPhaseND{i,j}(idx);
        errPhaseTmphBIND(:,j) = errorPhasehBIND{i,j}(idx);
        errRadiTmpND(:,j) = errorRadiND{i,j}(idx);
        errRadiTmphBIND(:,j) = errorRadihBIND{i,j}(idx);
    end
    e_ND(:,i) = mean(errTmpND.^2,2); 
    e_hBIND(:,i) = mean(errTmphBIND.^2,2);
    eP_ND(:,i) = mean(errPhaseTmpND.^2,2);
    eP_hBIND(:,i) = mean(errPhaseTmphBIND.^2,2);
    eR_ND(:,i) = mean(errRadiTmpND.^2,2);
    eR_hBIND(:,i) = mean(errRadiTmphBIND.^2,2);
end

%% BP prediction error distribution
figE = figure; hold on;
d = 0.01; edges = (-0.01:d:1.5);
drawCumErrFrac(figE, edges, e_ND, c2);
drawCumErrFrac(figE, edges, e_hBIND, [1 0 0]);
xlim([edges(1) edges(end)]);
ylim([-0.01 1.01]);

f = figure; f.Position(3) = 260;
    dat = e_ND;
dat1 = sum(dat,1)/(size(dat,1)-2);
    dat = e_hBIND;
dat2 = sum(dat,1)/(size(dat,1)-2);
disp(['ND vs. hBIND: BP prediction error: p = ', num2str(signrank(dat1, dat2))]);
ax = axes();
hold(ax);
boxchart(ones(size(dat1)), dat1,...
    'BoxFaceColor', c2, 'BoxFaceAlpha',1, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(dat2)), dat2,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',1, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [dat1 ; dat2], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1), dat1, 100,...
    'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor','none');
scatter(2*ones(10, 1), dat2, 100,...
    'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor','none');
ylim([0 0.3]);
xlim([0.5 2.5]);

%% Population state phase error distribution
figPhase = figure; hold on; figPhase.Position(3) = 480;
d = 0.01; edges = (-d:d:1);
drawCumErrFrac(figPhase, edges, eP_ND/pi, c2);
drawCumErrFrac(figPhase, edges, eP_hBIND/pi, [1 0 0]);
xlim([edges(1) edges(end)]);
ylim([-0.01 1.01]);

%
f = figure; f.Position(3) = 260;
    dat = eP_ND/pi;
dat1 = sum(dat,1)/(size(dat,1)-2);
    dat = eP_hBIND/pi;
dat2 = sum(dat,1)/(size(dat,1)-2);
disp(['ND vs. hBIND: phase error: p = ', num2str(signrank(dat1, dat2))]);
ax = axes();
hold(ax);
boxchart(ones(size(dat1)), dat1,...
    'BoxFaceColor', c2, 'BoxFaceAlpha',1, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(dat2)), dat2,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',1, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [dat1 ; dat2], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1), dat1, 100,...
    'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor','none');
scatter(2*ones(10, 1), dat2, 100,...
    'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor','none');
ylim([0 0.3]);
xlim([0.5 2.5]);

%% Population state radial error distribution
figRadi = figure; hold on; figRadi.Position(3) = 480;
d = 0.001; edges = (-d:d:0.15);
drawCumErrFrac(figRadi, edges, eR_ND,c2);
drawCumErrFrac(figRadi, edges, eR_hBIND,[1 0 0]);
xlim([edges(1) edges(end)]);
ylim([-d 1.01]);

%
f = figure; f.Position(3) = 260;
    dat = eR_ND;
dat1 = sum(dat,1)/(size(dat,1)-2);
    dat = eR_hBIND;
dat2 = sum(dat,1)/(size(dat,1)-2);
disp(['ND vs. hBIND: radial error: p = ', num2str(signrank(dat1, dat2))]);
ax = axes();
hold(ax);
boxchart(ones(size(dat1)), dat1,...
    'BoxFaceColor', c2, 'BoxFaceAlpha',1, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(dat2)), dat2,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',1, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [dat1 ; dat2], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1), dat1, 100,...
    'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor','none');
scatter(2*ones(10, 1), dat2, 100,...
    'MarkerFaceColor','k', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeColor','none');
ylim([0 0.03]);
xlim([0.5 2.5]);