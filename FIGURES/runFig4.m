clear;
close all;
clc;

addpath('..\FUNCTIONS\');
load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\expDataAnalysis.mat');

%%
color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];...
    [0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840];...
    [0.5 0 0];[0 0.5 0];[0 0 0.5]];

c1 = [42 183 202]/255;
c2 = [0.9290 0.6940 0.1250];
c3 = [254 74 73]/255;%[0.7 0.7 0.7];

%%
fig1 = figure; hold on;
fig2 = figure; hold on;
fig3 = figure; hold on;%xline(0); 
fig4 = figure; hold on;
fig5 = figure; hold on;
fig6 = figure; hold on;

fig = figure; hold on;

for ii = 1:10
    i = 11-ii;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure(fig);
%     % blood pressure (range: [-1 1])
%     bpTest = mean(bp{i},2);
%     bpTest = round(bpTest*12);
%     popAlgnPh = mean(algnRotPh{i},2)/pi;
%     popAlgnPh = popAlgnPh - mean(popAlgnPh(1:40));
%     stateBase = (-0.1:0.1:2);
%     for j = 1:numel(popAlgnPh)
%         [~, idx] = min(abs(stateBase - popAlgnPh(j)));
%         popAlgnPh(j) = stateBase(idx);
%     end
%     stateData = round(10*popAlgnPh+1);
%     [state, tuning] =  get_tuning_curve(stateData, bpTest);
%     [~, mnIdx] = min(tuning);
%     state = state - state(mnIdx);
%     plot(state, tuning/12);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Unaligned neural trajectory
    figure(fig1);
    scatter(iso1{i}(:), iso2{i}(:), 30, color(i,:));
    plot(iso1{i}(:), iso2{i}(:), 'color', color(i,:));

    % Unaligned dynamics of population state
    figure(fig2);
    plot3(iso1{i}(1:sidx(i),:), iso2{i}(1:sidx(i),:), repmat(t{i}(1:sidx(i)),1,size(raln1{i},2)), 'color', c1, 'LineWidth',0.5);hold on;
    plot3(iso1{i}(sidx(i):eidx(i),:), iso2{i}(sidx(i):eidx(i),:), repmat(t{i}(sidx(i):eidx(i),:),1,size(raln1{i},2)), 'color', c2, 'LineWidth',0.5);
    plot3(iso1{i}(eidx(i):end,:), iso2{i}(eidx(i):end,:), repmat(t{i}(eidx(i):end,:),1,size(raln1{i},2)), 'color', c3, 'LineWidth',0.5);
    scatter3(iso1{i}(1:sidx(i),:), iso2{i}(1:sidx(i),:), repmat(t{i}(1:sidx(i),:),1,size(raln1{i},2)), 30, c1);hold on;
    scatter3(iso1{i}(sidx(i):eidx(i),:), iso2{i}(sidx(i):eidx(i),:), repmat(t{i}(sidx(i):eidx(i),:),1,size(raln1{i},2)), 30, c2);
    scatter3(iso1{i}(eidx(i):end,:), iso2{i}(eidx(i):end,:), repmat(t{i}(eidx(i):end,:),1,size(raln1{i},2)), 30, c3);
    
    % Unaligned hemodynamic coupling
    figure(fig3);
    scatter(interp(iso1{i}(:),60), interp(iso2{i}(:),60), 1, interp(bp{i}(:),60), 'filled');hold on;
    scatter(iso1{i}(:), iso2{i}(:), 30, bp{i}(:));hold on;
    set(fig3, 'Colormap', colormap(viridis(300)));
    clim([-1 1]);
    axis equal;
    axis off;

    % Aligned neural trajectory
    figure(fig4);
    scatter(raln1{i}(:), raln2{i}(:), 30, color(i,:));
    plot(raln1{i}(:), raln2{i}(:), 'color', color(i,:));
    axis equal;

    % Aligned dynamics of population state
    figure(fig5);
    plot3(raln1{i}(1:sidx(i),:), raln2{i}(1:sidx(i),:), repmat(t{i}(1:sidx(i)),1,size(raln1{i},2)), 'color', c1, 'LineWidth',0.5);hold on;
    plot3(raln1{i}(sidx(i):eidx(i),:), raln2{i}(sidx(i):eidx(i),:), repmat(t{i}(sidx(i):eidx(i),:),1,size(raln1{i},2)), 'color', c2, 'LineWidth',0.5);
    plot3(raln1{i}(eidx(i):end,:), raln2{i}(eidx(i):end,:), repmat(t{i}(eidx(i):end,:),1,size(raln1{i},2)), 'color', c3, 'LineWidth',0.5);
    scatter3(raln1{i}(1:sidx(i),:), raln2{i}(1:sidx(i),:), repmat(t{i}(1:sidx(i),:),1,size(raln1{i},2)), 30, c1);hold on;
    scatter3(raln1{i}(sidx(i):eidx(i),:), raln2{i}(sidx(i):eidx(i),:), repmat(t{i}(sidx(i):eidx(i),:),1,size(raln1{i},2)), 30, c2);
    scatter3(raln1{i}(eidx(i):end,:), raln2{i}(eidx(i):end,:), repmat(t{i}(eidx(i):end,:),1,size(raln1{i},2)), 30, c3);
    
    % Aligned hemodynamic coupling
    figure(fig6);
    scatter(interp(raln1{i}(:),60), interp(raln2{i}(:),60), 1, interp(bp{i}(:),60), 'filled');hold on;
    scatter(raln1{i}(:), raln2{i}(:), 30, bp{i}(:));hold on;
    set(fig4, 'Colormap', colormap(viridis(300)));
    clim([-1 1]);
    axis equal;
    axis off;
end

%% Dynamics of population state
baseC = [0 0 0];
nBin = AnalysisStatTrans.nBin;

figure;
polarhistogram(unalgnPre, nBin, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k');
hold on;
h = polarhistogram(algnPre, nBin, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'FaceColor', c1, 'FaceAlpha', 0.8, 'EdgeColor', 'k');
m = max(h.Values);
r = 1.05*m;
p = (0:0.02:1.6); p = p(p>r); p = p(1);
polarplot(AnalysisDynmDist.tBmean(1)*[0 1], [0 0.9*p], 'k', 'LineWidth', 3);
polarplot(AnalysisDynmDist.tTmean(1)*[0 1], [0 0.9*p], 'g', 'LineWidth', 3);
h.Parent.RAxis.TickValues = 0;
h.Parent.ThetaAxis.TickValues = [0 90 180 270 360];

figure;
polarhistogram(unalgnStim, nBin, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k');
hold on;
h = polarhistogram(algnStim, nBin, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'FaceColor', c2, 'FaceAlpha', 0.8, 'EdgeColor', 'k');
m = max(h.Values);
r = 1.05*m;
p = (0:0.02:1.6); p = p(p>r); p = p(1);
polarplot(AnalysisDynmDist.tBmean(2)*[0 1], [0 0.9*p], 'k', 'LineWidth', 3);
polarplot(AnalysisDynmDist.tTmean(2)*[0 1], [0 0.9*p], 'g', 'LineWidth', 3);
h.Parent.RAxis.TickValues = 0;
h.Parent.ThetaAxis.TickValues = [0 90 180 270 360];

figure;
polarhistogram(unalgnPost, nBin, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'FaceColor', 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'k');
hold on;
h = polarhistogram(algnPost, nBin, 'Normalization', 'probability', 'DisplayStyle', 'bar', 'FaceColor', c3, 'FaceAlpha', 0.8, 'EdgeColor', 'k');
m = max(h.Values);
r = 1.05*m;
p = (0:0.02:1.6); p = p(p>r); p = p(1);
polarplot(AnalysisDynmDist.tBmean(3)*[0 1], [0 0.9*p], 'k', 'LineWidth', 3);
polarplot(AnalysisDynmDist.tTmean(3)*[0 1], [0 0.9*p], 'g', 'LineWidth', 3);
h.Parent.RAxis.TickValues = 0;
h.Parent.ThetaAxis.TickValues = [0 90 180 270 360];

%% State transition
figure;
ax = axes();
hold(ax);
boxchart(ones(size(AnalysisStatTrans.tranStd(1,:))), AnalysisStatTrans.tranStd(1,:),...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(AnalysisStatTrans.tranStd(2,:))), AnalysisStatTrans.tranStd(2,:),...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,nBin), AnalysisStatTrans.tranStd, 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(nBin, 1)+0.2*(rand(nBin,1)-0.5), AnalysisStatTrans.tranStd(1,:), 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(nBin, 1)+0.2*(rand(nBin,1)-0.5), AnalysisStatTrans.tranStd(2,:), 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');

disp('Standard deviation of state transition probability distribution for each population state phase (unaligned vs. aligned) (Wilcoxon signed-rank test)');
disp(['p = ', num2str(signrank(AnalysisStatTrans.tranStd(1,:), AnalysisStatTrans.tranStd(2,:)))]);

%% Hemodynamic coupling
figure;
ax = axes();
hold(ax);
boxchart(ones(size(AnalysisFuncDist.bpStd(1,:))), AnalysisFuncDist.bpStd(1,:),...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(AnalysisFuncDist.bpStd(2,:))), AnalysisFuncDist.bpStd(2,:),...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,nBin), AnalysisFuncDist.bpStd, 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(1,nBin)+0.2*(rand(1,nBin)-0.5), AnalysisFuncDist.bpStd(1,:), 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(1,nBin)+0.2*(rand(1,nBin)-0.5), AnalysisFuncDist.bpStd(2,:), 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');

disp('Standard deviation of blood pressure distribution for each population state phase (unaligned vs. aligned) (Wilcoxon signed-rank test)');
disp(['p = ', num2str(signrank(AnalysisFuncDist.bpStd(1,:), AnalysisFuncDist.bpStd(2,:)))]);

%% Representative decoding the hemodynamics
figure;
plot(t{1}, mean(bp{1},2), 'LineWidth', 3, 'Color', 'k');
hold on;
plot(t{1}, mean(YIndv{1},2), 'LineWidth', 3, 'Color', 'g');
plot(t{1}, mean(YComm{1},2), 'LineWidth', 7, 'Color', 'y');
plot(t{1}, mean(YComm{1},2), 'LineWidth', 3, 'Color', 'r');
box off;
xlim([-20 80]);
ylim([-1 1]);

%% Goodness-of-fitting in linear regressions
figure;
ax = axes();
hold(ax);
boxchart(ones(size(accuIndv)), accuIndv,...
    'BoxFaceColor', 'g', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(accuComm)), accuComm,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [accuIndv; accuComm], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1)+0.2*(rand(10,1)-0.5), accuIndv, 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(10, 1)+0.2*(rand(10,1)-0.5), accuComm, 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
ylim([0 1]);

disp('R^2 of rat-specific and across-rat linear regressions (Wilcoxon signed-rank test)');
disp(['p = ', num2str(signrank(accuIndv, accuComm))]);