clear;
close all;
clc;

addpath('..\FUNCTIONS\');
param = load('..\PARAMETER_DATA\params.mat');

% color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];...
%     [0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840];...
%     [0.5 0 0];[0 0.5 0];[0 0 0.5]];
color = zeros(10, 3);
for i = 1:10
    color(i,:) = rand(1,3);
end


c1 = [42 183 202]/255;
c2 = [0.9290 0.6940 0.1250];
c3 = [254 74 73]/255;%[0.7 0.7 0.7];


i = 1;
data = load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\expDataAnalysis.mat');
ratNumb = param.rat{i};

%% Representative single-neuron BP prediction
algn1 = data.algn1{i};
algn2 = data.algn2{i};

fr = data.fr{i};
cellN = size(fr,1);
bp = data.bp{i};
t = data.t{i};
for j = 1:cellN
    frR(j,:) = minmaxNorm(fr(j,:), [-1 1]);
end

figure; hold on;
[~, mxIdx] = max(frR, [], 2);
[~, srtIdx] = sort(mxIdx);
for j = 1:cellN
    imagesc(t, ones(numel(t),1)*j, frR(srtIdx(j),:));
end
ylim([0.5 cellN+0.5]);
xlim([t(1), t(end)]);

% representative BP(measured) vs. prediction by single-neurons
figure; hold on;
YSingle = data.YSingle{i};
% scatter(Y1, repmat(bp(:)',size(Y1,1),1), 30, 0.6*[1 1 1], 'filled', 'MarkerEdgeColor', 'k');
plot(t, YSingle, 'LineWidth', 0.5, 'Color', c1);
% plot(t, mean(bp,2), 'Color', 'r', 'LineWidth', 5);
plot(t, mean(bp,2), 'Color', 'k', 'LineWidth', 3);
ylim([-1 1]);
xlim([-20 80]);

% representative prediction accuracy by single-neuron firing rates
figure;
accuSing = data.accuSing{i};
bar((1:cellN),accuSing(srtIdx)); box off;
yline(mean(accuSing));
disp('Representative prediction accuracy with single neurons:');
disp(mean(accuSing));
xlim([0.5 cellN+0.5]);
ylim([0 1]);

%% Representative latent space

iso1 = data.iso1{i};
iso2 = data.iso2{i};

% neural trajectory
figure;
scatter(iso1(:), iso2(:), 100, bp(:), 'filled');
colormap viridis(100);
clim([-1 1]);
xlim([min(iso1(:)) max(iso1(:))]);
ylim([min(iso2(:)) max(iso2(:))]);
axis off;
axis equal;

% 2D latent space BP mapping
[b, Y1, coef, ~] = linRegress([ones(numel(bp(:)),1), iso1(:), iso2(:)], bp(:));
Y1 = reshape(Y1, size(bp,1), size(bp,2));
disp('Latent space decoding accuracy: ');
disp(coef^2);
figure; hold on;
[X, Y] = meshgrid((min(iso1(:)):0.1:max(iso1(:))), (min(iso2(:)):0.1:max(iso2(:))));
scatter(X(:), Y(:), 10, b(1) + b(2)*X(:) + b(3)*Y(:),'filled');
view(0, 90);
clim([-1 1]);
colormap viridis(100);
scatter(iso1(:), iso2(:), 80,'k', 'filled');
plot([mean(X(1,:))+4, mean(X(1,:))+4+100*b(2)], [mean(Y(:,1))+2, mean(Y(:,1))+2+100*b(3)], 'Color', 'w', 'LineWidth', 4);
xlim([min(iso1(:)) max(iso1(:))]);
ylim([min(iso2(:)) max(iso2(:))]);
axis equal;

% projection on decoding space
figure;
scatter([iso1(:), iso2(:)]*[b(2);b(3)], 0.5*(rand(numel(iso1(:)),1)-0.5), 80, bp(:), 'filled')
ylim([-3 3])
colormap viridis(100);
clim([-1 1]);
axis off;

% BP prediction by decoding of neural trajectory
figure;
plot(t, mean(bp,2), 'Color', 'k', 'LineWidth',3);
hold on;
plot(t, mean(Y1,2), 'Color', 'y', 'LineWidth',7);
plot(t, mean(Y1,2), 'Color', 'r', 'LineWidth',3);
xlim([-20 80]);
box off;

% XY graph for BP(measured) vs. BP prediction by decoding
figure; hold on;
scatter([iso1(:), iso2(:)]*[b(2);b(3)]+b(1), bp(:), 50, 'k', 'filled');
plot((-2:0.01:2), (-2:0.01:2), 'LineWidth',3, 'Color','r');
xlim([-1 1]);
ylim([-1 1]);

% prediction error distribution (single-neuron vs. decoding)
fig = figure; hold on;
edges = (0.0001:0.0001:2);
sqrErrSing = data.sqrErrSing{i}(:);
sqrErrIndv = data.sqrErrIndv{i};
drawCumErrFrac(fig, log10(edges), log10(sqrErrSing), [0 0 0]);
xline(log10(mean(sqrErrSing)), 'Color','k', 'LineWidth',3, 'LineStyle',':');
drawCumErrFrac(fig, log10(edges), log10(sqrErrIndv), [1 0 0]);
xline(log10(mean(sqrErrIndv)), 'Color','r', 'LineWidth',3, 'LineStyle',':');
ylim([-0.01 1.1]);

%%
fig1 = figure; hold on;
axis equal;
xlim([-1 1]);
ylim([-1 1]);
xline(0); yline(0);
axis off;

fig2 = figure; hold on;
axis equal; xline(0); yline(0);
xlim([-1 1]);
ylim([-1 1]);
fig3 = figure; hold on;
axis equal; xline(0); yline(0);
xlim([-1 1]);
ylim([-1 1]);

%% Statistics
for i = 1:numel(param.rat)
    %% parameter setting
    ratNumb = param.rat{i};
    dt = param.DT(i);
    tr = param.TRIAL(i);
    bDelay = param.bpDELAY(i);
    nDelay = param.ntsDELAY(i);
    pre = param.PRE(i);

    %% data setting
    t = data.t{i};
    iso1 = data.iso1{i};
    iso2 = data.iso2{i};
    bp = data.bp{i};
    fr = data.fr{i};

    %% latent sapce decoding
    % plot xy plot of latent space decoding
    figure(fig2);
    scatter(data.YIndv{i}, mean(bp,2), 30, [1 0.3 0.3], 'filled', 'MarkerEdgeColor', 'k');
    
    %% single neuron decoding
    %method 1: every single neuron vs. BP
    %method 2: closest estimation among single neurons at each time step vs. BP
    %method 3: average estimation across single neuorn vs. BP
    figure(fig3);
    scatter(mean(data.YSingle{i},1), mean(bp,2), 30, 0.6*[1 1 1], 'filled', 'MarkerEdgeColor', 'k');
    
    %% plot decoding space
    figure(fig1);
    quiver(0, 0, data.decVec(i,1), data.decVec(i,2), 0, 'Color', color(i,:), 'LineWidth', 3, 'MaxHeadSize',0.5);
end

%%
figure; hold on;
decVec = data.decVec;
ang = atan2(decVec(:,2), decVec(:,1));
ang(ang<0) = ang(ang<0) + 2*pi;
ang = ang/pi;
bar(1, mean(ang), 'BarWidth', 0.5, 'FaceColor', 'g', 'FaceAlpha',0.5);
errorbar(1, mean(ang), 0, std(ang), 'LineWidth', 2, 'Color', 'k');
scatter(ones(size(ang))+0.2*(rand(size(ang))-0.5), ang, 50, 'w', 'filled', 'MarkerEdgeColor','k');
ylim([0 2]);
xlim([0.5 1.5]);

%%
accuSing = zeros(1,numel(param.rat));
for i = 1:numel(param.rat)
    accuSing(i) = mean(data.accuSing{i});
end
accuIndv = data.accuIndv;

figure;
ax = axes();
hold(ax);
boxchart(ones(size(accuSing)), accuSing,...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(accuIndv)), accuIndv,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [accuSing; accuIndv], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1)+0.2*(rand(10,1)-0.5), accuSing, 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(10, 1)+0.2*(rand(10,1)-0.5), accuIndv, 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
ylim([0 1]);
disp('Accuracy: sinlge-neuron vs. decoding (Wilcoxon signed-rank test)');
disp(['p = ', num2str(ranksum(accuSing,accuIndv))]);

%%
fig = figure; hold on;
edges = (0.0001:0.0001:4);
seSing = [];
seIndv = [];
for i = 1:numel(param.rat)
    seSing = [seSing;mean(data.sqrErrSing{i},1)];
    seIndv = [seIndv;data.sqrErrIndv{i}'];
end
drawCumErrFrac(fig, log10(edges), log10(seSing'), [0 0 0]);
drawCumErrFrac(fig, log10(edges), log10(seIndv'), [1 0 0]);
mseSing = mean(seSing,2);
mseIndv = mean(seIndv,2);
xline(log10(mean(mseSing)), 'Color','k', 'LineWidth',3, 'LineStyle',':');
xline(log10(mean(mseIndv)), 'Color','r', 'LineWidth',3, 'LineStyle',':');

%%
figure;
ax = axes();
hold(ax);
mseSingLog = log10(mseSing);
mseIndvLog = log10(mseIndv);
boxchart(ones(size(mseSingLog)), mseSingLog,...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(mseIndvLog)), mseIndvLog,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [mseSingLog'; mseIndvLog'], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1)+0.2*(rand(10,1)-0.5), mseSingLog, 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(10, 1)+0.2*(rand(10,1)-0.5), mseIndvLog, 50,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');

disp('MESE: sinlge-neuron vs. decoding (Wilcoxon signed-rank test)');
disp(['p = ', num2str(ranksum(mseSingLog,mseIndvLog))]);