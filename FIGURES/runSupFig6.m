clear;
close all;
clc;

addpath('..\FUNCTIONS\');

params = load('..\PARAMETER_DATA\params.mat');

dat = load('..\POSTPROCESSED_MAT_CODE\DATA\overallPopulationAnalysis_matlab_save.mat');
allpop = dat.overallPopulationAnalysis;

%%
range = (-20:0.5:100);
tL = numel(range);
tr = 3;
cellN = size(allpop.firing_rate,2);

figure;
X = rand(cellN,1);
Y = rand(cellN,1);
scatter(X,Y, 80, 'r', 'filled', 'MarkerEdgeColor','k');
axis off;
iso = allpop.fr_isomap;
for i = 1:size(allpop.bp_data,2)
    bp(:,i) = minmaxNorm(allpop.bp_data(:,i),[-1 1]);
end

iso1 = dataFold(iso(:,1), tL, tr);
iso2 = dataFold(iso(:,2), tL, tr);

figure;
scatter(iso1(:),iso2(:), 80, mean(bp,2), 'filled');
colormap viridis(100);
axis equal;
axis off;

%
h1 = allpop.h1;
h2 = allpop.h2;
h1_len = h1(:,2) - h1(:,1);
h2_len = h2(:,2) - h2(:,1);

f = figure; f.Position(3) = 280;
plot(h1', repmat((1:size(h1,1)),2,1), 'LineWidth', 3, 'Color',[1 0 0 0.7]);
xlim([0 20]);box off;
f = figure; f.Position(3) = 280;
plot(h2', repmat((1:size(h2,1)),2,1), 'LineWidth', 3, 'Color',[0 0 0 0.7]);
xlim([10 20]);box off;

%
expVar = var(iso,[],1);
totVar = sum(expVar);
f = figure; f.Position(3) = 280;
plot(expVar/totVar*100, 'Marker', 'o', 'Color','k', 'LineWidth',3);
xlim([0.5 6.5]);
box off;

%
f = figure; f.Position(3) = 280; hold on;
    x = (-5:0.01:5);
    for j = 1:30
        plot(x, 1*x-4+0.2*j,'k');
    end
figure(f);
fitting = fitlm(allpop.th(7:end-7), allpop.corr_integral(7:end-7));
slope = table2array(fitting.Coefficients(2,1));
bias =  table2array(fitting.Coefficients(1,1));
plot(x, slope*x+bias, 'g', 'LineWidth',1);
scatter(allpop.th(7:end-7), allpop.corr_integral(7:end-7), 50, 'r', 'filled');
ylim([-0.5 1.5])
xlim([0.2 2])

%%
figure;
scatter(X(29:end),Y(29:end), 80, 'r', 'filled', 'MarkerEdgeColor','k');
axis off;

subsmp = dat.subsamplingAnalysis;
refpop = subsmp.reference;
tstpops = subsmp.test;

figure;
refiso = refpop.fr_isomap;
scatter(refiso(:,1), refiso(:,2), 80, mean(bp,2), 'filled');
colormap viridis(100);
axis equal
xlim([-3 3])
ylim([-3 3])

crr = zeros(numel(tstpops),1);
for i = 1:numel(tstpops)
    pop = tstpops{i};
    popiso = pop.fr_isomap;
    isoCrr1 = corrcoef(refiso(:,1), popiso(:,1));
    isoCrr2 = corrcoef(refiso(:,2), popiso(:,2));
    crr(i) = mean([isoCrr1(1,2),isoCrr2(1,2)].^2);
end
figure;
plot(crr);
[~, midx] = max(crr);

mxpop = tstpops{midx};

figure; hold on;
scatter(X(29:end),Y(29:end), 80, 'k', 'filled', 'MarkerEdgeColor','k');
scatter(X(mxpop.subpopIdx),Y(mxpop.subpopIdx), 80, 'r', 'filled', 'MarkerEdgeColor','k');
axis off;

figure;
mxiso = mxpop.fr_isomap;
scatter(mxiso(:,1), mxiso(:,2), 80, mean(bp,2), 'filled');
colormap viridis(100);
axis equal
xlim([-3 3])
ylim([-3 3])