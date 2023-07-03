clear;
close all;
clc;

addpath('..\FUNCTIONS\');

rat = {...
    '19', '21', '22',...
    'B2', 'B3', 'B4', 'B9_a1',...
    'B10_a1', 'B14_a1', 'B15_a1'};
DT = [0.6 0.5 0.75,...
    0.6 0.6 0.6 0.6,...
    0.6 0.6 0.6];
ntsDELAY = -1*[3 3 3,...
    3 3 3 3,...
    3 5 3];
bpDELAY = [2 3 2,...
    9 9 9 2,...
    0 0 0] - ntsDELAY;
PRE = [-20 -20 -20,...
    -20 -20 -20 -20,...
    -20 -20 -40];
POST = [90 90 90,...
    120 120 120 120,...
    120 100 160];
TRIAL = [4 3 3,...
    4 3 3 4,...
    4 4 3];
TrialIDX = {...
    [1 2 3 4], [1 2 3], [1 2 3],...
    [1 2 3 4],[10 11 12], [2 4], [3 4 5 6],...
    [3 4 5 6],[3 4 5 6],[2 3 4]};

color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];...
    [0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840];...
    [0.5 0 0];[0 0.5 0];[0 0 0.5]];

c1 = [42 183 202]/255;
c2 = [0.9290 0.6940 0.1250];
c3 = [254 74 73]/255;%[0.7 0.7 0.7];

figure;
i = 1;
    ratNumb = rat{i};
    data = load(['..\POSTPROCESSED_MAT_CODE\DATA\', ratNumb, '_matlab_save.mat']);

    dt = DT(i);
    tr = TRIAL(i);
    bDelay = bpDELAY(i);
    nDelay = ntsDELAY(i);
    pre = PRE(i);

    [t, bp, ~, fr, iso1, iso2] = dataPrepare(data, dt, tr, nDelay, bDelay, pre);
    [~, sidx] = min(abs(t+20));
    [~, eidx] = min(abs(t-80));

    hold on;
    tmpBP = bp(sidx:eidx,2:4);
    tmp1 = iso1(sidx:eidx,2:4);
    tmp2 = iso2(sidx:eidx,2:4);
    tmpC = t(sidx:eidx,2:4);
    scatter(tmp1(:), tmp2(:), 80, tmpC(:), 'filled');
    axis equal;
    axis off;
    colormap hot; 
    %colormap(viridis(100));

f = figure; hold on; f.Position(3:4) = [450 340];
plot(tmpC(:,1), mean(tmpBP,2), 'Color', 'k', 'LineWidth', 3);
% scatter(interp(tmpC(:,1),10), interp(mean(tmpBP,2),10), 50, interp(tmpC(:,1),10), 'filled');
colormap hot; %(viridis(100));
xlim([-20 80]);

% figure;
% scatter(t(1:L), bp_data(1:L), 50, t(1:L), 'filled');
% colormap cool;
% box off;
% xlim([0 60]);


%%
plot_prcnt = [0, 90, 20];
j = 1;
    ratNumb = rat{j};
    data = load(['..\POSTPROCESSED_MAT_CODE\DATA\', ratNumb, '_matlab_save.mat']);

    h0 = data.h0;
    h1 = data.h1;
    h2 = data.h2;
    
    h0_len = h0(:,2) - h0(:,1);
    h0_s = h0(h0_len > prctile(h0_len, plot_prcnt(1)),:);
    h1_len = h1(:,2) - h1(:,1);
    h1_s = h1(h1_len > prctile(h1_len, plot_prcnt(2)),:);
    h2_len = h2(:,2) - h2(:,1);
    h2_s = h2(h2_len > prctile(h2_len, plot_prcnt(3)),:);
    
%     subplot(3,1,1);hold on;
%     for i = 1:size(h0_s,1)
%         line([h0_s(i,1) h0_s(i,2)], [size(h0_s,1)-i size(h0_s,1)-i], 'LineWidth', 4, 'Color', [115 160 250]/255);
%     end
%     xlim([-0.01 35]); ylim([-2 25]);
    
figure; hold on;
for i = 1:size(h1_s,1)
    line([h1_s(i,1) h1_s(i,2)], [size(h1_s,1)-i size(h1_s,1)-i], 'LineWidth', 6, 'Color', [115 160 250]/255);
end
xlim([-0.01 60]); ylim([-2 18]);
    
figure; hold on;
for i = 1:size(h2_s,1)
    line([h2_s(i,1) h2_s(i,2)], [size(h2_s,1)-i size(h2_s,1)-i], 'LineWidth', 6, 'Color', [115 160 250]/255);
end
xlim([-0.01 60]); ylim([-2 12]);

