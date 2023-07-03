clc;
clear;
close all;

addpath('..\FUNCTIONS\');

load('..\PARAMETER_DATA\params.mat');

color = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];...
    [0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840];...
    [0.5 0 0];[0 0.5 0];[0 0 0.5]];

dt = 0.5;
k = 1;
frALL = cell(numel(rat),1);
bpALL = cell(numel(rat),1);
for i = 1:numel(rat)
    ratNumb = rat{i};
    
    pre = PRE(i);
    post = POST(i);
    tr = TRIAL(i);
    trIdx = TrialIDX{i};
    
    
    %% BP data processing
    bpData = load(['..\PROCESSED\BP_DATA\bp_data_', ratNumb, '_v7.mat']);

    % BP data synchronization
    bp = bpData.BP(:,2);
    

    % BP data filtering (get MBP)
%     [LF_num, LF_den] = butter(2, 0.5/(100*0.5), 'low');
%     bp = filtfilt(LF_num, LF_den, bp);
    %bp = gaussianFilt(bp, 0.01, 0.1, round(3/0.01)); % filter sigma = dt, filter length = round(3/dt)
    
    % BP data re-alignment
    len = (post-pre)*100;
    bp = reshape(bp, len, numel(bp)/len);
    bp = bp(:,trIdx);
    tBP = (0:size(bp,1)-1)*0.01 + pre;

    [~,idx1] = min(abs(tBP-(-20)));
    idx2 = idx1 + 100*100;
    % MBP data normalization
    for j = 1:size(bp,2)
        bp(:,j) = minmaxNorm(bp(:,j), [-1 1]);
    end    
    bpALL{i} = minmaxNorm(mean(bp(idx1:idx2,:),2), [-1 1]);
    %% Spike data processing
    spkData = load(['..\PROCESSED\SPIKE_DATA\spike_data_', ratNumb, '_v7.mat']);
    spk = spkData.spike;
    T = post-pre;
    edges = (0:dt:T*trALL(i));
    t = 0.5*(edges(2:end)+edges(1:end-1));
    t = t(1:numel(t)/trALL(i)) + pre;

    [~, idx1] = min(abs(t-(-20)));
    idx2 = idx1 + 100/dt;

    % Spike data screening
    fr = zeros(numel(spk), idx2-idx1+1);
    gaussianParams = struct('sigma', 0.5, 'length', round(3/dt));
    for j = 1:numel(spk)
        frTmp = reshape(calFiringRate(spk{j}, edges, gaussianParams), numel(t), trALL(i));
        fr(j,:) = minmaxNorm(mean(frTmp(idx1:idx2, trIdx),2), [-1 1]);
    end
    frALL{i} = fr;
end

%%
figure; hold on;
bp = [];
for j = 1:numel(bpALL)
    bp = [bp, bpALL{j}];
end
tBP = (0:size(bp,1)-1)*0.01 - 20;
plot(tBP, bp, 'LineWidth',1, 'Color',0.5*[1 1 1]);
plot(tBP, mean(bp,2), 'LineWidth',3, 'Color', 'r');

figure; hold on;
for j = 1:size(bp,2)
    ang(:,j) = asin(bp(:,j));
    [~, midx] = min(ang(:,j));
    ang(1:midx,j) = -pi - ang(1:midx,j);
end
plot(tBP, ang, 'LineWidth',1, 'Color',0.5*[1 1 1]);
plot(tBP, mean(ang,2), 'LineWidth',3, 'Color', 'r');


%%
fr = [];
for j = 1:numel(frALL)
    fr = [fr; frALL{j}];
end
tFR = (0:size(fr,2)-1)*dt - 20;

[~, mxIdx] = max(fr, [], 2);
mxIdx = mod(mxIdx - 20/dt,numel(tFR));
[~, srtIdx] = sort(mxIdx);
figure; hold on;
for j = 1:numel(srtIdx)
    imagesc(tFR, j*ones(size(tFR)), fr(srtIdx(j),:));
end

figure; hold on; % for the PSTH heatmap for all trials of excitatory neurons
sidx = 20/dt;
eidx = 80/dt;
for j = 1:size(fr,1)
    score_ex(j) = fr(j,:)*[-20*ones(sidx,1);(eidx-sidx:-1:1).^5';-2*ones(numel(tFR) - eidx,1)];
end
[~, srt_idx_ex] = sort(score_ex);

for j = 1:size(fr,1)
    imagesc(tFR, j*ones(size(tFR)), fr(srt_idx_ex(j),:));
end

ylim([0.5 size(fr,1)+0.5]);
xlim([-20 80]);
colormap viridis(100);

%%
ang = ang(1:50:end,:);

edges = (0:0.001:1);
for j = 1:numel(rat)
    corTmp = [];
    for k= 1:size(ang,2)
        tmp = corrcoef(ang(:,j), frALL{j}(k,:));
        corTmp = [corTmp, tmp(1,2)];
    end
    cor{j} = corTmp.^2;
    cnt(:,j) = histcounts(corTmp, edges);
    cnt(:,j) = cumsum(cnt(:,j))/sum(cnt(:,j));
    plot(log10(edges(2:end)), cnt(:,j), 'k');
end


%%
figure;
cor = zeros(numel(rat),1);
corR = zeros(numel(rat),1);
corALL = [];
corRALL = [];
edges = (0:0.001:1);
for i = 1:numel(frALL)
    for j = 1:size(frALL{i},1)
        corTmp = [];
        corShf = [];
        for k = 1:size(frALL{i},1)
            tmp = corrcoef(frALL{i}(j,:), frALL{i}(k,:));
            corTmp = [corTmp, abs(tmp(1,2))];
            
            L = size(frALL{i},2);
            tmp = corrcoef(frALL{i}(j,:), frALL{i}(k,randperm(L)));
            corShf = [corShf, abs(tmp(1,2))];
        end
    end
    cnt = histcounts(corTmp, edges);
    corALL = [corALL;cumsum(cnt)/sum(cnt)];
    cor(i) = mean(corTmp);

    cR = histcounts(corShf, edges);
    corRALL = [corRALL;cumsum(cR)/sum(cR)];
    corR(i) = mean(corShf);

end
hold on;

plot(log10(edges(2:end)), corALL', 'Color', [1 0 0 0.2], 'LineWidth',1);
plot(log10(edges(2:end)), mean(corALL,1), 'r', 'LineWidth',3);
plot(log10(edges(2:end)), corRALL', 'Color', [0 0 0 0.2], 'LineWidth',1);
plot(log10(edges(2:end)), mean(corRALL,1), 'k', 'LineWidth',3);

figure;
ax = axes();
hold(ax);
boxchart(ones(size(corR)), corR,...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(cor)), cor,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [corR' ; cor'], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1)+0.2*(rand(10,1)-0.5), corR, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(10, 1)+0.2*(rand(10,1)-0.5), cor, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
ylim([0 1]);
xlim([0.5 2.5]);

%%
var = zeros(numel(rat),1);
varR = zeros(numel(rat),1);
varALL = [];
varRALL = [];
edges = (0:0.001:1);
for i = 1:numel(frALL)
        [lambda,psi] = factoran(frALL{i}',6);
        varTmp = lambda*lambda'./(lambda*lambda' + diag(psi));
        varTmp = diag(varTmp);
        cnt = histcounts(varTmp, edges);
        varALL = [varALL;cumsum(cnt)/sum(cnt)];
        var(i) = mean(varTmp);
            
        L = size(frALL{i},2);
        frShf = zeros(size(frALL{i}));
        for j = 1:size(frShf,1)
            frShf(j,:) = frALL{i}(j,randperm(L));
        end
        [lambda,psi] = factoran(frShf',6);
        varShf = lambda*lambda'./(lambda*lambda' + diag(psi));
        varShf = diag(varShf);
        cR = histcounts(varShf, edges);
        varRALL = [varRALL;cumsum(cR)/sum(cR)];
        varR(i) = mean(varShf);
        
        clc;
end

figure;
hold on;
plot(edges(2:end), varALL, 'Color', [1 0 0 0.2], 'LineWidth',1);
plot(edges(2:end), mean(varALL,1), 'r', 'LineWidth',3);
plot(edges(2:end), varRALL, 'Color', [0 0 0 0.2], 'LineWidth',1);
plot(edges(2:end), mean(varRALL,1), 'k', 'LineWidth',3);

figure;
ax = axes();
hold(ax);
boxchart(ones(size(varR)), varR,...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(size(var)), var,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
 
plot(repmat([1;2],1,10), [varR' ; var'], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(10, 1)+0.2*(rand(10,1)-0.5), varR, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2*ones(10, 1)+0.2*(rand(10,1)-0.5), var, 100,...
    'MarkerFaceColor','w', 'MarkerEdgeColor','k');
% ylim([0 0.6]);
xlim([0.5 2.5]);
