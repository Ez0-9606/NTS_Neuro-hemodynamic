clear;
close all;
clc;

params = load('..\PARAMETER_DATA\params.mat');
ex = [];
in = [];
cmp = [];
c = ...
    [0.7609    0.2132    0.9060
    0.2166    0.7237    0.3695
    0.0230    0.2334    0.8859];

for i = 1:numel(params.rat)
    ratNumb = params.rat{i};
    dt = params.DT(i);
    pre = params.PRE(i);
    post = params.POST(i);
    tr = params.TRIAL(i);
    trIdx = params.TrialIDX{i};

    %%
    data = load(['..\PROCESSED\SPIKE_DATA\spike_data_', ratNumb, '_v7.mat']);
    actNumb(i) = numel(data.spike);
    nonactNumb(i) = numel(data.spikeOut);

    isiOut = [];
    for j = 1:numel(data.spikeOut)
        tmp = data.spikeOut{j};
        tmp = mod(tmp, post-pre);
        isiTmp = tmp(2:end)- tmp(1:end-1);
        isiTmp(isiTmp<0) = [];
        isiOut = [isiOut; isiTmp];
    end
    mISIout(i) = mean(isiOut);
    isiALLout{i} = isiOut;

    isi = [];
    for j = 1:numel(data.spike)
        tmp = data.spike{j};
        tmp = mod(tmp, post-pre);
        isiTmp = tmp(2:end)- tmp(1:end-1);
        isiTmp(isiTmp<0) = [];
        isi = [isi; isiTmp];
    end
    mISI(i) = mean(isi);
    isiALL{i} = isi;

    %%
    data = load(['..\POSTPROCESSED_MAT_CODE\DATA\', ratNumb, '_matlab_save.mat']);
    frTmp = data.firing_rate;
    fr = zeros(size(frTmp,2), floor((post-pre)/dt));
    for j = 1:size(frTmp,2)
        tmp = dataFold(frTmp(:,j), floor((post-pre)/dt), tr);
        fr(j,:) = mean(tmp,2);
    end
    fr = fr - mean(fr(:,1:round(dt*-pre)),2);
    fr = fr/max(max(abs(fr)));
    t = (0:floor((post-pre)/dt)-1)*dt+pre;

    range = (-20:0.5:80);
    idx = getNormIdx(t, range);
    t = t(idx);
    fr = fr(:,idx);


    stimTime = 60;
    if i == 4
        stimTime = 15;
    elseif i == 5 || i == 6
        stimTime = 30;
    end
    
    [~, maxidx] = max(fr,[],2);
    [~, minidx] = min(fr,[],2);

    exIdx = (maxidx >= (-pre)/dt) & (maxidx < (stimTime-pre)/dt);
    inIdx = (minidx >= (-pre)/dt) & (minidx < (stimTime-pre)/dt);
    exIdx(exIdx & inIdx) = false;
    inIdx(exIdx & inIdx) = false;
    cmpIdx = ~exIdx & ~inIdx;

    ex = [ex; mean(fr(exIdx,:),1)];
    wex(i) = sum(exIdx);
    in = [in; mean(fr(inIdx,:),1)];
    win(i) = sum(inIdx);
    cmp = [cmp; mean(fr(cmpIdx,:),1)];
    wcmp(i) = sum(cmpIdx);


    plot(t, fr(exIdx,:)', 'Color', c(1,:)); hold on;
    plot(t, fr(inIdx,:)', 'Color', c(2,:));
    if sum(cmpIdx) ~= 0
        plot(t, fr(cmpIdx,:)', 'Color', c(3,:));
    end
    xlim([-20 80]); box off;
end

%%
figure; hold on;
b = bar([nonactNumb', actNumb'], 'stacked');
b(1).FaceColor = 'k';
b(2).FaceColor = 'r';

%%
figure; hold on;
b = bar([wex', win', wcmp'], 'stacked');
b(1).FaceColor = c(1,:);
b(2).FaceColor = c(2,:);
b(3).FaceColor = c(3,:);

%%
figure;
repreTime = load('..\PROCESSED\SPIKE_DATA\spike_data_19_v7.mat');
repreWave = load('..\PROCESSED\SPIKE_DATA\wave_data_19_v7.mat');
subplot(2,1,1); hold on;
scatter(repreTime.spike{1}, zeros(numel(repreTime.spike{1}),1), 600, 'Marker', '|', 'MarkerEdgeColor', 'r');
for i = 1:4
    line([20 80]+110*(i-1), [1 1], 'LineWidth', 0.5);
end
ylim([-1 1]);
subplot(2,1,2);
scatter(repreTime.spikeOut{2}, zeros(numel(repreTime.spikeOut{2}),1), 600, 'Marker', '|', 'MarkerEdgeColor','k');
ylim([-1 1]);


figure;
subplot(2,1,1); hold on;
scatter(repreTime.spike{1}, zeros(numel(repreTime.spike{1}),1), 600, 'Marker', '|', 'MarkerEdgeColor', 'r', 'LineWidth',2);
ylim([-1 1]);xlim([129.75 130.25]);
for i = 1:5
    xline(130+0.05*(i-1), 'LineWidth', 0.5, 'Color', 'b', 'LineStyle', '--');
end
subplot(2,1,2); hold on;
scatter(repreTime.spikeOut{2}, zeros(numel(repreTime.spikeOut{2}),1), 600, 'Marker', '|', 'MarkerEdgeColor','k', 'LineWidth',2);
ylim([-1 1]);xlim([129.75 130.25]);
for i = 1:5
    xline(130+0.05*(i-1), 'LineWidth', 0.5, 'Color', 'b', 'LineStyle', '--');
end


figure;
subplot(2,1,1);hold on;
plot(repreWave.wave{1}', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
plot(mean(repreWave.wave{1},1), 'r', 'LineWidth', 3);
ylim([-200 100]);xlim([1 33]);
subplot(2,1,2);hold on;
plot(repreWave.waveOut{2}', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
plot(mean(repreWave.waveOut{2},1), 'k', 'LineWidth', 3);
ylim([-200 100]);xlim([1 33]);
%%
edges = (0.003:0.001:0.5);
fig1 = figure; hold on;
drawDist(fig1, edges, isiALL{1}, 'r');
fig2 = figure; hold on;
drawDist(fig2, edges, isiALLout{1}, 'k');

%%
fig = figure; hold on;
t = (-20:0.5:80);

x = [t, fliplr(t)];
wex = wex/sum(wex);
win = win/sum(win);
wcmp = wcmp/sum(wcmp);
mEx = weightedMean(ex,wex'); 
mIn = weightedMean(in,win');
mCmp = weightedMean(cmp,wcmp');
stdEx = 2*weightedStd(ex,wex');
stdIn = 2*weightedStd(in,win');
stdCmp = 2*weightedStd(cmp,wcmp');
yline(0);
fill(x, [mEx+stdEx, fliplr(mEx-stdEx)], c(1,:), 'FaceAlpha',0.2, 'LineStyle','none');
plot(t, mEx, 'Color', c(1,:), 'LineWidth', 3);
fill(x, [mIn+stdIn, fliplr(mIn-stdIn)], c(2,:), 'FaceAlpha',0.2, 'LineStyle','none');
plot(t, mIn, 'Color', c(2,:), 'LineWidth', 3);
fill(x, [mCmp+stdCmp, fliplr(mCmp-stdCmp)], c(3,:), 'FaceAlpha',0.2, 'LineStyle','none');
plot(t, mCmp, 'Color', c(3,:), 'LineWidth', 3);
