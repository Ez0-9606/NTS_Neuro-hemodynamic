clear;
close all;
clc;

addpath('..\FUNCTIONS\');
data = load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\expDataAnalysis.mat');
params = load('..\PARAMETER_DATA\params.mat');
N = numel(params.rat);
c = ...
    [0.7792    0.9340    0.1299
    0.0353    0.8510    0.9333
    247/255 110/255 202/255];

fr = cell(numel(params.rat),1);


corbp = cell(numel(params.rat),1);
corhr = cell(numel(params.rat),1);
corhrv = cell(numel(params.rat),1);
flag = cell(numel(params.rat),1);

for i = 1:N
    % parameter setting
    ratNumb = params.rat{i};
    dt = params.DT(i);
    pre = params.PRE(i);
    post = params.POST(i);
    tr = params.TRIAL(i);
    nDelay = params.ntsDELAY(i);
    bDelay = params.bpDELAY(i);
    trIdx = params.TrialIDX{i};
    
    %%
    delay(i) = bDelay*dt;

    %% BP data processing
    bp(:,i) = minmaxNorm(mean(data.bp{i},2),[-1 1]);
    hr(:,i) = minmaxNorm(mean(data.hr{i},2),[-1 1]);
    hrvLFpw(:,i) = data.hrv{i,1};
    hrvHFpw(:,i) = data.hrv{i,2};

    %%
    maccubp(i) = mean(data.accuSing{i});
    maccuhr(i) = mean(data.accuSingHR{i});
    maccuhrvLFpw(i) = mean(data.accuSingHRVlfpw{i});
    maccuhrvHFpw(i) = mean(data.accuSingHRVhfpw{i});
    stdaccubp(i) = std(data.accuSing{i});
    stdaccuhr(i) = std(data.accuSingHR{i});
    stdaccuhrvLFpw(i) = std(data.accuSingHRVlfpw{i});
    stdaccuhrvHFpw(i) = std(data.accuSingHRVhfpw{i});

end
minBP = data.minBP;
minHR = data.minHR;
minLF = data.minHRVlf;
minHF = data.minHRVhf;
maxBP = data.maxBP;
maxHR = data.maxHR;
maxLF = data.maxHRVlf;
maxHF = data.maxHRVhf;

%% BP delay plot
% representative plot
repreDat = load('..\POSTPROCESSED_MAT_CODE\DATA\19_matlab_save.mat');
[tR, bpR, ~, ~, ~, ~] = dataPrepare(repreDat, params.DT(1), params.TRIAL(1), params.ntsDELAY(1), 0, params.PRE(1));
figure; hold on;
plot(tR(:,1), mean(bpR,2), 'LineWidth',3);
xlim([-5 10]);
ylim([-20 0]);
yticks([-20 -10 0]);
yticklabels({'-20', '-10', '0'});
xline(delay(1), 'LineWidth', 3, 'LineStyle', '--', 'Color','k');
xline(0, 'LineWidth', 3, 'LineStyle', '--', 'Color','k');

% statistical plot
figure; hold on;
tmp = gcf; tmp.Position(3) = 270;
b = bar(1, mean(delay), 'BarWidth', 0.5, 'FaceAlpha',0.5);
b.FaceColor = 'c';
errorbar(1, mean(delay), 0, std(delay), 'LineWidth', 2, 'Color', 'k', 'LineStyle','none');
scatter(ones(N,1)+0.5*(rand(N,1)-0.5), delay, 50, 'w', 'filled', 'MarkerEdgeColor','k');
xlim([0 2])
ylim([0 8]);
yticks([0 4 8]);
yticklabels({'0', '4', '8'});

%% BP plot
% time courses
figure; hold on;
t = (-20:0.5:80);
plot(t, bp, 'LineWidth',1, 'Color',0.5*[1 1 1]);
plot(t, mean(bp,2), 'LineWidth',5, 'Color', c(1,:));
% min-max change
figure; hold on;
tmp = gcf; tmp.Position(3) = 160;
b = bar([1,2], mean([minBP,maxBP],1), 'BarWidth', 0.5, 'FaceAlpha',0.5);
b.FaceColor = c(1,:);
errorbar([1,2], mean([minBP,maxBP],1), [std(minBP) 0], [0 std(maxBP)], 'LineWidth', 2, 'Color', 'k', 'LineStyle','none');
x1 = ones(N,1)+0.3*(rand(N,1)-0.5);
x2 = 2*ones(N,1)+0.3*(rand(N,1)-0.5);
plot([x1';x2'], [minBP';maxBP'], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(x1, minBP, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(x2, maxBP, 50, 'w', 'filled', 'MarkerEdgeColor','k');
xlim([0.5 2.5]);
ylim([-40 5]);

% prediction accuracy with single-neurons
    % figure; hold on;
    % b = bar((1:N), maccubp);
    % b.FaceColor = c(1,:);
    % errorbar((1:N), maccubp, 0*stdaccubp, stdaccubp, 'Color', 'k', 'LineStyle', 'none');
    % for i = 1:N
    %     scatter(i*ones(numel(data.accuSing{i}),1)+0.2*(rand(numel(data.accuSing{i}),1)-0.5), data.accuSing{i}, 20, 'w', 'filled', 'MarkerEdgeColor','k');
    % end
    % ylim([0 1]);

%% HR plot
figure; hold on;
t = (-20:0.5:80);
plot(t, hr, 'LineWidth',1, 'Color',0.5*[1 1 1]);
plot(t, mean(hr,2), 'LineWidth',5, 'Color', c(2,:));
% min-max change
figure; hold on;
tmp = gcf; tmp.Position(3) = 160;
b = bar([1,2], mean([minHR,maxHR],1), 'BarWidth', 0.5, 'FaceAlpha',0.5);
b.FaceColor = c(2,:);
errorbar([1,2], mean([minHR,maxHR],1), [std(minHR) 0], [0 std(maxHR)], 'LineWidth', 2, 'Color', 'k', 'LineStyle','none');
x1 = ones(N,1)+0.3*(rand(N,1)-0.5);
x2 = 2*ones(N,1)+0.3*(rand(N,1)-0.5);
plot([x1';x2'], [minHR';maxHR'], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(x1, minHR, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(x2, maxHR, 50, 'w', 'filled', 'MarkerEdgeColor','k');
xlim([0.5 2.5]);
ylim([-20 10]);

% figure; hold on;
% avg = maccuhr;
% st = stdaccuhr;
% b = bar((1:N), avg);
% b.FaceColor = c(2,:);
% errorbar((1:N), avg, 0*st, st, 'Color', 'k', 'LineStyle', 'none');
% for i = 1:N
%     cellN = numel(data.accuSingHR{i});
%     scatter(i*ones(cellN,1)+0.2*(rand(cellN,1)-0.5), data.accuSingHR{i}, 20, 'w', 'filled', 'MarkerEdgeColor','k');
% end
% ylim([0 1]);
%% HRV plot
% low-freq.
figure; hold on;
plot(t, hrvLFpw, 'LineWidth',1, 'Color',0.5*[1 1 1]);
plot(t, mean(hrvLFpw,2), 'LineWidth',5, 'Color', c(3,:));
% min-max change
figure; hold on;
tmp = gcf; tmp.Position(3) = 160;
b = bar([1,2], mean([minLF,maxLF],1), 'BarWidth', 0.5, 'FaceAlpha',0.5);
b.FaceColor = c(3,:);
errorbar([1,2], mean([minLF,maxLF],1), [std(minLF) 0], [0 std(maxLF)], 'LineWidth', 2, 'Color', 'k', 'LineStyle','none');
x1 = ones(N,1)+0.3*(rand(N,1)-0.5);
x2 = 2*ones(N,1)+0.3*(rand(N,1)-0.5);
plot([x1';x2'], [minLF';maxLF'], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(x1, minLF, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(x2, maxLF, 50, 'w', 'filled', 'MarkerEdgeColor','k');
xlim([0.5 2.5]);
% ylim([-20 10]);

% figure; hold on;
% avg = maccuhrvLFpw;
% st = stdaccuhrvLFpw;
% b = bar((1:N), avg);
% b.FaceColor = c(3,:);
% errorbar((1:N), avg, 0*st, st, 'Color', 'k', 'LineStyle', 'none');
% for i = 1:N
%     cellN = numel(data.accuSingHRVlfpw{i});
%     scatter(i*ones(cellN,1)+0.2*(rand(cellN,1)-0.5), data.accuSingHRVlfpw{i}, 20, 'w', 'filled', 'MarkerEdgeColor','k');
% end
% ylim([0 1]);

% high-freq.
figure; hold on;

plot(t, hrvHFpw, 'LineWidth',1, 'Color',0.5*[1 1 1]);
plot(t, mean(hrvHFpw,2), 'LineWidth',5, 'Color', 'y');
% min-max change
figure; hold on;
tmp = gcf; tmp.Position(3) = 160;
b = bar([1,2], mean([minHF,maxHF],1), 'BarWidth', 0.5, 'FaceAlpha',0.5);
b.FaceColor = 'y';
errorbar([1,2], mean([minHF,maxHF],1), [std(minHF) 0], [0 std(maxHF)], 'LineWidth', 2, 'Color', 'k', 'LineStyle','none');
x1 = ones(N,1)+0.3*(rand(N,1)-0.5);
x2 = 2*ones(N,1)+0.3*(rand(N,1)-0.5);
plot([x1';x2'], [minHF';maxHF'], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(x1, minHF, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(x2, maxHF, 50, 'w', 'filled', 'MarkerEdgeColor','k');
xlim([0.5 2.5]);
% ylim([-20 10]);

% figure; hold on;
% avg = maccuhrvHFpw;
% st = stdaccuhrvHFpw;
% b = bar((1:N), avg);
% b.FaceColor = 'y';
% errorbar((1:N), avg, 0*st, st, 'Color', 'k', 'LineStyle', 'none');
% for i = 1:N
%     cellN = numel(data.accuSingHRVhfpw{i});
%     scatter(i*ones(cellN,1)+0.2*(rand(cellN,1)-0.5), data.accuSingHRVhfpw{i}, 20, 'w', 'filled', 'MarkerEdgeColor','k');
% end
% ylim([0 1]);

%% Accuracy comparison
figure;
ax = axes();
hold(ax);
boxchart(ones(N,1), maccuhrvLFpw,...
    'BoxFaceColor', c(3,:), 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(N,1), maccuhrvHFpw,...
    'BoxFaceColor', 'y', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(3*ones(N,1), maccuhr,...
    'BoxFaceColor', c(2,:), 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(4*ones(N,1), maccubp,...
    'BoxFaceColor', c(1,:), 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10); 
plot(repmat([1;2;3;4],1,N), [maccuhrvLFpw;maccuhrvHFpw;maccuhr;maccubp], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(N,1), maccuhrvLFpw, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(2*ones(N,1), maccuhrvHFpw, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(3*ones(N,1), maccuhr, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(4*ones(N,1), maccubp, 50, 'w', 'filled', 'MarkerEdgeColor','k');
xlim([0.5 4.5]);
ylim([0 1]);

disp(['(single) hrv LF vs hrv HF: ', num2str(ranksum(maccuhrvLFpw, maccuhrvHFpw))]);
disp(['(single) hrv HF vs hr: ', num2str(ranksum(maccuhrvHFpw, maccuhr))]);
disp(['(single) bp vs hr: ', num2str(ranksum(maccubp, maccuhr))]); 

%%
figure;
ax = axes();
hold(ax);
decbp = data.accuIndv;
dechr = data.accuIndvHR;
dechf = data.accuIndvHRVhfpw;
declf = data.accuIndvHRVlfpw;

boxchart(ones(N,1), declf,...
    'BoxFaceColor', c(3,:), 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(N,1), dechf,...
    'BoxFaceColor', 'y', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(3*ones(N,1), dechr,...
    'BoxFaceColor', c(2,:), 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(4*ones(N,1), decbp,...
    'BoxFaceColor', c(1,:), 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10); 
plot(repmat([1;2;3;4],1,N), [declf;dechf;dechr;decbp], 'Color', [0 0 0 0.2], 'LineWidth', 1.5);
scatter(ones(N,1), declf, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(2*ones(N,1), dechf, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(3*ones(N,1), dechr, 50, 'w', 'filled', 'MarkerEdgeColor','k');
scatter(4*ones(N,1), decbp, 50, 'w', 'filled', 'MarkerEdgeColor','k');
xlim([0.5 4.5]);
ylim([0 1]);

disp(['(decoding) hrv LF vs hrv HF: ', num2str(ranksum(declf, dechf))]);
disp(['(decoding) hrv HF vs hr: ', num2str(ranksum(dechf, dechr))]);
disp(['(decoding) bp vs hr: ', num2str(ranksum(dechr, decbp))]); 

return

%%
figure;
subplot(2,1,1); hold on;
errorbar(mean(dlfpw), std(dlfpw), 'Color', 'k', 'LineStyle', 'none');
bar(1, mean(dlfpw), 'FaceColor','r');
scatter(ones(1,10)+0.2*(rand(1,10)-0.5), dlfpw, 50, 'w', 'filled', 'MarkerEdgeColor','k');
subplot(2,1,2); hold on;
errorbar(mean(dhfpw), std(dhfpw), 'Color', 'k', 'LineStyle', 'none');
bar(1, mean(dhfpw), 'FaceColor','g');
scatter(ones(1,10)+0.2*(rand(1,10)-0.5), dhfpw, 50, 'w', 'filled', 'MarkerEdgeColor','k');


figure;hold on;
errorbar(mean(dbp), std(dbp), 'Color', 'k', 'LineStyle', 'none');
bar(1, mean(dbp), 'FaceColor',c(1,:));
scatter(ones(1,10)+0.2*(rand(1,10)-0.5), dbp, 50, 'w', 'filled', 'MarkerEdgeColor','k');

figure;hold on; dhr([5 8]) = [];
errorbar(mean(dhr), std(dhr), 'Color', 'k', 'LineStyle', 'none');
bar(1, mean(dhr), 'FaceColor',c(2,:));
scatter(ones(1,8)+0.2*(rand(1,8)-0.5), dhr, 50, 'w', 'filled', 'MarkerEdgeColor','k');

%%
tp = tp-20;
% figure;hold on;
% plot(tp, (lfpw)', 'LineWidth',1, 'Color',0.5*[1 1 1]);
% plot(tp, mean(lfpw,1), 'LineWidth',3, 'Color', 'r');

%%
figure;
subplot(2,1,1); hold on;
plot(tp, (lfpw)', 'LineWidth',1, 'Color',0.5*[1 1 1]);
plot(tp, mean(lfpw,1), 'LineWidth',3, 'Color', 'r');
subplot(2,1,2); hold on;
plot(tp, (hfpw)', 'LineWidth',1, 'Color',0.5*[1 1 1]);
plot(tp, mean(hfpw,1), 'LineWidth',3, 'Color', 'r');
figure; hold on;
plot(tp, ratio', 'LineWidth',1, 'Color',0.5*[1 1 1]);
plot(tp, mean(ratio,1), 'LineWidth',3, 'Color', 'r');
%%

%%
figure;hold on;
errorbar(mean(delay), std(delay), 'Color', 'k', 'LineStyle', 'none');
bar(1, mean(delay));
scatter(ones(1,10)+0.7*(rand(1,10)-0.5), delay, 50, 'w', 'filled', 'MarkerEdgeColor','k');

