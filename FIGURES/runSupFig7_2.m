clear;
close all;
clc;

addpath('..\FUNCTIONS\');

params = load('..\PARAMETER_DATA\params.mat');

%%
i = 1;
dt = params.DT(i);
% Normalized neural trajectory
data = load(['..\POSTPROCESSED_MAT_CODE\DATA\', params.rat{i}, '_matlab_save.mat']);
[t, bp, ~, fr,iso1, iso2] = dataPrepare(data, dt, params.TRIAL(i), params.ntsDELAY(i), params.bpDELAY(i), params.PRE(i));
t(:,1) = [];
frR = zeros(size(t,1), numel(fr));
for j = 1:numel(fr)
    frR(:,j) = mean(fr{j}(:,2:end),2);
end
bp(:,1) = [];
iso1(:,1) = [];
iso2(:,1) = [];

% set time range
t = t(:,1);
% Get aligned latent space
[b, ~, ~, ~] = linRegress([ones(size(iso1,1),1), mean(iso1,2), mean(iso2,2)], mean(bp,2));
[~, algn1,algn2] = alignLatentSlide(iso1, iso2, 2);
[~, algn1,algn2] = alignLatentScale(algn1, algn2, 2);
[~, algn1,algn2] = alignLatentRotate(algn1, algn2, b);

tra1 = mean(algn1,2);
tra2 = mean(algn2,2);
% i = 1;
% exp = load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\expDataAnalysis.mat');
% fr = exp.fr{i};
% t = exp.t{i};

[~, mxIdx] = max(frR, [], 1);
mxIdx = mod(mxIdx,numel(t));
[~, srtIdx] = sort(mxIdx);
frR = frR(:,srtIdx);
figure; hold on;
for j = 1:numel(srtIdx)
    frR(:,j) = minmaxNorm(frR(:,j), [-1 1]);
    imagesc(t',j*ones(size(t')), frR(:,j)');
end
ylim([0.5 size(frR,2)+0.5]);
xlim([t(1) t(end)]);
xline(0, 'Color', 'k', 'LineWidth',3);
xline(60, 'Color', 'k', 'LineWidth',3);
colormap viridis(100);

f = figure; hold on; f.Position(4) = 180;
repidx = [3 11 19 24];
plot(t, frR(:,repidx)' - mean(frR(1:40,repidx),1)', 'LineWidth', 3);
xlim([t(1) t(end)]);

%%
rnnData = load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\LowRankRNNTestRepresentative.mat');
rnn = rnnData.TrainedRNN;
actF = @(x) tanh(x);
trainInfo = rnnData.trainInfo;
[x, ~] = hBINDsimulate(trainInfo.N, trainInfo.tau, trainInfo.tau, 0, zeros(trainInfo.N), actF, t,...
    rnn.z, rnn.u, rnn.Wz, rnn.Wu, rnn.Wr);

frRNN = actF(x);
[~, mxIdx] = max(frRNN, [], 2);
mxIdx = mod(mxIdx,numel(t));
[~, srtIdx] = sort(mxIdx);
frRNN = frRNN(srtIdx,:);

figure; hold on;
for i = 1:numel(srtIdx)
    frRNN(i,:) = minmaxNorm(frRNN(i,:), [-1 1]);
    imagesc(t, i*ones(size(t)), frRNN(i,:));
end
ylim([0.5 size(frRNN,1)+0.5]);
xlim([t(1) t(end)]);
xline(0, 'Color', 'k', 'LineWidth',3);
xline(60, 'Color', 'k', 'LineWidth',3);
colormap viridis(100);

%
closest = zeros(1, numel(repidx));
for i = 1:numel(repidx)
    accuracyTmp = zeros(trainInfo.N,1);
    for j = 1:trainInfo.N
        tmp = corrcoef(frR(:,repidx(i)), frRNN(j,:));
        accuracyTmp(j) = sign(tmp(1,2))*tmp(1,2)^2;
    end
    [~, closest(i)] = max(accuracyTmp);
end
f = figure; hold on; f.Position(4) = 180;
plot(t, frRNN(closest,:) - mean(frRNN(closest,1:40),2), 'LineWidth', 3);
xlim([t(1) t(end)]);

%%
figure;
scatter(tra1, tra2, 80, t, 'filled', 'MarkerEdgeColor', 'k');
axis equal;
xlim([-1 1]);
ylim([-1 1]);
colormap hot;

figure;
scatter(rnn.Wr(:,1)'*actF(x), rnn.Wr(:,2)'*actF(x), 80, t, 'filled', 'MarkerEdgeColor', 'k');
axis equal;
xlim([-1 1]);
ylim([-1 1]);
colormap hot;

%%
figure;hold on;
plot(t, tra1, 'LineWidth', 3', 'Color', 'k');
plot(t, rnn.Wr(:,1)'*actF(x), 'LineWidth', 3', 'Color', 'r');
xlim([t(1) t(end)]);
ylim([-1 1]);
colormap hot;
tmp = corrcoef(tra1, rnn.Wr(:,1)'*actF(x));
disp(tmp(1.2)^2)

figure;hold on;
plot(t, tra2, 'LineWidth', 3', 'Color', 'k');
plot(t, rnn.Wr(:,2)'*actF(x), 'LineWidth', 3', 'Color', 'r');
xlim([t(1) t(end)]);
ylim([-1 1]);
colormap hot;
tmp = corrcoef(tra2, rnn.Wr(:,2)'*actF(x));
disp(tmp(1.2)^2)