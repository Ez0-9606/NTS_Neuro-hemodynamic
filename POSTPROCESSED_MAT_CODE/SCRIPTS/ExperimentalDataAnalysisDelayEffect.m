clear;
close all;
clc;

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
TRIAL = [4 3 3,...
    4 3 3 4,...
    4 4 3];
PRE = [-20 -20 -20,...
    -20 -20 -20 -20,...
    -20 -20 -40];
POST = [90 90 90,...
    120 120 120 120,...
    120 100 160];

f = figure; f.Position = [1, 41, 1920, 962];
for i = 1:numel(rat)
    if i == 9
        c = [1 0 0];
    else
        c = [0 0 0];
    end
    %i = numel(rat) - ii + 1;
    ratNumb = rat{i};
    data = load(['..\DATA\', ratNumb, '_matlab_save.mat']);

    dt = DT(i);
    tr = TRIAL(i);
    bDelay = bpDELAY(i);
    nDelay = ntsDELAY(i);
    pre = PRE(i);
    %%
    [t, bp, ~, iso1, iso2] = dataPrepare(data, dt, tr, nDelay, bDelay, pre);
    
    %%
    [algn1,algn2] = alignLatentSlide(iso1, iso2, 2);
    [algn1,algn2] = alignLatentScale(algn1, algn2, 2);

%     hold off
    plot(t(:,1), -mean(algn2,2), 'g');
    [~, idx] = min(abs(t(:,1)));
    xline(0);
    hold on;
    mBP = mean(bp,2);
    plot(t(:,1), mBP - mean(mBP(1:idx)),'Color',c);
    
end
xlim([-10 15]);