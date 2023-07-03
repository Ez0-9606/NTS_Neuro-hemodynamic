clear;
clc;
close all;

addpath('..\FUNCTIONS');
preprcssPath = '..\PREPROCESSED\';
params = load('..\PARAMETER_DATA\params.mat');

if ~exist(['.\', 'SPIKE_DATA'], 'dir')
    mkdir(['.\', 'SPIKE_DATA']);
end
if ~exist(['.\', 'BP_DATA'], 'dir')
    mkdir(['.\', 'BP_DATA']);
end
if ~exist(['.\', 'FR_DATA'], 'dir')
    mkdir(['.\', 'FR_DATA']);
end

for m = 1:numel(params.rat)
% neuron screening
% find only stimulation-driven neurons by using autocorrelation.
    % parameter initialization
    ratNumb = params.rat{m};
    dt = params.DT(m);
    tr = params.TRIAL(m);
    bDelay = params.bpDELAY(m);
    nDelay = params.ntsDELAY(m);
    pre = -params.PRE(m);
    post = params.POST(m);
    L = post + pre;
    
    plvTh = 0.1;
    gaussianParams = struct('sigma', dt, 'length', round(3/dt));

    % Data load
    if contains(ratNumb,'_')
        idx = strrep(ratNumb,'_', '\');
    else
        idx = ratNumb;
    end
    load([preprcssPath, idx, '\stimInfo']);
    load([preprcssPath, idx, '\spkTime']);
    load([preprcssPath, idx, '\spkWave']);
    load([preprcssPath, idx, '\BP_preprocessed']);

    % Spike data clipping
    [spkTimeClip, spkWaveClip] = spkClip(spkTime, spkWave, pre, post, stimOnTime);
    edges = (0:dt:L*numel(stimOnTime));
    t = 0.5*(edges(2:end)+edges(1:end-1));

    % Spike data screening
    frTmp = zeros(numel(spkTime), numel(t));
    for i = 1:numel(spkTime)
        frTmp(i,:) = calFiringRate(spkTimeClip{i}, edges, gaussianParams);
    end

    regSig = cos(2*pi*1/L*t);
    plv = zeros(numel(spkTime),1);
    for i = 1:numel(spkTime)
        plv(i) = calSigPLV(frTmp(i,:), regSig);
    end
    
    spike = spkTimeClip(plv > plvTh)'; 
    wave = spkWaveClip(plv > plvTh)'; 
    fr = frTmp(plv > plvTh, :);

    spikeOut = spkTimeClip(plv <= plvTh)';
    waveOut = spkWaveClip(plv <= plvTh)'; 
    frOut = frTmp(plv <= plvTh, :);


    % BP data clipping
    try
        BP = bpClip(bp, pre, post, stimOnTime-delay); % the 'delay' is measurment error!!!
        HR = hrClip(hr, pre, post, stimOnTime-delay);
        HRV = hrvClip(hrv, pre, post, stimOnTime-delay);
    catch
        fprintf(['Clip range error, rat number: ', ratNumb, ' (', num2str(m), ')\n']);
    end
    tBP = (0:size(BP,1)-1)*0.01;
    


%     figure;
%     subplot(2,1,1);
%     hold on;
%     plot(tBP, BP);
%     for n = 1:numel(stimOnTime)
%         xline(pre + (n-1)*L, '-', [num2str(stimFreq(n)), 'Hz, ', num2str(stimDur(n)), 's'], 'Color', 'r');
%         xline(stimDur(n) + pre + (n-1)*L, 'Color', 'k');
%     end
%     subplot(2,1,2);
%     hold on;
%     for i = 1:size(fr,1)
%         plot(t, fr(i,:));
%     end
%     for n = 1:numel(stimOnTime)
%         xline(pre + (n-1)*L, '-', [num2str(stimFreq(n)), 'Hz, ', num2str(stimDur(n)), 's'], 'Color', 'r');
%         xline(stimDur(n) + pre + (n-1)*L, 'Color', 'k');
%     end


    % Data save
%     save(['.\', 'SPIKE_DATA\', 'spike_data_', ratNumb, '_', 'v7.mat'], 'spike', 'spikeOut', '-v7');
%     save(['.\', 'SPIKE_DATA\', 'wave_data_', ratNumb, '_', 'v7.mat'], 'wave', 'waveOut', '-v7');
    save(['.\', 'BP_DATA\', 'bp_data_', ratNumb, '_', 'v7.mat'], 'BP', 'HR', 'HRV', '-v7');
%     save(['.\', 'FR_DATA\', 'fr_data_', ratNumb, '_', 'v7.mat'], 'fr', 'frOut', '-v7');

    % Re-initialize
    clearvars 'spike' 'spkTimeSave';
end


