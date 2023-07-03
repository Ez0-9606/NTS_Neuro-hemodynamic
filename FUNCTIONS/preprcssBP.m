function [bp, hr, hrv, tBP, figBP] = preprcssBP(ratPath, rawDir, fs)
    % Load data
    bpAddrs = findBPfile(rawDir);
    bpAddrs = bpAddrs.addrs;
    d = load(bpAddrs);
    bp = zeros(numel(d.b1(:,2)),2);
%     bp(:,1) = double(d.b1(:,2));
    
    % Define time axis
    tBP = getTimeAxis(size(bp,1), fs);
    
    % High pass filtering
    [HF_num, HF_den] = butter(2, 0.001/(fs*0.5), 'high');
    [LF_num, LF_den] = butter(2, 0.5/(fs*0.5), 'low');
%     bp(:,2) = filtfilt(LF_num, LF_den, bp(:,1));
%     bp(:,3) = filtfilt(HF_num, HF_den, bp(:,1));
%     bp(:,4) = filtfilt(LF_num, LF_den, bp(:,3));
    bp(:,1) = filtfilt(HF_num, HF_den, double(d.b1(:,2)));
    bp(:,2) = filtfilt(LF_num, LF_den, bp(:,1));

    % HR calculation
    [~, HR] = heartRate(tBP', bp(:,1));
%     hr(:,1) = HR;
%     hr(:,2) = filtfilt(LF_num, LF_den, hr(:,1));
%     hr(:,3) = filtfilt(HF_num, HF_den, hr(:,1));
    hr(:,1) = filtfilt(HF_num, HF_den, HR);
    hr(:,2) = filtfilt(LF_num, LF_den, hr(:,1));

    plot(bp)
    hold on;
    plot(hr)

    % HRV calculation
    Fs = 100;
    dt = 15; %[s]
    window = round(dt*Fs);
    g = hann(window,"periodic");
    overlapRatio = 0.8;% 80 percent overlap of window
    noverlap = round(overlapRatio*dt*Fs); 
    nfft = max(256,2^ceil(log2(window)));
    
    [~, fHRV, tTmp, psp] = spectrogram(hr(:,1)', g, noverlap, nfft, Fs, 'power');
    flimit = 3; %[Hz]
    n = round(nfft/Fs*flimit);
    fHRV = fHRV(1:n);
    psp = psp(1:n,:); % unit [second^2]

    % HRV data interpolation
    % get temporal pre-interpolation data array (time & power data)
    preInterpIdx = (0.5*dt*Fs+1:(1-overlapRatio)*window:0.5*dt*Fs+1+(1-overlapRatio)*window*(numel(tTmp)-1));
    postInterpIdx = (1:numel(tBP));
    %tInterpFun = griddedInterpolant(preInterpIdx, tTmp);
    %tHRV = tInterpFun(postInterpIdx);
    %%%%% !!! tHRV is equal to tBP !!! (error < 10^-12)

    [F, T] = ndgrid((1:numel(fHRV)), preInterpIdx);
    [aF, aT] = ndgrid((1:numel(fHRV)), postInterpIdx);
    pwInterpFun = griddedInterpolant(F, T, psp);
    pwHRV = pwInterpFun(aF, aT);
    
    hrv = struct();
    hrv.fHRV = fHRV;
    hrv.tHRV = tBP;
    hrv.pwHRV = pwHRV;

    % Save BP, HR data
    save([ratPath, 'BP_preprocessed.mat'], 'bp', 'hr', 'hrv');

    figBP = 0;
%     figBP = figure; hold on;
%     plot(tBP, bp);
%     tmp = fieldnames(d);
%     for i = 1:numel(tmp)-2
%         xline(d.(tmp{i+2}).time, '-', {d.(tmp{i+2}).value});
%     end
%     saveas(figBP, [ratPath, 'BP_timeline.fig']);

end