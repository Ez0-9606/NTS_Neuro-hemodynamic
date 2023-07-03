function [bp, hr, hrv, tBP, figBP] = preprcssBPwithoutLoad(bpR,ratPath, fs)
    % Define time axis
    tBP = getTimeAxis(size(bpR,1), fs);
    
    % High pass filtering
    [HF_num, HF_den] = butter(2, 0.004/(fs*0.5), 'high');
    [LF_num, LF_den] = butter(2, 0.5/(fs*0.5), 'low');
    %bp(:,1) = filtfilt(HF_num, HF_den, bpR);
    bp(:,1) = bpR;
    bp(:,2) = filtfilt(LF_num, LF_den, bpR);
    
    % HR calculation
    [~, HR] = heartRate(tBP', bp(:,1));
    hr(:,1) = filtfilt(HF_num, HF_den, HR);
    hr(:,2) = filtfilt(LF_num, LF_den, hr(:,1));


    plot(bpR)
    hold on;
    plot(hr)

    % HRV calculation
    Fs = 100;
    dt = 15; %[s]
    window = round(dt*Fs);
    g = hann(window,"periodic");
    overlapRatio = 0.8;% 80 percent overlap of window
    noverlap = round(overlapRatio*dt*Fs); % 80 percent overlap of window
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
    pwHRV(pwHRV<0) = 10^-5;
    
    hrv = struct();
    hrv.fHRV = fHRV;
    hrv.tHRV = tBP;
    hrv.pwHRV = pwHRV;

    % Save BP, HR data
    save([ratPath, 'BP_preprocessed.mat'], 'bp', 'hr', 'hrv');
    
    figBP = 0;
%     figBP = figure; hold on;
%     plot(tBP, bp);
%     saveas(figBP, [ratPath, 'BP_timeline.fig']);
end