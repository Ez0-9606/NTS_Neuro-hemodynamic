function hrv = heartRateVar(tHR, hr, fs, dt)
    
    % HRV calculation
    Fs = fs;
    window = round(dt*Fs);
    g = hann(window,"periodic");
    overlapRatio = 0.5;% 80 percent overlap of window
    noverlap = round(overlapRatio*dt*Fs); 
    nfft = max(256,2^ceil(log2(window)));
    
    [~, fHRV, tTmp, psp] = spectrogram(hr, g, noverlap, nfft, Fs, 'power');
    flimit = 10; %[Hz]
    n = round(nfft/Fs*flimit);
    fHRV = fHRV(1:n);
    psp = psp(1:n,:); % unit [second^2]
    pwHRV = psp;

    % HRV data interpolation
    % get temporal pre-interpolation data array (time & power data)
    preInterpIdx = (0.5*dt*Fs+1:(1-overlapRatio)*window:0.5*dt*Fs+1+(1-overlapRatio)*window*(numel(tTmp)-1));
    postInterpIdx = (1:numel(tHR));
    %tInterpFun = griddedInterpolant(preInterpIdx, tTmp);
    %tHRV = tInterpFun(postInterpIdx);
    %%%%% !!! tHRV is equal to tBP !!! (error < 10^-12)

    [F, T] = ndgrid((1:numel(fHRV)), preInterpIdx);
    [aF, aT] = ndgrid((1:numel(fHRV)), postInterpIdx);
    pwInterpFun = griddedInterpolant(F, T, psp, 'linear', 'nearest');
    pwHRV = pwInterpFun(aF, aT);
    pwHRV(pwHRV<0) = 0;

    hrv = struct();
    hrv.fHRV = fHRV;
    hrv.tHRV = tHR;
    hrv.pwHRV = pwHRV;
end