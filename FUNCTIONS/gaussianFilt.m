function sigFilt = gaussianFilt(data, dt, sigma, length)
    sig = sigma;
    L = length;
    x = dt*((0:L-1)- floor(L*0.5));
    gWin = exp(-(x.^2) / (2 * sig^2)) / (sig * sqrt(2 * pi));
    sigFilt = filtfilt(gWin, sum(gWin), data);
    
end