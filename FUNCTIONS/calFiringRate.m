function fr = calFiringRate(spkTime, edges, gaussianParams)

dt = mean(edges(2:end)-edges(1:end-1));
sig = gaussianParams.sigma;
L = gaussianParams.length;
x = dt*((0:L-1)- floor(L*0.5));
gWin = exp(-(x.^2) / (2 * sig^2)) / (sig * sqrt(2 * pi));

tmp = histcounts(spkTime, edges);
fr = filtfilt(gWin, sum(gWin), tmp);

errIdx = tmp<0;
fr(errIdx) = 0;

end