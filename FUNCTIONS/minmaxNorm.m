function out = minmaxNorm(data, range)
    minData = min(data);
    maxData = max(data);
    minR = range(1);
    maxR = range(2);
    out = (data-minData)./(maxData-minData)*(maxR-minR) + minR;
end