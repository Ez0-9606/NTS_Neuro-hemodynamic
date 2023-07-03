function out = weightedStd(data, weight)
    if size(weight,2) ~= 1
        error('weightedMean error:: weight should be n x 1 vector');
    end

    stdTmp = std(data,[],1);
    out = stdTmp*sum(weight.^2)/sum(weight)^2;
end