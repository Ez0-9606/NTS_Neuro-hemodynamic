function out = weightedMean(data, weight)
    if size(weight,2) ~= 1
        error('weightedMean error:: weight should be n x 1 vector');
    end
    out = sum(data.*weight,1);
end