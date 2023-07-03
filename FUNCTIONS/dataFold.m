function out = dataFold(data, rowN, colN)
    out = zeros(rowN, colN);
    for i = 1:colN
        out(:,i) = data(1+(i-1)*rowN:i*rowN);
    end
end