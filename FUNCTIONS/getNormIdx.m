function idx = getNormIdx(inputArray, referenceArray)
    idx = zeros(size(referenceArray));
    for i = 1:numel(referenceArray)
        [~, idx(i)] = min(abs(inputArray-referenceArray(i)));
    end
end