function [spkTimeClip, spkWaveClip] = spkClip(spkTimeIn, spkWaveIn, preTime, postTime, stimOnTime)

    s = stimOnTime - preTime;
    e = stimOnTime + postTime;
    L = postTime + preTime;
    
    cellNumb = numel(spkTimeIn);
    
    spkTimeClip = cell(1,cellNumb);
    for i = 1:cellNumb
        tmp1 = [];
        tmp2 = [];
        for j = 1:numel(stimOnTime)
            idx = (spkTimeIn{i} < e(j)) & (spkTimeIn{i} >= s(j));
            tmp1 = [tmp1; spkTimeIn{i}(idx) - s(j) + L*(j-1)];
            tmp2 = [tmp2; spkWaveIn{i}(idx,:)];
        end
        spkTimeClip{i} = tmp1;
        spkWaveClip{i} = tmp2;
    end

end