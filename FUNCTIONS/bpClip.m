function bpDataClip = bpClip(bp, preTime, postTime, stimOnTime)

s = stimOnTime - preTime;
e = stimOnTime + postTime;

tBP = (0:size(bp,1)-1)*0.01;

tmp = [];
for i = 1:numel(stimOnTime)

    if s(i) < 0
        error('bpClip error the preTime is too short\n');
    end
    if e(i) > tBP(end)
        error('bpClip error the postTime is too long\n');
    end

    idx = (tBP < e(i)) & (tBP >= s(i));
    tmp = [tmp; bp(idx,:)];
end
bpDataClip = tmp;
%bpDataClip = (bpDataClip(:,4));

end