function hrDataClip = hrClip(hr, preTime, postTime, stimOnTime)

s = stimOnTime - preTime;
e = stimOnTime + postTime;

tHR = (0:size(hr,1)-1)*0.01;

tmp = [];
for i = 1:numel(stimOnTime)

    if s(i) < 0
        error('hrClip error the preTime is too short\n');
    end
    if e(i) > tHR(end)
        error('hrClip error the postTime is too long\n');
    end

    idx = (tHR < e(i)) & (tHR >= s(i));
    tmp = [tmp; hr(idx,:)];
end
hrDataClip = tmp;
%hrDataClip = (hrDataClip(:,3));

end