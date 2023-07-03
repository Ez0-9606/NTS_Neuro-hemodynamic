function hrvDataClip = hrvClip(hrv, preTime, postTime, stimOnTime)

s = stimOnTime - preTime;
e = stimOnTime + postTime;

tHRV = hrv.tHRV;

pwTmp = [];
tTmp = [];
for i = 1:numel(stimOnTime)

%     if s(i) < 0
%         error('hrClip error the preTime is too short\n');
%     end
%     if e(i) > tHRV(end)
%         error('hrClip error the postTime is too long\n');
%     end

    idx = (tHRV < e(i)) & (tHRV >= s(i));
    tTmp = [tTmp, hrv.tHRV(idx)];
    pwTmp = [pwTmp, hrv.pwHRV(:,idx)];

end
hrvDataClip.pwHRV = pwTmp;
hrvDataClip.fHRV = hrv.fHRV;
hrvDataClip.tHRV = tTmp;

end