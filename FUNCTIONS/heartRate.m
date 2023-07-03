function [raw_HR, filt_HR] = heartRate(BP_TIME, BP_RAW)

moving_avg_window = 300;
moving_avg1 = 1/moving_avg_window*ones(moving_avg_window,1);

[LF_num, LF_den] = butter(2, 0.5/(100*0.5), 'low');
[HF_num, HF_den] = butter(2, 1/(100*0.5), 'High');


FILT = filtfilt(HF_num, HF_den, BP_RAW);
%plot(b2(:,1), filter(LF_num, LF_den, filter(HF_num, HF_den, b2(:,2))));
[~, locs] = findpeaks(FILT, BP_TIME, 'MinPeakDistance', 0.1);

test = abs(hilbert(FILT));

%scatter(locs, pks);
%[pks, locs] = findpeaks(filtfilt(HF_num, HF_den, b2(:,2)), b2(:,1), 'MinPeakDistance', 0.12);
%[pks, locs] = findpeaks(filter(LF_num, LF_den, filter(HF_num, HF_den, b2(:,2))), b2(:,1), 'MinPeakDistance', 0.12);
%scatter(locs, pks);

EDGES = (BP_TIME(1:end-1)+BP_TIME(2:end))*0.5;
EDGES = [2*BP_TIME(1)-EDGES(1); EDGES; 2*BP_TIME(end)-EDGES(end)];
[COUNTS, ~] = histcounts(locs, EDGES);
HR = COUNTS'./(EDGES(2:end) - EDGES(1:end-1))*60;

%%
raw_HR = HR;
% filt_HR = filtfilt(LF_num, LF_den, HR);
filt_HR = filtfilt(moving_avg1, 1, HR);

% padding to reject edge error
pad = moving_avg_window;
if size(filt_HR,1) < size(filt_HR,2)
    filt_HR(1:pad) = fliplr(filt_HR(pad+1:2*pad));
    filt_HR(end-pad+1:end) = fliplr(filt_HR(end-2*pad+1:end-pad));
else
    filt_HR(1:pad) = flipud(filt_HR(pad+1:2*pad));
    filt_HR(end-pad+1:end) = flipud(filt_HR(end-2*pad+1:end-pad));
end
end