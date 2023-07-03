function DATA = tdt2klusta(data)

Fs = double(typecast(uint32(1186905120), 'single'));
dt = 1/Fs;
kernel = gausswin(round(0.002/dt));
kernel = kernel/sum(kernel);

%% artifact remove

[num, den] = butter(3, 20/(Fs*0.5), 'high');

for i = 1:16
    TEMP = double(data(i,:));
    TEMPTEMP = filtfilt(num, den, TEMP);
    FILT = filtfilt(kernel, 1, abs(TEMPTEMP));
    RMS = mean(abs(FILT).^2).^0.5;
    STD = std(abs(FILT));
    minpeak = (max(FILT)+min(FILT))*0.5;
    
    [pks, locs] = findpeaks(FILT, Fs, 'MinPeakHeight', minpeak, 'MinPeakDistance', 0.0005);
    
    locs_idx = round(locs*Fs) + 1;
    
    %%%%%%%% remove the artifact tails (template substraction) %%%%%%%%
    src = (20:1:300)';
    
    platform = repmat(src, 1, length(locs_idx));

    locs_idx_platform = repmat(locs_idx, length(src), 1) + platform;
    try
        stim_artifact = TEMP(locs_idx_platform);
    catch
        fprintf('artifact peak detection was wrong \n');
%         figure;
%         plot(TEMP);
%         yline(minpeak);
%         hold on; 
%         scatter(locs*fs, pks);
        break;
    end

    stim_artifact = reshape(stim_artifact, length(src), length(locs_idx));
    mean_stim_artifact = mean(stim_artifact, 2);
    TEMP(locs_idx_platform) = TEMP(locs_idx_platform) - repmat(mean_stim_artifact, 1, length(locs_idx));
    
   %% Discard and copy
    neg = -30;
    pos = 30;

    src_neg = (neg:1:0)';
    src_pos = (1:1:pos)';

    platform_neg = repmat( src_neg, 1, length(locs_idx));
    platform_pos = repmat( src_pos, 1, length(locs_idx));

    locs_idx_platform_neg = repmat(locs_idx, length(src_neg), 1) + platform_neg;
    locs_idx_platform_pos = repmat(locs_idx, length(src_pos), 1) + platform_pos;

    TEMP(locs_idx_platform_neg) = TEMP(locs_idx_platform_neg + neg);
    TEMP(locs_idx_platform_pos) = TEMP(locs_idx_platform_pos + pos);
    
    data(i,:) = TEMP;
    
    clear TEMP;
    
end

%figure;plot(data.streams.R_NT.data(1,:)); hold on; plot(data1.streams.R_NT.data(1,:));

DATA = data;

