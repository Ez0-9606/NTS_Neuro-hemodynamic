function [ntsData, stimData, tNTS, stimOnTime, fig] = preprcssNTS(ratPath, rawDir, fs, ntsCh, stimCh, artrmInfo)
    
    % Load data
    [ntsData, stimData] = loadNTSraw(rawDir, ntsCh, stimCh);
    
    % Define time axis
    tNTS = getTimeAxis(size(ntsData,2), fs);
    
    % Re-arrange
    ntsData = reArrange16(ntsData);
    
    % remove artifact
    ntsDataTmp = ntsData;
    [~, stimOnTime] = findpeaks(stimData(1,:), fs, 'MinPeakHeight', 10^4, 'MinPeakDistance', artrmInfo.trialInterval);
    for i = 1:numel(stimOnTime)
        [~, idx1] = min(abs(tNTS - stimOnTime(i) + artrmInfo.preDur));
        [~, idx2] = min(abs(tNTS - stimOnTime(i) - artrmInfo.postDur));
        try
            ntsData(:,idx1:idx2) = tdt2klusta(ntsData(:,idx1:idx2));
        catch
            id = fopen('errLog.txt', 'a');
            fprintf(id, '%s', ['Stimulus artifact removal error for ', ratNumb, '\n']); fclose(id);
            break;
        end
    end

    % verify artifact removal
    fig = figure;
    subplot(2,1,1);
    hold on;
    plot(downsample(tNTS, 10), downsample(ntsDataTmp(1,:),10));
    plot(downsample(tNTS, 10), downsample(ntsData(1,:),10));
    for i = 1:numel(stimOnTime)
        [~, idx] = min(abs(tNTS - stimOnTime(i)));
        xline(tNTS(idx));
    end
    
    subplot(2,1,2);
    hold on;
    t = seconds(tNTS);
    t.Format = 'hh:mm:ss';
    plot(downsample(t, 10), downsample(ntsDataTmp(1,:),10));
    plot(downsample(t, 10), downsample(ntsData(1,:),10));
    for i = 1:numel(stimOnTime)
        tmp = seconds(stimOnTime(i));
        tmp.Format = 'hh:mm:ss';
        [~, idx] = min(abs(t - tmp));
        xline(t(idx));
    end

    % save NTS recording data
    fileID = fopen([ratPath, 'NTS_preprocessed.dat'],'w');
    fwrite(fileID, ntsData,'int16');
    fclose(fileID);
    saveas(fig, [ratPath, 'NTS_timeline.fig']);
end