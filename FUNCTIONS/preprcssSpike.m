function [spkTime, spkWave] = preprcssSpike(ratPath, fs, ntsData, tNTS)


    % load cluster info
    id = fopen([ratPath, '\cluster_info.tsv']);
    srt = textscan(id, '%s %s %s %s %s %s', 'Headerlines', 1);
    
    % load spike data
    name = [ratPath, '\NTS_preprocessed.kwik'];
    time_sample = double(h5read(name, '/channel_groups/0/spikes/time_samples'))/fs;
    cluster = double(h5read(name, '/channel_groups/0/spikes/clusters/main'));
    
    % structurize spike points and spike waveforms
    spkTime = {};
    spkWave = {};
    k = 1;
    for m = 1: numel(srt{5}) % numel(srt{5}): number of clusters
        if ~strcmp(cell2mat(srt{5}(m)), 'good') % select only for 'good' case
            continue;
        end
        idx = cluster == str2double(cell2mat(srt{1}(m))); 
        spike_points = time_sample(idx);% find spike points included in a cluster
    
        % exclude start&end spike points only if their waveforms (-16/fs ~ +16/fs)are excluded out of data range
        while(1)
            if (spike_points(end) + 16/fs) >  tNTS(end)
                spike_points(end) = [];
            else
                break;
            end
        end
    
        while(1)
            if (spike_points(1) - 16/fs) <  tNTS(1)
                spike_points(1) = [];
            else
                break;
            end
        end
    
        % calculate an averaged waveform to exclude amplitude jitters
        waveforms = ntsData(str2double(cell2mat(srt{2}(m)))+1,:);%ntsData(probe_channel_number,:)
        IDX = round(spike_points*fs) + 1;
        waveforms = waveforms(repmat((-16:1:16),numel(IDX),1) + repmat(IDX, 1, 33)); % get each spike waveform
        waveforms = waveforms - repmat(mean(waveforms, 2), 1, 33);% unbias each spike waveform
    
        % save spike points and spike waveforms            
        spkTime{1,k} = spike_points;
        spkWave{1,k} = waveforms;
    
        k = k+1;
        clear spike_points waveforms IDX idx;
    end
end