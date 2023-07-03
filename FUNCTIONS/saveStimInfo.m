function stimInfo = saveStimInfo()

tmp = input('Input stim on time values:', 's');
stimInfo.stimOnTime = str2double(tmp);

tmp = input('Input stim duration values:', 's');
stimInfo.stimDur = str2double(tmp);

tmp = input('Input stim frequency values:', 's');
stimInfo.stimFreq = str2double(tmp);

end