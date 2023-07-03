function [status, cmdout] = runSpikeSorting(scriptPath, ratPath)

[status, cmdout] = system([scriptPath, 'spikesorting.bat ', '"', ratPath, '"'], '-echo');

end