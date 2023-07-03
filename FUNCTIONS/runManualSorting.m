function [status, cmdout] = runManualSorting(scriptPath, ratPath)

[status, cmdout] = system([scriptPath, 'manualsorting.bat ', '"', ratPath, '" ', '"NTS_preprocessed.kwik"'], '-echo');

end