function ntsData = loadNTSdata(path, ntsCh)

ntsID = fopen(path, 'r');
ntsData = fread(ntsID, 'int16');
ntsData = reshape(ntsData, ntsCh, numel(ntsData)/ntsCh);
fclose(ntsID);

end