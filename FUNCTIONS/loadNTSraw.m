function [ntsData, stimData] = loadNTSraw(varargin)

if numel(varargin) < 1
    error('loadNTSdata error: no raw directory');
end

rawDir = varargin{1};

switch numel(varargin)
    case 1
        ntsCh = 16;
        stimCh = 4;
    case 2
        ntsCh = varargin{2};
        stimCh = 4;
    case 3
        ntsCh = varargin{2};
        stimCh = varargin{3};
    otherwise
        error('loadNTSdata error: too many input')
end

% NTS recording data
ntsAddrs = findNTSfile(rawDir);
ntsData = loadNTSdata(ntsAddrs.addrs, ntsCh);

% Stimulation data
stimAddrs = findStimfile(rawDir);
stimData = loadNTSdata(stimAddrs.addrs, stimCh);

end