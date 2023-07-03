function t = getTimeAxis(varargin)
% Jiho Lee, 2022-12-12
% input:
% (1) Length of time array (L) (0 ~ L-1)
% (2) Sampling rate (fs) (selective): 24414.0625 Hz (default)
% 
% output: time array (0:1/fs:L-1)
% 
if (numel(varargin) > 2) || (numel(varargin) == 0)
    error('getTimeAxis error: not proper input');
end
L = varargin{1};
switch numel(varargin)
    case 1
        fs = getSamplingRate();
        t =(0:L-1)/fs; 
    case 2
        fs = varargin{2};
        t =(0:L-1)/fs; 
end

end