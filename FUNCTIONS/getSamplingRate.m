function fs = getSamplingRate(varargin)
% Jiho Lee, 2022-12-12
% input: (none)
% 
% output: time array
% 
if numel(varargin) > 1
    error('getSamplingRate error: not proper input');
end
fs = double(typecast(uint32(1186905120), 'single'));

end