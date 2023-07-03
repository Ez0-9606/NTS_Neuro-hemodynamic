function functionCallTest(varargin)

% runSpikeSorting()
varargin{1}();
if numel(varargin) > 1
    varargin{2}();
end

end