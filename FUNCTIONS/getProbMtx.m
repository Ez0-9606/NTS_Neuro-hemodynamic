function probMtx = getProbMtx(varargin)

    if numel(varargin) ~= 2
        error('getProbMtx:: Wrong argument input\n');
    end

    Mtx = varargin{1};
    axis = varargin{2};

    sumMtx = sum(Mtx,axis);
    switch axis
        case 1
            probMtx = Mtx./repmat(sumMtx, size(Mtx,1), 1);
        case 2
            probMtx = Mtx./repmat(sumMtx, 1, size(Mtx,2));
    end
end