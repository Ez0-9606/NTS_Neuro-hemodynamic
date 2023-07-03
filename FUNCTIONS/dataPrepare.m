function [t, bp, hr, fr, iso1, iso2] = dataPrepare(postMatData, dt, tr, ntsDelay, bpDelay, pre)
    
    fr_rate = postMatData.firing_rate;
    fr_isomap = postMatData.fr_isomap;
    
    len = floor(size(fr_rate,1)/tr);

    %%
    if isfield(postMatData, 'bp_data')
        bp = postMatData.bp_data;
        bp = dataFold(bp, len, tr);
    else
        bp = [];
    end

    if isfield(postMatData, 'hr_data')
        hr = postMatData.hr_data;
        hr = dataFold(hr, len, tr);
    else
        hr = [];
    end

    %%
    fr = cell(size(fr_rate,2),1);
    for i = 1:size(fr_rate,2)
        fr{i} = dataFold(fr_rate(:,i), len, tr);
    end
    %%
    fr_isoTmp1 = fr_isomap(:,1);
    fr_isoTmp2 = fr_isomap(:,2);
    iso1 = dataFold(fr_isoTmp1, len, tr);
    iso2 = dataFold(fr_isoTmp2, len, tr);
    
    %%
    t = repmat(pre + (0:len-1)'*dt, 1, tr);
    if ntsDelay > 0
        iso1 = iso1(ntsDelay+1:end,:);
        iso2 = iso2(ntsDelay+1:end,:);
        for i = 1:numel(fr)
            tmp = fr{i};
            fr{i} = tmp(ntsDelay+1:end,:);
        end
        t = t(1:end-ntsDelay,:);
        if isfield(postMatData, 'bp_data')
            bp = bp(ntsDelay+1:end,:);
        end
        if isfield(postMatData, 'hr_data')
            hr = hr(ntsDelay+1:end,:);
        end

    elseif ntsDelay < 0
        iso1 = iso1(1:end+ntsDelay,:);
        iso2 = iso2(1:end+ntsDelay,:);
        for i = 1:numel(fr)
            tmp = fr{i};
            fr{i} = tmp(1:end+ntsDelay,:);
        end
        t = t(-ntsDelay+1:end,:);
        if isfield(postMatData, 'bp_data')
            bp = bp(1:end+ntsDelay,:);
        end
        if isfield(postMatData, 'hr_data')
            hr = hr(1:end+ntsDelay,:);
        end
    end

    if bpDelay > 0
        iso1 = iso1(1:end-bpDelay,:);
        iso2 = iso2(1:end-bpDelay,:);
        for i = 1:numel(fr)
            tmp = fr{i};
            fr{i} = tmp(1:end-bpDelay,:);
        end
        t = t(1:end-bpDelay,:);
        if isfield(postMatData, 'bp_data')
            bp = bp(bpDelay+1:end,:);
        end
        if isfield(postMatData, 'hr_data')
            hr = hr(bpDelay+1:end,:);
        end
    elseif bpDelay <0
        iso1 = iso1(-bpDelay+1:end,:);
        iso2 = iso2(-bpDelay+1:end,:);
        for i = 1:numel(fr)
            tmp = fr{i};
            fr{i} = tmp(-bpDelay+1:end,:);
        end
        t = t(-bpDelay+1:end,:);
        if isfield(postMatData, 'bp_data')
            bp = bp(1:end+bpDelay,:);
        end
        if isfield(postMatData, 'hr_data')
            hr = hr(1:end+bpDelay,:);
        end
    end
end