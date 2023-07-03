function ang = bp2PopStatPhase(bp,sidx)
    latent = bp;
    if ~isempty(find(abs(latent)>1, 1))
        tmp = find(abs(latent)>1);
        for i = 1:numel(tmp)
            if abs(latent(tmp(i))) < (1+10e-3)
                latent(tmp(i)) = double(sign(latent(tmp(i))));
                continue;
            else
                error('error::bp2PopStatPhase, input bp value is > 1 (or < -1)\n');
            end
        end
    end

    ang = asin(latent);
    
    idxUp = find(abs(ang(sidx:end) - (0.5*pi)) < 0.01); idxUp = idxUp + sidx-1;
    idxDn = find(abs(ang(sidx:end) - (-0.5*pi)) < 0.01); idxDn = idxDn + sidx-1;
    
    if ~(isempty(idxUp) || isempty(idxDn))
        if idxUp < idxDn
            error('error::bp2PopStatPhase, unexpected waveform\n');
        end
    end

    ang(1:idxDn) = -pi - ang(1:idxDn);
    ang(idxUp:end) = pi - ang(idxUp:end);

    idxUp2 = find(abs(ang(sidx:end) - (1.5*pi)) < 0.1); idxUp2 = idxUp2 + sidx-1;
    idxDn2 = find(abs(ang(sidx:end) - (-1.5*pi)) < 0.1); idxDn2 = idxDn2 + sidx-1;
    ang(idxUp2:end) = 3*pi - ang(idxUp2:end);
    ang(1:idxDn2) = -3*pi - ang(1:idxDn2);

end