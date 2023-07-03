function out = calSingleRotPhase(phase, baseIdx)
    out = zeros(size(phase));
    for i = 1:size(phase,2)
        tmp = phase(:,i);
        for j = 2:size(phase,1)
            if (tmp(j) - tmp(j-1))> pi
                tmp(j:end) = tmp(j:end) - 2*pi;
            elseif (tmp(j) - tmp(j-1)) < -pi
                tmp(j:end) = tmp(j:end) + 2*pi;
            end
        end
        tmp = tmp - mean(tmp(1:baseIdx));
        out(:,i) = tmp;
        clearvars tmp;
    end
end