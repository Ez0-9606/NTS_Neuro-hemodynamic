function transition = calStateTranMtx(src, edges)

    N = numel(edges)-1;

    transition = zeros(N,N);

    tmp = zeros(size(src)); % downsampling the source (dt = 0.5 => dt = 0.5*dsr)
    for i = 1:N
        tmp((edges(i)<=src)&(edges(i+1)>src)) = i; % digitizing the source
    end
    
    % count the number of state transition
    for i = 1:N
        tmp2 = circshift(tmp == i,1,1);
        tmp2(1,:) = false;
        tmp3 = tmp(tmp2); % get the posterior states which have previous state of 'k'
        for l = 1:numel(tmp3)
            transition(i, tmp3(l)) = transition(i, tmp3(l)) + 1; % add a count for the case of transition from previous state 'k' to posterior state 'tmp3(l)'
            % thus, transition{i}(x, y) -> x:previous state, y:posterior state
        end
        clearvars tmp2 tmp3;
    end
    clearvars tmp;
end

