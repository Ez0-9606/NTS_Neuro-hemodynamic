function  [state, tuning] = get_tuning_curve(x, fr)
    tuning = zeros(max(x) - min(x) +1, size(fr,2));
    count = zeros(max(x) - min(x) + 1, size(fr,2));
    state = (0:max(x) - min(x));

    for i = 1:size(fr,1)
        for j = 1:size(fr,2)
            idx = x(i)-min(x)+1;
            tuning(idx, j) = tuning(idx, j) + fr(i,j);
            count(idx, j) = count(idx, j) + 1;
        end
    end
    
    tuning = tuning./count;
    tuning = interp_missing_in_tuning(tuning);    
end