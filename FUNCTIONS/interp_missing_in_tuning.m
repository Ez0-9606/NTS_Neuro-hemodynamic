function tuning_intrp = interp_missing_in_tuning(tuning)
    x = (0:size(tuning,1)-1);
    tuning_intrp = zeros(size(tuning));
    for i = 1:size(tuning, 2)
        idx = ~isnan(tuning(:,i));
        tuning_intrp(:,i) = interp1(x(idx), tuning(idx,i), x);
    end
end