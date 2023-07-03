function tuning_n = norm_tuning(tuning)
    min_mtx = repmat(min(tuning,[],1), size(tuning, 1), 1);
    max_mtx = repmat(max(tuning,[],1), size(tuning, 1), 1);
    tuning_n = (tuning - min_mtx)./(max_mtx - min_mtx);
end