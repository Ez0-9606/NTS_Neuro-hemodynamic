function fov = exp_var_tuning(tuning, fr)
% % %     % Method # 1: Using squared sum
% % %     ss_tot = sum((fr - repmat(mean(fr, 1), size(fr,1), 1)).^2, 1)';
% % %     ss_res = sum((fr - tuning).^2, 1)';
% % %     fov = 1 - ss_res./ss_tot;  

    % Method # 2: Direct calculation of explained variance
    var_tot = var(fr,0 , 1);
    var_exp = var(tuning, 0, 1);
    fov = var_exp./var_tot;

end