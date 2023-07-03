function [sx, aligned1, aligned2] = alignLatentScale(latent1, latent2, axis)
    mx(:,1) = mean(latent1,axis);
    mx(:,2) = mean(latent2,axis);

    sx(1) = 0.5*(max(mx(:,1))-min(mx(:,1)));
    sx(2) = 0.5*(max(mx(:,2))-min(mx(:,2)));

    aligned1 = latent1/sx(1);
    aligned2 = latent2/sx(2);
end



