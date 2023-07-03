function [dx, aligned1, aligned2] = alignLatentSlide(latent1, latent2, axis)
    mx(:,1) = mean(latent1,axis);
    mx(:,2) = mean(latent2,axis);

    dx(1) = 0.5*(max(mx(:,1))+min(mx(:,1)));
    dx(2) = 0.5*(max(mx(:,2))+min(mx(:,2)));

    aligned1 = latent1-dx(1);
    aligned2 = latent2-dx(2);
end