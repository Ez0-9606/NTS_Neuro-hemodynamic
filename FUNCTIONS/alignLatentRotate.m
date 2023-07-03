function [decodingAng, Raligned1, Raligned2] = alignLatentRotate(latent1, latent2, b)
    
    decodingAng = atan2(-b(3), -b(2));
    dtheta = -pi/2 - decodingAng;
    aligned = [cos(dtheta), -sin(dtheta);sin(dtheta), cos(dtheta)]*[latent1(:),latent2(:)]';
    
    Raligned1 = dataFold(aligned(1,:), size(latent1,1), size(latent1,2));
    Raligned2 = dataFold(aligned(2,:), size(latent2,1), size(latent2,2));
end