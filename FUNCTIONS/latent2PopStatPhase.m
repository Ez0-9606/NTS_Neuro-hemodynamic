function ang = latent2PopStatPhase(latentX, latentY, SIDX, EIDX)
    ang = atan2(latentY, latentX);
    
    idx = [];
    sidx = SIDX; % screening range from sidx ~ end, + is to allow non-strict monitoring
    eidx = EIDX; % screening range from sidx ~ end, - is to allow non-strict monitoring

%     for i = 1:numel(ang)-sidx
%         if ang(sidx+i) < (ang(sidx+i-1) - pi) 
%             idx = [idx, sidx+i-1];
%         end
%     end

    for i = 1:sidx
        if ang(i) > 0.25*pi
            ang(i) = ang(i) - 2*pi;
        end
    end

    for i = eidx:numel(ang)
        if ang(i) < -0.5*pi
            ang(i) = ang(i) + 2*pi;
        end
    end

%     switch numel(idx)
%         case 0
% %             figure;
% %             plot(ang);
% %             xline(sidx);
% %             error('error::latent2PopStatPhase, no circulation.\n');
%         case 1
%             ang(1:idx) = ang(1:idx) - 2*pi;
%         case 2
%             ang(1:idx(1)) = ang(1:idx(1)) - 2*pi;
%             ang(idx(2):end) = ang(idx(2):end) + 2*pi;
%         otherwise
%             error('error::latent2PopStatPhase, unexpected waveform.\n');
%     end
end