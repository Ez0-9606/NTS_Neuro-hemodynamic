clc;
clear;
close all;

actF = @(x) tanh(x);
bp2latent = @(x) x;
latent2bp = @(x) x;

addpath("..\FUNCTIONS\");

rat = {...
    '19', '21', '22',...
    'B2', 'B3', 'B4', 'B9_a1',...
    'B10_a1', 'B14_a1', 'B15_a1'};

DT = [0.6 0.5 0.75,...
    0.6 0.6 0.6 0.6,...
    0.6 0.6 0.6];

PRE = [-20 -20 -20,...
    -20 -20 -20 -20,...
    -20 -20 -40];

POST = [90 90 90,...
    120 120 120 120,...
    120 100 160];

ntsDELAY = -1*[3 3 3,...
    3 3 3 3,...
    3 5 3];
bpDELAY = [2 3 2,...
    9 9 9 2,...
    0 0 0] - ntsDELAY;

c1 = [42 183 202]/255;
c2 = [0.9290 0.6940 0.1250];
c3 = [254 74 73]/255;%[0.7 0.7 0.7];

load('..\OUTPUTS\hBINDTestKFold.mat');
disp(corND');
disp(corhBIND');
errorND = cell(numel(rat),1);
errorPhaseND = cell(numel(rat),1);
errorRadiND = cell(numel(rat),1);
errorhBIND = cell(numel(rat),1);
errorPhasehBIND = cell(numel(rat),1);
errorRadihBIND = cell(numel(rat),1);
tt = cell(numel(rat,1),2);
for i = 1:numel(rat)
    %% Experimental data loading
    %i = test(ii);
    ratNumb = rat{i};
    data = load(['..\DATA\', ratNumb, '_matlab_save.mat']);
    
    dt = DT(i);
    tr = double(data.trial);
    bDelay = bpDELAY(i);
    nDelay = ntsDELAY(i);
    pre = PRE(i);
    
    % Prepare experimental data
    if i == 6
        len = floor(numel(data.bp_data)/3);
        data.bp_data(len+1:2*len) = [];
        tr = tr-1;
    end
    stimTime = 60;
    if (i == 4) || (i == 5)
        stimTime = 30;
    end
    [tTmp, bp, ~, iso1, iso2] = dataPrepare(data, dt, tr, nDelay, bDelay, pre);
    tTmp = tTmp(:,1);
    [~, sidx] = min(abs(tTmp));
    [~, eidx] = min(abs(tTmp-stimTime));
    tt{i,1} = tTmp;
    tt{i,2} = [sidx, eidx];

    % Get aligned latent space
    [b, Y, ~, ~] = linRegress([ones(numel(bp(:)),1), iso1(:), iso2(:)], bp(:));
    [algn1,algn2] = alignLatentSlide(iso1, iso2, 2);
    [algn1,algn2] = alignLatentScale(algn1, algn2, 2);
    [algn1,algn2] = alignLatentRotate(algn1, algn2, b);

    if contains(['B2', 'B3', 'B4', 'B14_a1', 'B15_a1'], ratNumb)
        algn1 = alignLatentFlip(algn1);
    end

    %% BP dynamics prediction (H-BIND vs Neural Decoding)
    for j = 1:numel(hBIND{i})
        trueY = bp(:,j);
        yAng = bp2PopStatPhase(trueY, sidx);

        %% BP dynamics prediction via Neural Decoding
        testND = latent2bp(algn2(:,j));
        coefND = corCoef(testND, trueY);
        errorND{i,j} = trueY-testND;
        errorPhaseND{i,j} = latent2PopStatPhase(algn1(:,j), algn2(:,j), sidx, eidx)-yAng;
            errorPhaseND{i,j}(errorPhaseND{i,j}>pi) = errorPhaseND{i,j}(errorPhaseND{i,j}>pi)-2*pi;
            errorPhaseND{i,j}(errorPhaseND{i,j}<-pi) = errorPhaseND{i,j}(errorPhaseND{i,j}<-pi)+2*pi;
        errorRadiND{i,j} = sqrt(algn1(:,j).^2+algn2(:,j).^2)-1;

        %% BP dynamics prediction via H-BIND
        hBINDtmp = hBIND{i}{j};
        dt = hBINDtmp.dt;
        u = hBINDtmp.u;
        Wr = hBINDtmp.Wr;
        Wz = hBINDtmp.Wz;
        Wu = hBINDtmp.Wu;
        z = hBINDtmp.z;
        trainInfo = hBINDtmp.trainInfo;
            N = trainInfo.N;
            tau = trainInfo.tau;
            g = 0;
            J = zeros(N);
        [xNew, ~] = hBINDsimulate(N, dt, tau, g, J, actF, tTmp, z, u, Wz, Wu, Wr);
        latent = Wr'*actF(xNew);
        testhBIND = latent2bp(latent(2,:));
        coefhBIND = corCoef(testhBIND,trueY);
        errorhBIND{i,j} = trueY-testhBIND';
        errorPhasehBIND{i,j} = latent2PopStatPhase(latent(1,:)', latent(2,:)',sidx, eidx)-yAng;
            errorPhasehBIND{i,j}(errorPhasehBIND{i,j}>pi) = errorPhasehBIND{i,j}(errorPhasehBIND{i,j}>pi)-2*pi;
            errorPhasehBIND{i,j}(errorPhasehBIND{i,j}<-pi) = errorPhasehBIND{i,j}(errorPhasehBIND{i,j}<-pi)+2*pi;
        errorRadihBIND{i,j} = sqrt(sum((Wr'*actF(xNew)).^2,1)')-1;

        %% Save representative K-fold feedforward simulation output 
        if (i == 7) && (j == 3)
            f1 = figure; hold on;
            t = tTmp;
            buff{1} = algn1;
            buff{2} = algn2;
            algn1 = algn1(:,j);
            algn2 = algn2(:,j);
            plot(t, trueY, 'LineWidth', 3, 'Color','k');
            plot(t, testND, 'LineWidth', 3, 'Color',c2);
            plot(t, testhBIND, 'LineWidth', 3, 'Color','r');
            xlim([-20 80]);
            ylim([-1.2 1.2]);

            save('..\OUTPUTS\hBINDTestKFoldFeedforwardRepresentative1.mat',...
                'actF', 'dt', 'yAng', 't',...
                'xNew', 'algn1', 'algn2', 'Wr');
            savefig(f1, '..\OUTPUTS\RepresentativeFig1.fig');
            algn1 = buff{1};
            algn2 = buff{2};
        end
        
        if (i == 8) && (j == 2)
            f1 = figure; hold on;
            t = tTmp;
            buff{1} = algn1;
            buff{2} = algn2;
            algn1 = algn1(:,j);
            algn2 = algn2(:,j);
            plot(t, trueY, 'LineWidth', 3, 'Color','k');
            plot(t, testND, 'LineWidth', 3, 'Color',c2);
            plot(t, testhBIND, 'LineWidth', 3, 'Color','r');
            xlim([-20 80]);
            ylim([-1.2 1.2]);

            save('..\OUTPUTS\hBINDTestKFoldFeedforwardRepresentative2.mat',...
                'actF', 'dt', 'yAng', 't',...
                'xNew', 'algn1', 'algn2', 'Wr');
            savefig(f1, '..\OUTPUTS\RepresentativeFig2.fig');
            algn1 = buff{1};
            algn2 = buff{2};
        end
    end
    close all;
end

save('..\OUTPUTS\hBINDTestKFoldFeedforwardStat.mat', 'tt', 'errorND', 'errorhBIND', 'errorRadiND', 'errorRadihBIND', 'errorPhaseND', 'errorPhasehBIND');