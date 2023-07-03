clc;
clear;
close all;

actF = @(x) tanh(x);
bp2latent = @(x) x;
latent2bp = @(x) x;

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
params = load('..\..\PARAMETER_DATA\params.mat');
hBIND = cell(numel(rat),1);
testY = cell(numel(rat),1);
corND = zeros(numel(rat),1);
corhBIND = zeros(numel(rat),1);
for i = 1:numel(rat)
    %% Experimental data loading
    ratNumb = rat{i};
    data = load(['..\DATA\', ratNumb, '_matlab_save.mat']);
    
    dt = DT(i);
    tr = params.TRIAL(i);
    bDelay = bpDELAY(i);
    nDelay = ntsDELAY(i);
    pre = PRE(i);
    
    % Prepare experimental data
    if i == 6
        len = floor(numel(data.bp_data)/3);
        data.bp_data(len+1:2*len) = [];
        tr = tr-1;
    end

    %%
    [tt, bp, ~, ~, iso1, iso2] = dataPrepare(data, dt, tr, nDelay, bDelay, pre);
    % set time range
    tt = tt(:,1);
    % Get aligned latent space
    [b, Y, ~, ~] = linRegress([ones(numel(bp(:)),1), iso1(:), iso2(:)], bp(:));
    [~, algn1,algn2] = alignLatentSlide(iso1, iso2, 2);
    [~, algn1,algn2] = alignLatentScale(algn1, algn2, 2);
    [~, algn1,algn2] = alignLatentRotate(algn1, algn2, b);

    if contains(['B2', 'B3', 'B4', 'B14_a1', 'B15_a1'], ratNumb)
        algn1 = alignLatentFlip(algn1);
    end

    %% BP dynamics prediction (H-BIND vs Neural Decoding)
    hBINDkFold = cell(tr,1);
    testYkFold = cell(tr,1);
    corNDtmp = zeros(1,tr);
    corhBINDtmp = zeros(1,tr);

    fprintf([num2str(i), ':', ratNumb, '\n']);
    for j = 1:tr
        % Test BP data
        testYtmp = bp(:,j);
        testYkFold{j} = testYtmp;
        
        %% BP dynamics prediction via Neural Decoding
        % ND results
        testX_ND = latent2bp(algn2(:,j));
        corNDtmp(j) = corCoef(testX_ND, testYtmp);
        %errNDtmp(j) = testYtmp - testX_ND;
        
        clearvars testX_ND;

        %% BP dynamics prediction via H-BIND
        bpNoiseLevel = 0.04;
        trainInfo = struct('alpha', 50, 'itr', 20, 'tau', dt, 'N', 500, 'bpNoiseLevel', bpNoiseLevel);
        
        % set stimulation 
        stimTime = 60;
        if (i == 4) || (i == 5)
            stimTime = 30;
        end
        u = zeros(1,size(bp,1));
        [~, sidx] = min(abs(tt));
        [~, eidx] = min(abs(tt-stimTime));
        u(sidx:eidx) = 1;
        
        % Strart training section
        fprintf(['--- k-Fold: ', num2str(j), '/', num2str(tr), '\n']);


        trRP = 5;
        ffRP = 3;

        cortmpZ = zeros(trRP,1);
        cortmpY = zeros(trRP,1);
        hBINDtmp = cell(trRP,1);
        % Repeat training
        for m = 1:trRP
            %fprintf(['------ Training repeat... ', 'm: ', num2str(m), '/', num2str(trRP), '\n']);
            
            % Generate training data
            [~, TrainData] = trainPrepare(data.bp_data, dt, tr, nDelay, bDelay, pre, trainInfo.bpNoiseLevel);
            TrainData(:,j) = [];
            
            % Run training
            hBINDtmp{m} = trainhBIND(dt, trainInfo, TrainData, u);
            
            % Feedforward simulation
            corffZ = zeros(ffRP,1);
            corffY = zeros(ffRP,1);
            for n = 1:ffRP
                [xNew, ~] = hBINDsimulate(trainInfo.N, dt, trainInfo.tau, 0, zeros(trainInfo.N), actF, tt,...
                    hBINDtmp{m}.z, hBINDtmp{m}.u, hBINDtmp{m}.Wz, hBINDtmp{m}.Wu, hBINDtmp{m}.Wr);
                test_hBIND = hBINDtmp{m}.Wr'*actF(xNew);
                corffZ(n) = mean([corCoef(test_hBIND(1,:), hBINDtmp{m}.z(1,:)),...
                                  corCoef(test_hBIND(2,:), hBINDtmp{m}.z(2,:))]);
                corffY(n) = corCoef(latent2bp(test_hBIND(2,:)), testYtmp);
                
                %fprintf(['--------- Feedforward coef.: ',...
                %    '(vs. trained data) ', num2str(corffZ(n)),...
                %    ';   (vs. test data) ', num2str(corffY(n)), '  (n:', num2str(n), '/', num2str(ffRP), ')\n']);
            end

            % Get a mean performance per each repeated training
            cortmpZ(m) = mean(corffZ);
            cortmpY(m) = mean(corffY);
            fprintf(['------ ', 'Training coef.:',...
                '(vs. trained data) ', num2str(cortmpZ(m)),...
                ';   (vs. test data) ', num2str(cortmpY(m)), '  (m:', num2str(m), '/', num2str(trRP), ')\n']);            
            clearvars corffZ corffY TrainData;
        end

        % H-BIND results
        idx = find(cortmpZ > 0.9);
        if isempty(idx)
            fprintf('warning::training fail, no coef. vs. trained data is > 0.9\n');
            corhBINDtmp(j) = NaN;
            hBINDkFold{j} = [];
        else
            [~, srtIdx] = sort(cortmpY(idx));
            if cortmpY(idx(srtIdx(end))) < 0.8
                fprintf('warning::training fail, no coef. vs. test data is > 0.8\n');
                corhBINDtmp(j) = NaN;
                hBINDkFold{j} = [];
            else
                corhBINDtmp(j) = cortmpY(idx(srtIdx(end)));
                hBINDkFold{j} = hBINDtmp{idx(srtIdx(end))};
            end
        end
        clearvars testYtmp u cortmpZ cortmpY hBINDtmp  idx;

        fprintf(['--- k-Fold coef.: ', num2str(corhBINDtmp(j)), '  (j:', num2str(j), '/', num2str(tr), ')\n']);
        fprintf(['--- corresponding neural decoding coef.: ', num2str(corNDtmp(j)), '\n']);
        fprintf('---\n---\n');
    end

    % Save k-Fold cross validation for each rat
    hBIND{i} = hBINDkFold;
    testY{i} = testYkFold;
    corND(i) = mean(corNDtmp);
    corhBIND(i) = mean(corhBINDtmp);

    clearvars tt hBINDkFold testYkFold corNDtmp corhBINDtmp

    fprintf('\n\n\n'); 
end

figure; hold on;
for i = 1:numel(rat)
    line([1,2], [corND(i), corhBIND(i)], 'LineWidth', 3, 'Color', 'k');
end
scatter(ones(1,numel(rat)), corND);
scatter(2*ones(1,numel(rat)), corhBIND);

% save("hBINDTestKFold.mat", "hBIND", "corND", "corhBIND", "testY");