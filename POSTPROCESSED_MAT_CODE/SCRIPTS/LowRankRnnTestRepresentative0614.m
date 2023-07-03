clc;
clear;
close all;

addpath('..\..\FUNCTIONS\');

actF = @(x) tanh(x);

params = load('..\..\PARAMETER_DATA\params.mat');

i = 1;
    %% Experimental data loading
    dt = params.DT(i);
    % Normalized neural trajectory
    data = load(['..\DATA\', params.rat{i}, '_matlab_save.mat']);
    [t, bp, ~, ~,iso1, iso2] = dataPrepare(data, dt, params.TRIAL(i), params.ntsDELAY(i), params.bpDELAY(i), params.PRE(i));
    t(:,1) = [];
    bp(:,1) = [];
    iso1(:,1) = [];
    iso2(:,1) = [];
    % set time range
    t = t(:,1);
    % Get aligned latent space
    [b, ~, ~, ~] = linRegress([ones(size(iso1,1),1), mean(iso1,2), mean(iso2,2)], mean(bp,2));
    [~, algn1,algn2] = alignLatentSlide(iso1, iso2, 2);
    [~, algn1,algn2] = alignLatentScale(algn1, algn2, 2);
    [~, algn1,algn2] = alignLatentRotate(algn1, algn2, b);

    tra1 = mean(algn1,2);
    tra2 = mean(algn2,2);
    
    %% Low rank RNN training
        bpNoiseLevel = 0.1;
        
        trainInfo = struct('alpha', 50, 'itr', 20, 'tau', dt, 'N', 500, 'bpNoiseLevel', bpNoiseLevel);

        % set stimulation duration
        stimTime = 60;
        
        % set stimulation input array
        u = zeros(1,size(t,1));
        [~, sidx] = min(abs(t));
        [~, eidx] = min(abs(t-stimTime));
        u(sidx:eidx) = 1;
        
        % Strart training section
        trRP = 5; % training repeat
        ffRP = 3; % feedforward test repeat

        cortmpZ = zeros(trRP,1);
        RNN = cell(trRP,1);
        % Repeat training
        for m = 1:trRP
            fprintf(['------ Training repeat... ', 'm: ', num2str(m), '/', num2str(trRP), '\n']);
            
            % Generate training data
            TrainData = [cos(atan2(tra2, tra1)), sin(atan2(tra2, tra1))];
            % Run training
            RNN{m} = trainLowRankRNN(dt, trainInfo, TrainData, u);
            
            % Feedforward simulation
            corffZ = zeros(ffRP,1);
            corffY = zeros(ffRP,1);
            for n = 1:ffRP
                [xNew, ~] = hBINDsimulate(trainInfo.N, dt, trainInfo.tau, 0, zeros(trainInfo.N), actF, t,...
                    RNN{m}.z, RNN{m}.u, RNN{m}.Wz, RNN{m}.Wu, RNN{m}.Wr);
                test_hBIND = RNN{m}.Wr'*actF(xNew);
                corffZ(n) = mean([corCoef(test_hBIND(1,:), TrainData(:,1)),...
                                  corCoef(test_hBIND(2,:), TrainData(:,2))]);
                
                fprintf(['--------- Feedforward coef.: ',...
                   '(vs. trained data) ', num2str(corffZ(n)),...
                   ';  (n:', num2str(n), '/', num2str(ffRP), ')\n']);
            end

            % Get a mean performance per each repeated training
            cortmpZ(m) = mean(corffZ);
            fprintf(['------ ', 'Training coef.:',...
                '(vs. trained data) ', num2str(cortmpZ(m)),...
                ';  (m:', num2str(m), '/', num2str(trRP), ')\n']);
        end

        % H-BIND results
        idx = find(cortmpZ > 0.9);
        if isempty(idx)
            fprintf('warning::training fail, no coef. vs. trained data is > 0.9\n');
            corRNN = NaN;
            TrainedRNN = [];
        else
            [~, srtIdx] = sort(cortmpZ(idx));
            corRNN = cortmpZ(idx(srtIdx(end)));
            TrainedRNN = RNN{idx(srtIdx(end))};
        end

    fprintf('\n\n\n'); 

save("..\OUTPUTS\LowRankRNNTestRepresentative.mat", "trainInfo", "TrainedRNN", "corRNN", "TrainData");