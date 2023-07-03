clear;
close all;
clc;

addpath('..\..\FUNCTIONS\');

params = load('..\..\PARAMETER_DATA\params.mat');

%% initialization
N = numel(params.rat);

t = cell(N,1); % time courses
bp = cell(N,1); % recorded blood pressure data
hr = cell(N,1);
hrv = cell(N,2);

fr = cell(N,1);

minBP = zeros(N,1);
maxBP = zeros(N,1);
minHR = zeros(N,1);
maxHR = zeros(N,1);
minHRVlf = zeros(N,1);
maxHRVlf = zeros(N,1);
minHRVhf = zeros(N,1);
maxHRVhf = zeros(N,1);

iso1 = cell(N,1);% principal component 1 of Isomap
algn1 = cell(N,1); % tologically aligned principal component 1 of Isomap
raln1 = cell(N,1); % functionally aligned principal component 1 of Isomap
iso2 = cell(N,1); % principal component 2 of Isomap
algn2 = cell(N,1); % tologically aligned principal component 2 of Isomap
raln2 = cell(N,1); % functionally aligned principal component 2 of Isomap

dx = zeros(N,2);
sx = zeros(N,2);

regCoef = zeros(N,3); % linear regression coefficient for BP = a1 + a2*iso1 + a3*iso2

decVec = zeros(N,2); % unit vector of decoding space

YSingle = cell(N,1); % representative neural decoding through rat-specific linear regression
YIndv = cell(N,1); % representative BP prediction by single-neurons
YComm = cell(N,1); % representative neural decoding through across-rat linear regression

accuSing = cell(N,1);
accuSingHR = cell(N,1);
accuSingHRVlfpw = cell(N,1);
accuSingHRVhfpw = cell(N,1);
accuIndv = zeros(1, N);
accuIndvHR = zeros(1, N);
accuIndvHRVlfpw = zeros(1, N);
accuIndvHRVhfpw = zeros(1, N);
accuComm = zeros(1, N);

sqrErrSing = cell(N,1);
sqrErrIndv = cell(N,1);
sqrErrComm = cell(N,1);

unalgnPh = cell(N,1); % 'unaligned' population state phase
algnPh = cell(N,1); % 'aligned' population state phase

sidx = zeros(1, N);
eidx = zeros(1, N);

unalgnPre = []; % unaligned population state phase during pre-stimulation
unalgnStim = []; % unaligned population state phase during stimulation
unalgnPost = []; % unaligned population state phase during post-stimulation
algnPre = []; % aligned population state phase during pre-stimulation
algnStim = []; % aligned population state phase during stimulation
algnPost = []; % aligned population state phase during post-stimulation

unalgnRotPh = cell(1,N); % quantification of unaligned rotation
algnRotPh = cell(1,N); % quantification of aligned rotation

%% calculation (data load, alignment)
for i = 1:N
    % parameter setting
    ratNumb = params.rat{i};
    
    %% load data
    data = load(['..\DATA\', ratNumb, '_matlab_save.mat']);
    dataCar = load(['..\..\PROCESSED\BP_DATA\bp_data_', ratNumb, '_v7.mat']);
    syncBase = (0:size(dataCar.BP,1)-1)*0.01;
    syncIdx = getNormIdx(syncBase, data.t + (params.TrialIDX{i}(1)-1)*(params.POST(i)-params.PRE(i)));
    data.bp_data = dataCar.BP(syncIdx,2);
    data.hr_data = dataCar.HR(syncIdx,2);

    %% HRV calculation
    hrD = reshape(dataCar.HR(:,1), 100*(params.POST(i)-params.PRE(i)), []);
    hrD = hrD(:, params.TrialIDX{i});
    if i == 1
        hrD(:,1) = [];
    end
    df = 30;
    dr = 1/df;
    rangeHRV = (-20:dr:80);
    x = (0:size(hrD,1)-1)*0.01 + params.PRE(i);
    idx = getNormIdx(x, rangeHRV);
    hrD = hrD(idx+(params.ntsDELAY(i)+params.bpDELAY(i))*params.DT(i)*100,:);
    for j = 1:size(hrD,2)
        hrD(:,j) = minmaxNorm(hrD(:,j), [-1 1]);
    end   
    hrD = minmaxNorm(mean(hrD,2), [-1 1]);
    hrvD = heartRateVar(rangeHRV, hrD, df, 10);
    [~, lim] = min(abs(hrvD.fHRV - 0.15));
    hrvLFpwTmp = sum(hrvD.pwHRV(1:lim,:), 1);
    hrvHFpwTmp = sum(hrvD.pwHRV(lim+1:end,:), 1);

    %%
    [tTmp, bpTmp, hrTmp, frTmp, iso1Tmp, iso2Tmp] = dataPrepare(data,...
        params.DT(i), params.TRIAL(i), params.ntsDELAY(i), params.bpDELAY(i), params.PRE(i));
    %%%% !!! exclude one trial from rat#1 !!! %%%%
    if i == 1
        tTmp(:,1) = [];
        bpTmp(:,1) = [];
        hrTmp(:,1) = [];
        %hrvTmp(:,:,1) = [];
        for j = 1:numel(frTmp)
            frTmp{j}(:,1) = [];
        end
        iso1Tmp(:,1) = [];
        iso2Tmp(:,1) = [];
    end

    %% Hemodynamic data processing
    % normalization of data sampling 
    range = (-20:0.5:80);
    tL = numel(range);
    tr = size(bpTmp,2);
    idx = getNormIdx(tTmp(:,1), range);
    tTmp = tTmp(idx,1);
    bpTmp = bpTmp(idx,:);
    hrTmp = hrTmp(idx,:);
    %hrvTmp = hrvTmp(:,idx,:);
    for j = 1:numel(frTmp)
        frTmp{j} = frTmp{j}(idx,:);
    end
    iso1Tmp = iso1Tmp(idx,:);
    iso2Tmp = iso2Tmp(idx,:);

    %%
    %[~, fLimIdx1] = min(abs(data.hrv_data.fHRV - flimit1));
    %[~, fLimIdx2] = min(abs(data.hrv_data.fHRV - flimit2));
    %[~, fLimIdx3] = min(abs(data.hrv_data.fHRV - flimit3));
    %hrvLFpwTmp = reshape(sum(hrvTmp(2:fLimIdx2,:,:),1),tL, tr);
    %hrvHFpwTmp = reshape(sum(hrvTmp(fLimIdx2+1:end,:,:),1),tL, tr);

    [~, bpSidx] = min(abs(range));
    minBPTmp = mean(min(bpTmp,[],1) - mean(bpTmp(1:bpSidx,:),1));
    maxBPTmp = mean(max(bpTmp,[],1)- mean(bpTmp(1:bpSidx,:),1));
    minHRTmp = mean(min(hrTmp,[],1) - mean(hrTmp(1:bpSidx,:),1));
    maxHRTmp = mean(max(hrTmp,[],1) - mean(hrTmp(1:bpSidx,:),1));
    

    [~, hrvsidx] = min(abs(rangeHRV));
    minHRVlfTmp = mean(min(hrvLFpwTmp) - mean(hrvLFpwTmp(1:hrvsidx)));
    maxHRVlfTmp = mean(max(hrvLFpwTmp) - mean(hrvLFpwTmp(1:hrvsidx)));
    minHRVhfTmp = mean(min(hrvHFpwTmp) - mean(hrvHFpwTmp(1:hrvsidx)));
    maxHRVhfTmp = mean(max(hrvHFpwTmp) - mean(hrvHFpwTmp(1:hrvsidx)));

    for j = 1:size(bpTmp,2)
        bpTmp(:,j) = minmaxNorm(bpTmp(:,j), [-1 1]);
        hrTmp(:,j) = minmaxNorm(hrTmp(:,j), [-1 1]);
    end
    hrvLFpwTmp = minmaxNorm(hrvLFpwTmp, [-1 1])';
    hrvHFpwTmp = minmaxNorm(hrvHFpwTmp, [-1 1])';

    %% Hemodynamic prediction by single-neurons
    frR = zeros(numel(frTmp), tL);
    YSingleTmp = zeros(numel(frTmp), tL);
    sqrErrSingTmp = zeros(numel(frTmp), tL);
    corSingTmp = zeros(numel(frTmp),1);
    corSingTmpHR = zeros(numel(frTmp),1);
    corSingTmpHRVlfpw = zeros(numel(frTmp),1);
    corSingTmpHRVhfpw = zeros(numel(frTmp),1);
    
    hrvidx = getNormIdx(rangeHRV, range);
    for j = 1:numel(frTmp)
        % calculate firing rate
        frR(j,:) = mean(frTmp{j},2);
        % BP vs. Single-neurons
        [~, YSingleTmp(j,:), corSingTmp(j), ~] = linRegress([ones(size(frR(j,:)')), frR(j,:)'], mean(bpTmp,2));
        sqrErrSingTmp(j,:) = (YSingleTmp(j,:)-mean(bpTmp,2)').^2;
        % HR vs. Single-neurons
        [~, ~, corSingTmpHR(j), ~] = linRegress([ones(size(frR(j,:)')), frR(j,:)'], mean(hrTmp,2));
        
        % HRV(Low-freq power) vs. Single-neurons
        [~, ~, corSingTmpHRVlfpw(j), ~] = linRegress([ones(size(frR(j,:)')), frR(j,:)'], hrvLFpwTmp(hrvidx));
        % HRV(High-freq power) vs. Single-neurons
        [~, ~, corSingTmpHRVhfpw(j), ~] = linRegress([ones(size(frR(j,:)')), frR(j,:)'], hrvHFpwTmp(hrvidx));
    end

    %% rat-specific latent space decoding (linear regression between latent space and bp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% method 1: not averaging across trials %%%
    %[b, YTmp, corTmp, ~] = linRegress([ones(numel(bpTmp(:)),1), iso1Tmp(:), iso2Tmp(:)], bpTmp(:));
    %YTmp = reshape(YTmp, size(bpTmp,1), size(bpTmp,2));
    %sqrrErrIndvTmp = (YTmp(:) -  bpTmp(:)).^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% method 2: averaging across trials %%%%%%%
    [regCoefTmp, YIndvTmp, corIndvTmp, ~] = linRegress([ones(size(bpTmp,1),1), mean(iso1Tmp,2), mean(iso2Tmp,2)], mean(bpTmp,2));
    sqrErrIndvTmp = (YIndvTmp -  mean(bpTmp,2)).^2;
    % HR vs. decoding
    [~, ~, corIndvTmpHR, ~] = linRegress([ones(tL,1), mean(iso1Tmp,2), mean(iso2Tmp,2)], mean(hrTmp,2));
    % HRV(Low-freq power) vs. decoding
    [~, ~, corIndvTmpHRVlfpw, ~] = linRegress([ones(tL,1), mean(iso1Tmp,2), mean(iso2Tmp,2)], hrvLFpwTmp(hrvidx));
    % HRV(High-freq power) vs. decoding
    [~, ~, corIndvTmpHRVhfpw, ~] = linRegress([ones(tL,1), mean(iso1Tmp,2), mean(iso2Tmp,2)], hrvHFpwTmp(hrvidx));

    % decoding space unit vector
    decVec(i,:) = regCoefTmp(2:3)/norm(regCoefTmp(2:3),2);

    %% across-rat decoding
    % latent space alignment
    [dxTmp, algn1Tmp,algn2Tmp] = alignLatentSlide(iso1Tmp, iso2Tmp, 2);
    [sxTmp, algn1Tmp,algn2Tmp] = alignLatentScale(algn1Tmp, algn2Tmp, 2);
    [angTmp, raln1Tmp,raln2Tmp] = alignLatentRotate(algn1Tmp, algn2Tmp, regCoefTmp);
    if contains(['B2', 'B3', 'B4', 'B14_a1', 'B15_a1'], ratNumb)
        raln1Tmp = alignLatentFlip(raln1Tmp);
    end

    % calculate population state phases (unaligned vs. aligned)
    unalgnPhTmp = atan2(algn2Tmp, algn1Tmp);
    algnPhTmp = atan2(raln2Tmp, raln1Tmp);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% method 1: not averaging across trials %%%
    %YTmp2 = sin(algnPhTmp(:));
    %corCommTmp = corrcoef(YTmp2, bpTmp(:));
    %YTmp2 = reshape(YTmp2, size(bpTmp,1), size(bpTmp,2));
    %sqrErrCommTmp = (YTmp2(:) - bpTmp(:)).^2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% method 2: averaging across trials %%%%%%%
    YCommTmp = mean(sin(algnPhTmp),2);
    corCommTmp = corrcoef(YCommTmp, mean(bpTmp,2));
    corCommTmp = corCommTmp(1,2);
    disp(corCommTmp^2);
    sqrErrCommTmp = (YCommTmp - mean(bpTmp,2)).^2;

    % find idx of stim-on & -off
    [~, sidx(i)] = min(abs(tTmp));
    if (i == 4) || (i == 5)
        [~, eidx(i)] = min(abs(tTmp-30));
    else
        [~, eidx(i)] = min(abs(tTmp-60));
    end

    % quantify the rotation via the population state phases
    unalgnRotPhTmp = calSingleRotPhase(unalgnPhTmp, sidx(i));
    algnRotPhTmp = calSingleRotPhase(algnPhTmp, sidx(i));

    %% data save
    t{i} = tTmp;
    bp{i} = bpTmp;
    hr{i} = hrTmp;
    hrv{i,1} = hrvLFpwTmp(hrvidx);
    hrv{i,2} = hrvHFpwTmp(hrvidx);
    %hrv{i,3} = hrvD.tHRV;

    minBP(i) = minBPTmp;
    maxBP(i) = maxBPTmp;
    minHR(i) = minHRTmp;
    maxHR(i) = maxHRTmp;
    minHRVlf(i) = minHRVlfTmp;
    maxHRVlf(i) = maxHRVlfTmp;
    minHRVhf(i) = minHRVhfTmp;
    maxHRVhf(i) = maxHRVhfTmp;

    fr{i} = frR;

    iso1{i} = iso1Tmp;
    algn1{i} = algn1Tmp;
    raln1{i} = raln1Tmp;

    dx(i,:) = dxTmp;
    sx(i,:) = sxTmp;
    
    iso2{i} = iso2Tmp;
    algn2{i} = algn2Tmp;
    raln2{i} = raln2Tmp;

    decVec;
    
    YSingle{i} = YSingleTmp;% representative neural decoding through rat-specific linear regression
    YIndv{i} = YIndvTmp;% representative BP prediction by single-neurons
    YComm{i} = YCommTmp; % representative neural decoding through across-rat linear regression

    accuSing{i} = corSingTmp.^2; % single-neuronal prediction
    accuSingHR{i} = corSingTmpHR.^2;% HR
    accuSingHRVlfpw{i} = corSingTmpHRVlfpw.^2;% HRV low-freq. power
    accuSingHRVhfpw{i} = corSingTmpHRVhfpw.^2;% HRV high-freq. power

    accuIndv(i) = corIndvTmp^2; % rat-specific prediction
    accuIndvHR(i) = corIndvTmpHR^2;
    accuIndvHRVlfpw(i) = corIndvTmpHRVlfpw^2;
    accuIndvHRVhfpw(i) = corIndvTmpHRVhfpw^2;
    
    accuComm(i) = corCommTmp^2; % across-rat prediction

    sqrErrSing{i} = sqrErrSingTmp;
    sqrErrIndv{i} = sqrErrIndvTmp;
    sqrErrComm{i} = sqrErrCommTmp;

    unalgnPh{i} = unalgnPhTmp;
    algnPh{i} = algnPhTmp;

    unalgnPre = [unalgnPre;reshape(unalgnPhTmp(1:sidx,:),[],1)];
    unalgnStim = [unalgnStim;reshape(unalgnPhTmp(sidx:eidx,:),[],1)];
    unalgnPost = [unalgnPost;reshape(unalgnPhTmp(eidx:end,:),[],1)];
    algnPre = [algnPre;reshape(algnPhTmp(1:sidx,:),[],1)];
    algnStim = [algnStim;reshape(algnPhTmp(sidx:eidx,:),[],1)];
    algnPost = [algnPost;reshape(algnPhTmp(eidx:end,:),[],1)];

    unalgnRotPh{i} = unalgnRotPhTmp;
    algnRotPh{i} = algnRotPhTmp;
end

%% Analysis - Topological distribution of population state phase  (unaligned vs. aligned)
name = {'Pre', 'Stim', 'Post'};
tBmean = zeros(3,1);
tBstd = zeros(3,1);
tTmean = zeros(3,1);
tTstd = zeros(3,1);
plvB = zeros(3,1);
plvT = zeros(3,1);
for i = 1:3
    src = eval(['unalgn',name{i}]);
    tBmean(i) = phaseMean(src); % mean of unaligned population state phase across rat
    tBstd(i) = phaseStd(src); % standard deviation of unaligned population state phase across rat
    plvB(i) = calPhPLV(src,tBmean(i)); % phase-locking value of unaligned population state phase

    src = eval(['algn',name{i}]);
    tTmean(i) = phaseMean(src);
    tTstd(i) = phaseStd(src);
    plvT(i) = calPhPLV(src,tTmean(i));
end

% Save analysis result
AnalysisDynmDist.tBmean = tBmean;
AnalysisDynmDist.tBstd = tBstd;
AnalysisDynmDist.tTmean = tTmean;
AnalysisDynmDist.tTstd = tTstd;
AnalysisDynmDist.plvB = plvB;
AnalysisDynmDist.plvT = plvT;

%% Analysis - State transition matrix (unaligned vs. aligned)
nBin = 30;
edges = (-pi:2*pi/nBin:pi);
bin = edges(1:end-1);
tranStd = zeros(2, nBin);
dsr = 6;
transition = cell(2,1);

% source: population state phase distribution
src = cell(2,1);
src{1} = unalgnPh;
src{2} = algnPh;
for i = 1:2
    transition{i} = zeros(nBin,nBin);
    
    % count state transition
    for j = 1:N
        transition{i} = transition{i} + calStateTranMtx(downsample(src{i}{j}(:),dsr), edges);
    end
    
    % get the probability of state transition for each previous state
    transition{i} = getProbMtx(transition{i},2); 

    % calculate the standard deviation of state transition probability
    for j = 1:nBin
        xVal = sum(sin(edges(1:end-1)).*transition{i}(j,:));
        yVal = sum(cos(edges(1:end-1)).*transition{i}(j,:));
        tranStd(i,j) = (-log(xVal^2 + yVal^2))^0.5;
    end
end

% Save analysis result
AnalysisStatTrans.nBin = nBin;
AnalysisStatTrans.edges = edges;
AnalysisStatTrans.transition = transition;
AnalysisStatTrans.tranStd = tranStd;

%% Analysis - Functional distribution of population state phase (unaligned vs. aligned)
bpTmp = [];
raln2Tmp = [];
bAng = [];
tAng = [];
for i = 1:N
    bpTmp = [bpTmp;bp{i}(:)];
    raln2Tmp = [raln2Tmp;raln2{i}(:)];
    bAng = [bAng;unalgnPh{i}(:)];
    tAng = [tAng;algnPh{i}(:)];
end

bpDist = cell(2,nBin);
bpMean = zeros(2,nBin);
bpStd = zeros(2,nBin);
for i = 1:2
    if i == 1
        src = bAng;
    else
        src = tAng;
    end
    for j = 1:nBin
        bpDist{i,j} = bpTmp((edges(j)<=src)&(edges(j+1)>src));
        %bpMean(i,j) = mean(bpDist{i,j});
        bpStd(i,j) = std(bpDist{i,j});
    end
end

AnalysisFuncDist.bpStd = bpStd;

%% Save the results
save('..\OUTPUTS\expDataAnalysis.mat',...
    "t", "bp", "hr", "hrv",...
    "minBP","maxBP","minHR","maxHR","minHRVlf","maxHRVlf","minHRVhf","maxHRVhf",...
    "fr",...
    "iso1", "algn1", "raln1",...
    "iso2", "algn2", "raln2",...
    "dx", "sx",...
    "decVec",...
    "YIndv", "YSingle", "YComm",...
    "accuSing", "accuSingHR", "accuSingHRVlfpw", "accuSingHRVhfpw",...
    "accuIndv", "accuIndvHR", "accuIndvHRVlfpw", "accuIndvHRVhfpw",...
    "accuComm",...
    "sqrErrSing", "sqrErrIndv", "sqrErrComm",...
    "unalgnPh", "algnPh",...
    "sidx", "eidx",...
    "unalgnPre", "unalgnStim", "unalgnPost",...
    "algnPre", "algnStim", "algnPost",...
    "unalgnRotPh", "algnRotPh",...
    "AnalysisDynmDist","AnalysisStatTrans", "AnalysisFuncDist");