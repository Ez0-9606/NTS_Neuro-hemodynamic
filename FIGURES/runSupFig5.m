clear;
close all;
clc;

addpath('..\FUNCTIONS\');

params = load('..\PARAMETER_DATA\params.mat');

orD = load('..\POSTPROCESSED_MAT_CODE\DATA\19_matlab_save.mat');
orDpost = load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\expDataAnalysis.mat');
%% initialization
tr = params.TRIAL(1);

tL = numel(orDpost.t{1});
bp = orDpost.bp{1};
frO = orDpost.fr{1};
iO1 = orDpost.iso1{1};
iO2 = orDpost.iso2{1};
dVecO = orDpost.decVec(1,:);
accO = orDpost.accuIndv(1);
sidx = orDpost.sidx(1);

figure; hold on;
scatter(iO1(:), iO2(:), 80, bp(:), 'filled');
colormap viridis(100);clim([-1 1]);
quiver(0, 0, dVecO(1), dVecO(2), 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize',0.5);
axis equal;
axis off;
h1O = orD.h1;
% figure; hold on;
% for i = 1:size(frO,1)
%     plot((-20:0.5:80), minmaxNorm(frO(i,:), [-0.2 0.2]) + i, 'k', 'LineWidth', 2);
% end
figure;
[X,Y] = meshgrid((1:7),(1:4));
X = X+(rand(4,7)-0.5);
Y = Y+(rand(4,7)-0.5);
scatter(X(:),Y(:), 200, 'r', 'filled', 'MarkerEdgeColor','k');
axis off;


subD = load('..\POSTPROCESSED_MAT_CODE\DATA\subpopulationAnalysis_matlab_save_50%.mat');
    subD = subD.subpopulationAnalysis;

    dxfig = figure; hold on; axis equal;
    xlim([-4 14]);ylim([-4 14]);xline(0, 'LineWidth', 3', 'Color', 'g'); yline(0, 'LineWidth', 3', 'Color', 'g');
    sxfig = figure; hold on; axis equal;
    xlim([0 30]);ylim([0 30]);plot((0:0.01:30), (0:0.01:30), 'LineWidth', 3, 'Color', 'g');
    axfig = figure; hold on;
    xlim([-20 80]); yticks([-2*pi 0 2*pi]); yticklabels({'-2\pi', '0', '2\pi'});
    axfigProp = figure; hold on;

    fig = figure; hold on;

    algfig = figure; hold on;
    

    axCount = 0;
    for i = 1:numel(subD)
        subpopIdx = subD{i}.subpopIdx+1;
    
        [t, ~, ~, fr, iS1, iS2] = dataPrepare(subD{i},...
            params.DT(1), params.TRIAL(1), params.ntsDELAY(1), params.bpDELAY(1), params.PRE(1));
        range = (-20:0.5:80);
        idx = getNormIdx(t(:,1), range);
        
        cellN = numel(fr);
        frS = zeros(tL, cellN);
        for j = 1:cellN
            frS(:,j) = mean(fr{j}(idx,2:end),2);
        end
        iS1 = iS1(idx,2:end);
        iS2 = iS2(idx,2:end);

        h1S = subD{i}.h1;
        thS = subD{i}.th;
        ngS = subD{i}.avg_neigh_num;

        %%
        [dx, tmp1, tmp2] = alignLatentSlide(iS1, iS2, 2);
        [sx, tmp1, tmp2] = alignLatentScale(tmp1, tmp2, 2);
        ax = atan2(tmp2, tmp1);
        for k = 1:size(ax,2)
            for l = 2:size(ax,1)
                if ax(l,k) > ax(l-1,k) + 0.7*pi
                    ax(l:end,k) = ax(l:end,k) - 2*pi;
                elseif ax(l,k) < ax(l-1,k) - 0.7*pi
                    ax(l:end,k) = ax(l:end,k) + 2*pi;
                end
            end
        end
        ax = mean(ax,2);
        [~, sidx] = min(abs(range));
        ax = ax - mean(ax(1:sidx));
        if ax(end) > 0
            axCount = axCount +1;
        end
        figure(axfig);
        plot(range, ax, 'k');

        figure(dxfig);
        scatter(dx(1), dx(2), 80, 'k', 'LineWidth', 1);

        figure(sxfig);
        scatter(sx(1), sx(2), 80, 'k', 'LineWidth', 1);
    
        %% Alignment
        [~, test1, cTmp, ~] = linRegress([ones(tL,1), mean(iS1,2), mean(iS2,2)], mean(bp,2));
        [rTmp, ~, ~, ~] = linRegress([ones(numel(iS1(:)),1), iS1(:), iS2(:)], bp(:));
        dVecS(i,:) = rTmp(2:3)/norm(rTmp(2:3),2);

        [~, iSalg1, iSalg2] = alignLatentRotate(tmp1, tmp2, rTmp);

        if ax(end) < 0
            iSalg1 = -iSalg1;
        end

        figure(algfig);
        scatter(mean(iSalg1,2), mean(iSalg2,2), 5, mean(bp,2));
        colormap viridis(100);
        
        %% BP prediction performance
        [~, test2, cTmpAlg, ~] = linRegress([ones(tL,1), mean(iSalg1,2), mean(iSalg2,2)], mean(bp,2));
        
        accS(i) = cTmp^2;
        cTmpAlg = corrcoef(mean(iSalg2,2), mean(bp,2));
        accSalg(i) = cTmpAlg(1,2)^2;
        
        if accS(i) > accSalg(i) + 0.2
            figure;
            plot(mean(bp,2)); hold on;
            plot(test1);
            plot(test2);
            print('test\n');
        end
%         if ismember(i, [23 40 46 52]) && r == 2
%             figure;
%             scatter(iS1(:), iS2(:), 80, bp(:), 'filled'); hold on;
%             quiver(0, 0, dVecS(i,1), dVecS(i,2), 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize',10); hold off;
%             colormap viridis(100);clim([-1 1]);
%             axis equal;
%             axis off;
%     
%             figure;
%             c = zeros(28,1);
%             c(subpopIdx) = 1;
%             scatter(X(:),Y(:), 200, c, 'filled', 'MarkerEdgeColor','k');colormap([0 0 0;1 0 0 ]);
%             axis off;
%         end
   
        
%         figure(fig);
%         quiver(0, 0, dVecS(i,1), dVecS(i,2), 0, 'Color', 'k', 'LineWidth', 1, 'MaxHeadSize',10);
%         axis equal;
%         xlim([-1 1]); ylim([-1 1]);
%         axis off;
        %% H1 longest lifetime distribution
        hL = h1S(:,2)-h1S(:,1);
        zhLong(i) = (max(hL) - mean(hL))/std(hL);
    end
figure(axfigProp);
b = bar(1, [axCount; 100 - axCount], 'stacked');

figure(algfig);
ylim([-1.5 1.5])
    xlim([-1.5 1.5])
    axis equal
    plot(cos(2*pi*(0:100)/100), sin(2*pi*(0:100)/100), 'r', 'LineWidth', 0.5, 'LineStyle', '--');

figure;
histogram(accS, (0:0.05:1), 'Normalization','probability'); box off;
ylim([0 0.25]);

figure;
ax = axes();
hold(ax);
boxchart(ones(100,1), accS,...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(100,1), accSalg,...
    'BoxFaceColor', 'r', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
plot([ones(1,100);ones(1,100)+1], [accS; accSalg], 'Color', [0 0 0 0.2]);
xlim([0.5 2.5]);
ylim([0 1]);

return;
%%
figure;
ax = axes();
hold(ax);
boxchart(ones(100,1), zhLong{1},...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(2*ones(100,1), zhLong{2},...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
boxchart(3*ones(100,1), zhLong{3},...
    'BoxFaceColor', 'k', 'BoxFaceAlpha',0.2, 'BoxMedianLineColor','k', 'LineWidth', 1.5, 'MarkerStyle', '+', 'MarkerColor','k', 'MarkerSize', 10);
scatter(1+0.3*(rand(100,1)-0.5), zhLong{1}, 30, 'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(2+0.3*(rand(100,1)-0.5), zhLong{2}, 30, 'MarkerFaceColor','w', 'MarkerEdgeColor','k');
scatter(3+0.3*(rand(100,1)-0.5), zhLong{3}, 30, 'MarkerFaceColor','w', 'MarkerEdgeColor','k');
ylim([0 20]);