clear;
close all;
clc;

addpath('..\FUNCTIONS\');

c1 = [42 183 202]/255;
c2 = [0.9290 0.6940 0.1250];
c3 = [254 74 73]/255;%[0.7 0.7 0.7];

load('..\POSTPROCESSED_MAT_CODE\OUTPUTS\hBINDTestKFoldFeedforwardRepresentative1.mat'); 

[~, sidx] = min(abs(t)); sidx = sidx -1;
stimTime = 60;
[~, eidx] = min(abs(t-stimTime)); eidx = eidx -1;


figure; hold on;
scatter(cos(yAng), sin(yAng), 100, 'k', 'LineWidth', 1);
axis equal;
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
yRst = mean([cos(yAng(sidx-round(3/dt):sidx-1)), sin(yAng(sidx-round(3/dt):sidx-1))],1);
yAdp = mean([cos(yAng(eidx-round(3/dt):eidx-1)), sin(yAng(eidx-round(3/dt):eidx-1))],1);
% yRcv = mean([cos(yAng(end-round(3/dt):end-1)), sin(yAng(end-round(3/dt):end-1))],1)*0.9;
quiver(yRst(1)*0.5,yRst(2)*0.5,yRst(1)*0.4,yRst(2)*0.4,0, 'Color', 'g', 'LineWidth', 3, 'MaxHeadSize',1.5);
quiver(yAdp(1)*0.5,yAdp(2)*0.5,yAdp(1)*0.4,yAdp(2)*0.4,0, 'Color', 'b', 'LineWidth', 3, 'MaxHeadSize',1.5);
scatter(yRst(1), yRst(2), 80, 'g', 'filled');
scatter(yAdp(1), yAdp(2), 80, 'b', 'filled');
% quiver(0,0,yRcv(1),yRcv(2),0, 'Color', 'm', 'LineWidth', 3, 'MaxHeadSize',0.7);
time = [0 20 40 60 80 120];
for i = 1:numel(time)
    [~, tIdx] = min(abs(t-time(i)));
    tIdx = tIdx -1;
    line([1.1*cos(yAng(tIdx)) 1.2*cos(yAng(tIdx))], [1.1*sin(yAng(tIdx)) 1.2*sin(yAng(tIdx))], 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
end

figure; hold on;
plot(cos(2*pi*(1:100)/100), sin(2*pi*(1:100)/100), 'k', 'LineWidth', 1);
quiver(yRst(1)*0.5,yRst(2)*0.5,yRst(1)*0.4,yRst(2)*0.4,0, 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'MaxHeadSize',1.5);
quiver(yAdp(1)*0.5,yAdp(2)*0.5,yAdp(1)*0.4,yAdp(2)*0.4,0, 'Color', [0.5 0.5 0.5], 'LineWidth', 3, 'MaxHeadSize',1.5);
scatter(algn1, algn2, 100, c2, 'LineWidth', 1);
scatter(Wr(:,1)'*actF(xNew), Wr(:,2)'*actF(xNew), 100, 'r', 'LineWidth', 1);
NDAng = latent2PopStatPhase(algn1, algn2, sidx, eidx);
hBINDAng = latent2PopStatPhase(Wr(:,1)'*actF(xNew), Wr(:,2)'*actF(xNew), sidx, eidx);
axis equal;
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
axis off;


f = figure; f.Position(3:4) = [200 150]; hold on;
plot(t, (NDAng-yAng)/pi, 'Color', c2, 'LineWidth', 2);
plot(t, (hBINDAng'-yAng)/pi, 'r', 'LineWidth', 2);
xlim([-20 120]);

f = figure; f.Position(3) = 400; hold on;
plot(t, 1-(algn1.^2+algn2.^2).^0.5, 'Color', c2, 'LineWidth', 3);
plot(t, 1-((Wr(:,1)'*actF(xNew)).^2+(Wr(:,2)'*actF(xNew)).^2).^0.5, 'r', 'LineWidth', 3);

base = mean(yAng(1:sidx));
f = figure; f.Position(3) = 400; hold on;
plot(t, (yAng-base)/pi, 'k', 'LineWidth', 3);
plot(t, (NDAng-base)/pi, 'Color', c2, 'LineWidth', 3);
plot(t, (hBINDAng-base)/pi, 'r', 'LineWidth', 3);
ylim([-0.1 2.1]);
% ylim([-1 1]);

