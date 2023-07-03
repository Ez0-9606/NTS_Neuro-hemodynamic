close all;

%% Analysis - training performance

errFig = figure;
plot(t(2:end), sqrt(sum(ePlus.^2,1))./sqrt(sum(eMinus.^2,1)));
xlabel('Time (s)');
ylabel('Ratio : SE(t+1)/SE(t)');
title('Training error');

% %% Analysis
% figure; hold on;
% 
nFr = actF(x);
[~, mxIdx] = max(nFr, [], 2);
mxIdx = mod(mxIdx - numel(tPre),numel(t));
[~, srtIdx] = sort(mxIdx);
% 
% for i = 1:numel(srtIdx)
%     imagesc(t, i*ones(size(t)), nFr(srtIdx(i),:));
% end
% 
% ylim([0.5 N+0.5]);
% xlim([t(1) t(end)]);
% xline(tStim(1), 'Color', 'k', 'LineWidth',3);
% xline(tStim(end), 'Color', 'k', 'LineWidth',3);
% % colormap viridis(100);


%% Analysis - feedforward simulation of trained network
f1 = figure;f1.Position(1:2) = [70 50];
f2 = figure;f2.Position(2) = 50;

xNew = zeros(size(x)); % with stimulation
xNew2 = zeros(size(x)); % without stimulation

for j = 1:5
    if j == 1
        xNew(:,1) = x(:,end);
        xNew2(:,1) = x(:,end);
    else
        xNew(:,1) = xNew(:,end);
        xNew2(:,1) = xNew2(:,end);
    end
    for i = 1:numel(t)-1
        % update network states
        dx = (dt./tau).*(-xNew(:,i) + g*J*actF(xNew(:,i)) + Wz*Wr'*actF(xNew(:,i)) + 1*Wu*u(i) + randn(N,1)/sqrt(N));
        xNew(:,i+1) = xNew(:,i) + dx;
    
        dx2 = (dt./tau).*(-xNew2(:,i) + g*J*actF(xNew2(:,i)) + Wz*Wr'*actF(xNew2(:,i)) + 0*Wu*u(i) + randn(N,1)/sqrt(N));
        xNew2(:,i+1) = xNew2(:,i) + dx2;
    end
    figure(f1);
    subplot(1,2,1);
    scatter(z(1,:), z(2,:), 50, 'k');hold on;
    scatter(Wr(:,1)'*actF(xNew), Wr(:,2)'*actF(xNew), 50, t, 'filled');
    hold off;
    subplot(1,2,2); hold on;
    scatter(z(1,:), z(2,:), 50, 'k');hold on;
    scatter(Wr(:,1)'*actF(xNew2), Wr(:,2)'*actF(xNew2), 50, t, 'filled');
    hold off;
    
    figure(f2);
    subplot(2,1,1); hold on;
    plot(t+(j-1)*dt*numel(t), z(1,:), 'Color', 'k');
    plot(t+(j-1)*dt*numel(t), Wr(:,1)'*actF(xNew), 'Color', 'r');
    plot(t+(j-1)*dt*numel(t), Wr(:,1)'*actF(xNew2), 'Color', 'b');
    xline((j-1)*dt*numel(t), '-', {['simulation:', num2str(j)]});
    subplot(2,1,2); hold on;
    plot(t+(j-1)*dt*numel(t), z(2,:), 'k');
    plot(t+(j-1)*dt*numel(t), Wr(:,2)'*actF(xNew), 'Color', 'r');
    plot(t+(j-1)*dt*numel(t), Wr(:,2)'*actF(xNew2), 'Color', 'b');
    xline((j-1)*dt*numel(t), '-', {['simulation:', num2str(j)]});
end

figure(f1);
subplot(1,2,1);
xlim([-1.01 1.01]);
ylim([-1.01 1.01]);
xlabel('PC_1');
ylabel('PC_2');
legend({'Original', 'Trained'});
title('Neural trajectory (w Stim)');
subplot(1,2,2);
xlim([-1.01 1.01]);
ylim([-1.01 1.01]);
xlabel('PC_1');
ylabel('PC_2');
legend({'Origianl', 'Trained'});
title('Neural trajectory (w/o Stim)');

figure(f2);
subplot(2,1,1);
ylabel('PC_1');
legend({'Orig(w.Stim)', 'w Stim.', 'w/o Stim.'});
subplot(2,1,2);
ylabel('PC_2');
legend({'Orig(w.Stim)', 'w Stim.', 'w/o Stim.'});

%% Analysis - network activity
netFig = figure; hold on; netFig.Position(1:2) = [1300 50];
maxMtx = repmat(max(actF(xNew),[],2),1,size(xNew,2));
minMtx = repmat(min(actF(xNew),[],2),1,size(xNew,2));
nFr = (actF(xNew) - minMtx)./(maxMtx-minMtx);
for i = 1:numel(srtIdx)
    imagesc(t, i*ones(size(t)), actF(xNew(srtIdx(i),:)));
end
ylim([0.5 N+0.5]);
xlim([t(1) t(end)]);
xline(tStim(1), '-', {'Stim-On'}, 'Color', 'r', 'LineWidth',3);
xline(tStim(end), '-', {'Stim-On'}, 'Color', 'r', 'LineWidth',3);
xlabel('Time (s)');
ylabel('# of neuron in network');

clearvars xNew xNew2