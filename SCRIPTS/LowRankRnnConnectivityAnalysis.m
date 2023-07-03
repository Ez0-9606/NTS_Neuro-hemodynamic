clear;
close all;
clc;
load('..\OUTPUTS\LowRankRnnTestData.mat');

actF = @(x) tanh(x);

figure;
scatter(Wz(:,1), N*Wr(:,2),10, [0.9290 0.6940 0.1250]);
hold on;
scatter(Wz(:,1), N*Wr(:,1),10, [0.5, 0.5 0.5]);
xlim([-5 5]);

figure;
scatter(Wz(:,2), N*Wr(:,1),10, [0.5, 0.5 0.5]);
hold on;
scatter(Wz(:,2), N*Wr(:,2), 10, [0.9290 0.6940 0.1250]);
xlim([-5 5]);

c1 = [42 183 202]/255;
c2 = [0.9290 0.6940 0.1250];
c3 = [254 74 73]/255;%[0.7 0.7 0.7];


figure; hold on;
idx = find((abs(Wz(:,2))<0.2)&(abs(Wu)<0.2));
scatter(Wz(idx,1), N*Wr(idx,1),50,c1, 'filled');
idx = find((abs(Wz(:,1)-1)<0.2)&(abs(Wu)<0.2));
scatter(Wz(idx,2), N*Wr(idx,1),50,c2, 'filled');
idx = find((abs(Wz(:,1)-1)<0.2)&(abs(Wz(:,2))<0.2));
scatter(Wu(idx), N*Wr(idx,1),50,c3, 'filled');

figure; hold on;
idx = find((abs(Wz(:,1))<0.2)&(abs(Wu)<0.2));
scatter(Wz(idx,2), N*Wr(idx,2),50,c1, 'filled');
idx = find((abs(Wz(:,2)-1)<0.2)&(abs(Wu)<0.2));
scatter(Wz(idx,1), N*Wr(idx,2),50,c2, 'filled');
idx = find((abs(Wz(:,2)-1)<0.2)&(abs(Wz(:,1))<0.2));
scatter(Wu(idx), N*Wr(idx,2),50,c3, 'filled');

f = figure; hold on; f.Position(3) = 450;
scatter3(Wz(:,2), Wu, Wz(:,1), 10, Wr(:,1))
xlabel('Wz2')
ylabel('Wu')
zlabel('Wz1')
view(45, 45);
scatter3(Wz(:,2)*0-6, Wu, Wz(:,1), 10, Wr(:,1));
scatter3(Wz(:,2), Wu*0+7, Wz(:,1), 10, Wr(:,1));
xlim([-6 3]);
ylim([-4 7]);

f = figure; hold on; f.Position(3) = 450;
scatter3(Wz(:,1), Wu, Wz(:,2), 10, Wr(:,2))
xlabel('Wz1')
ylabel('Wu')
zlabel('Wz2')
view(45, 45);
scatter3(Wz(:,1)*0-6, Wu, Wz(:,2), 10, Wr(:,2));
scatter3(Wz(:,1), Wu*0+7, Wz(:,2), 10, Wr(:,2));
xlim([-6 3]);
ylim([-4 7]);
return;


xx = 2*rand(3,50000)-1;
scatter(xx(1,:), Wr(:,1)'*actF(Wz*xx(1:2,:) + Wu*xx(3,:)));
scatter3(xx(1,:), xx(2,:), xx(3,:), 10, Wr(:,1)'*actF(Wz*xx(1:2,:) + Wu*xx(3,:)));
xlabel('z1')
ylabel('z2')
zlabel('u')
figure;
scatter3(xx(1,:), xx(2,:), xx(3,:), 10, Wr(:,2)'*actF(Wz*xx(1:2,:) + Wu*xx(3,:)));
xlabel('z1')
ylabel('z2')
zlabel('u')
scatter3(Wz(:,1), Wz(:,2), Wu, 10, Wr(:,1));
xlabel('Wz1')
ylabel('Wz2')
zlabel('Wu')
return;

figure;subplot(1,2,1);
scatter(Wr(:,1)'*actF(x), Wr(:,2)'*actF(x), 50, t);
subplot(1,2,2);
[~, score, ~] = pca(actF(x'));
scatter(score(:,1), score(:,2), 50, t);

N = size(x,1);
% for i = 1:50
% idx = randperm(N, round(N*0.01));
% [~, score, ~] = pca(actF(x(idx,:)'));
% scatter(score(:,1), score(:,2), 50, t); hold on;
% xlim([-6 6]);
% ylim([-6 6]);
% end

zz = 2*rand(2, 50000)-1;
for i = 1:size(zz,2)
fun = @(a,b,c) exp(-a.^2)/sqrt(2*pi).*exp(-b.^2)/sqrt(2*pi).*exp(-c.^2)/sqrt(2*pi)...
    .*(6*tanh(0.7*a)).*(4*exp(-b.^2/5)).*cos(1.4*c)...
    .*tanh(zz(1,i).*a + zz(2,i).*b + 0*c);
vOn(i) = integral3(fun, -10, 10, -10, 10, -10, 10);
end

for i = 1:size(zz,2)
fun = @(a,b,c) exp(-a.^2)/sqrt(2*pi).*exp(-b.^2)/sqrt(2*pi).*exp(-c.^2)/sqrt(2*pi)...
    .*(6*tanh(0.7*a)).*(4*exp(-b.^2/5)).*cos(1.4*c)...
    .*tanh(zz(1,i).*a + zz(2,i).*b + 1*c);
vOff(i) = integral3(fun, -10, 10, -10, 10, -10, 10);
end