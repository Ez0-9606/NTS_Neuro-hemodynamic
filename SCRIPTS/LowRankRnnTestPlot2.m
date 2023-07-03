close all

figure;

xx = cos(2*pi*(1:N)/N)+0.5*(rand(1,N));
yy = sin(2*pi*(1:N)/N)+0.5*(rand(1,N));
test = 1-[zeros(100,1), (1:100)'/100, (1:100)'/100];
for i = 100:numel(t)
    scatter(xx(:),yy(:), 60, actF(x(srtIdx,i)), 'filled', 'MarkerEdgeColor', 'k');
%     clim([-1 1]);
%     xlim([-1.5 1.5]);
%     ylim([-1.5 1.5]);
    
    colormap(viridis(100));
%     colormap(test);
%     hold off;
end
return;

zx = cos(2*pi*(1:300)/300);
zy = sin(2*pi*(1:300)/300);
% zx = (-2:)
% [X,Y] = meshgrid(zx,zy);
X = zx; Y = zy;
zNew = [X(:)';Y(:)'];
v = -zNew + Wr'*actF(Wz*zNew + Wu*0);
v = sum(v.^2,1).^0.5;
scatter(X(:),Y(:), 30, v, 'filled'); hold on;
% scatter(zx, zy, 1, 'k');
axis equal;
xlim([-1 1])
ylim([-1 1])
% return;
colormap viridis(100);
% figure;
zx = cos(2*pi*(1:30)/30);
zy = sin(2*pi*(1:30)/30);
X = zx;
Y = zy;
zx = (-1.1:0.1:1.1);
zy = (-1.1:0.1:1.1);
% [X,Y] = meshgrid(zx,zy);

% for i = 1:numel(zx)
%     for j = 1:numel(zy)
%         xTmp = zeros(N,1);
%         for k = 1:10
%             dx = (dt./tau).*(-xTmp + g*J*actF(xTmp) + Wz*[zx(i);zy(j)] + 1*Wu);
%             xTmp = xTmp + dx;
%         end
%         mtx(:,j,i) = g*J*xTmp;
%     end
% end
Vx = -X(:)' + Wr(:,1)'*actF(Wz*[X(:)';Y(:)'] + Wu*1);% + reshape(mtx,N, numel(mtx(:))/N));
Vx = reshape(Vx, size(X,1), size(X,2));
Vy = -Y(:)' + Wr(:,2)'*actF(Wz*[X(:)';Y(:)'] + Wu*1);% + reshape(mtx,N, numel(mtx(:))/N));
Vy = reshape(Vy, size(Y,1), size(Y,2));
% vx = -zx + Wr(:,1)'*actF(Wz*[zx;zy] + Wu*1);
% vy = -zy + Wr(:,2)'*actF(Wz*[zx;zy] + Wu*1);
hold on;
quiver(X,Y, Vx, Vy,2, 'k', 'LineWidth', 2);
colormap viridis(100);
plot(cos(2*pi*(1:200)/200), sin(2*pi*(1:200)/200));
axis equal;
xlim([-1 1])
ylim([-1 1])
% quiver(X,Y, sign(Vx).*log(abs(Vx)), sign(Vy).*log(abs(Vy)), 5);
% figure;
% hold on;
% scatter(X(:),Y(:),10,log(sqrt(Vx(:).^2+Vy(:).^2)), 'filled');
% scatter(cos(2*pi*(1:100)/100), sin(2*pi*(1:100)/100), 10, 'k', 'filled');
% axis equal;
% xlim([-1 1])
% ylim([-1 1])
