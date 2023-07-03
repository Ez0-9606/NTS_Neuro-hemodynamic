function [xNew, xNew2] = hBINDsimulate(N, dt, tau, g, J, actF, t, z, u, Wz, Wu, Wr)

    %%
    
    tPre = t(1:find(u==1, 1)-1);

    xNew = zeros(N, numel(t));
    xNew2 = zeros(N, numel(t));

    f2 = figure;
    for j = 1:5
        if j == 1
            xNew(:,1) = Wz*mean(z(:,1:numel(tPre)),2);
            xNew2(:,1) = Wz*mean(z(:,1:numel(tPre)),2);
        else
            xNew(:,1) = xNew2(:,end);
            xNew2(:,1) = xNew2(:,end);
        end

        
        for i = 1:numel(t)-1
            % update network states
            try
                dx = (dt./tau).*(-xNew(:,i) + g*J*actF(xNew(:,i)) + Wz*Wr'*actF(xNew(:,i)) + 1*Wu*u(i) + randn(N,1)/sqrt(N));
            catch
                disp(j);
                disp(i);
                disp(t(i));
            end
            xNew(:,i+1) = xNew(:,i) + dx;
        
            dx2 = (dt./tau).*(-xNew2(:,i) + g*J*actF(xNew2(:,i)) + Wz*Wr'*actF(xNew2(:,i)) + 0*Wu*u(i) + randn(N,1)/sqrt(N));
            xNew2(:,i+1) = xNew2(:,i) + dx2;
        end
        
        %%
        figure(f2);
        subplot(2,1,1); hold on;
        plot(t+(j-1)*dt*numel(t), z(1,:), 'Color', 'k');
        plot(t+(j-1)*dt*numel(t), Wr(:,1)'*actF(xNew), 'Color', 'r');
        plot(t+(j-1)*dt*numel(t), Wr(:,1)'*actF(xNew2), 'Color', 'b');
        plot(t+(j-1)*dt*numel(t), u, 'Color', 'g');
        subplot(2,1,2); hold on;
        plot(t+(j-1)*dt*numel(t), z(2,:), 'k');
        plot(t+(j-1)*dt*numel(t), Wr(:,2)'*actF(xNew), 'Color', 'r');
        plot(t+(j-1)*dt*numel(t), Wr(:,2)'*actF(xNew2), 'Color', 'b');
        plot(t+(j-1)*dt*numel(t), u, 'Color', 'g');
    end

end