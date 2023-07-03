function RNN = trainLowRankRNN(dt, trainInfo, latent, u)

    actF = @(x) tanh(x);

    %%
    sidx = find(u==1,1);
    preTime = -dt*sidx;
    t = (0:numel(u)-1)*dt + preTime;
    T = numel(t);

    %%
    latent1 = latent(:,1);
    latent2 = latent(:,2);
    ang = atan2(latent2, latent1);

    %%
    tau = dt;
    alpha = trainInfo.alpha;
    itr = trainInfo.itr;
    N = trainInfo.N;
    p = eye(N)/alpha;

    J = zeros(N,N);
    g = 0;
    
    x = zeros(N,T);
    
    Wz = randn(N,2);
    Wu = randn(N,1);
    
    Wr = zeros(N,2);
    
    z = [cos(ang)';sin(ang)'] + 2*(rand(2,numel(ang))-0.5)*trainInfo.bpNoiseLevel;

    for m = 1:itr
        x(:,1) = Wz*mean(z(:,1:sidx-1),2);
        x(:,2:end) = 0;
        for k = 1:numel(t)-1
            % update recursive least square(RLS) states - FORCE learning
%             dx = (dt./tau).*(-x(:,k) + g*J*actF(x(:,k)) + Wz*z(:,k) + Wu*u(k) + randn(N,1)/sqrt(N));
            dx = (dt./tau).*(-x(:,k) + g*J*actF(x(:,k)) + Wz*z(:,k) + randn(N,1)/sqrt(N));
            x(:,k+1) = x(:,k) + dx;
    
            r = actF(x(:,k+1));
            p = p - p*(r*r')*p/(1 + r'*p*r);
            eMinus(:,k) = Wr'*r - z(:, k+1);
            Wr = Wr - p*r*eMinus(:,k)';
            %ePlus(:,k) = Wr'*r - z(:,k+1);
        end
        %fprintf(['Iteration: ', num2str(m), '/', num2str(itr), '\n']);
    end

    %%
    RNN.dt = dt;
    RNN.u = u;
    RNN.Wr = Wr;
    RNN.Wz = Wz;
    RNN.Wu = Wu;
    RNN.z = z;
    RNN.trainInfo = trainInfo;
end