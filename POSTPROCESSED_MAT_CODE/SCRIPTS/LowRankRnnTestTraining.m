clc;
clear;
close all;

% 1st order kinetics calculator
firstOrderKin = @(x0, xInf, t, tau) repmat(xInf,1,numel(t))-repmat((xInf-x0),1,numel(t)).*exp(-(t-t(1))/tau);

% function: population state phase (angle) to the population state in latent space (Z)
phiToZ = @(phi) [cos(phi);sin(phi)];

% activation fucntion
actF = @(x) tanh(x);


%% Parameter Initialization
dt = 0.1;% time step [s]
Tpre = 30;% pre-stimulation duration
Tstim = 60;% stimulation duration
Tpost = 30;% post-stimulation duration
T = Tpre+Tstim+Tpost; % total time
nT = round(T/dt);% total sample length

N = 500;% number of neurons in a network
Rank = 1; % rank of network
tau = 0.1;% passive time constant of neuron membrane[s]

x = zeros(N, nT); % zero initialization of network state (biophysically means the membrane current of a neuron)

g = 0.0;% gain of heterogeneous connectivity
J = randn(N)/sqrt(N); % random heterogeneous connectivity

Wr = zeros(N,Rank); % projection (neural population in the network -> latent space)  = synaptic weight from original neural population to artificial neurons in latent space

Wz = randn(N,2); % projection (latent space -> neural population) = synaptic weight from artificial neuron in latent space to original neural population

Wu = randn(N,1); % synaptic weight for external stimulus input

%% Training data preparation
tPre = dt*(1:round(Tpre/dt)); 
tStim = dt*(round(Tpre/dt)+1:round((Tpre+Tstim)/dt));
tPost = dt*(round((Tpre+Tstim)/dt)+1:round((Tpre+Tstim+Tpost)/dt));
t = [tPre, tStim, tPost];

% pre-determined neural trajectory in latent space
phiOff = -1.3*pi; % pre-determined population state phase in resting
phiOn = 0.3*pi; % pre-determined population state phase during stimulation

% dynamics of population state phase
phiPre = phiOff*ones(1,numel(tPre)); % resting (pre-stimulation)
phiStim = firstOrderKin(phiOff, phiOn, tStim, 7); % during stimulation
phiPost = firstOrderKin(phiOn, 2*pi+phiOff, tPost, 7); % recovery (post-stimulation)
phi = [phiPre, phiStim, phiPost];

% dynamics of population state in latent space
zPre = phiToZ(phiPre);
zStim = phiToZ(phiStim);
zPost = phiToZ(phiPost);
z = [zPre, zStim, zPost];
trainDataFig = figure; trainDataFig.Position(1) = 70;
subplot(1,2,1);
scatter(z(1,:), z(2,:), 50, t);
title('Neural trajectory in latent space');
xlabel('PC_1');
ylabel('PC_2');
subplot(1,2,2);
scatter(t, z(2,:), 50, t);
title('Neural decoding of latent space');
xline(Tpre, '-', {'Stim-On'}, 'Color', 'r');
xline(Tpre+Tstim, '-', {'Stim-Off'}, 'Color', 'r');
ylabel('BP dynamics (a.u)');
xlabel('Time (s)');

% external input
uPre = zeros(1, numel(tPre));
uStim = ones(1, numel(tStim)); % stimulation normalized as 1
uPost = zeros(1, numel(tPost));
u = [uPre, uStim, uPost];

%% Training the neural network 
alpha = 5;
p = eye(N)/alpha;

itr = 20; % iteration

for j = 1:itr
    x(:,1) = x(:,end);
    x(:,2:end) = 0;
    for i = 1:numel(t)-1
        % update recursive least square(RLS) states - FORCE learning
        dx = (dt./tau).*(-x(:,i) + g*J*actF(x(:,i)) + Wz*z(:,i) + Wu*u(i));
        x(:,i+1) = x(:,i) + dx;

        r = actF(x(:,i+1));
        p = p - p*(r*r')*p/(1 + r'*p*r);
        eMinus(:,i) = Wr'*r - z(:, i+1);
        Wr = Wr - p*r*eMinus(:,i)';
        ePlus(:,i) = Wr'*r - z(:,i+1);
    end
    fprintf(['Iteration: ', num2str(j), '/', num2str(itr), '\n']);
end

%% Results analyses 
%mtx = [Wz, N*Wr, Wu];
%covmtx = mtx'*mtx/N;
run('LowRankRnnPlot.m');