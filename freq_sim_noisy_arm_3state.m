function sim = freq_sim_noisy_arm_3state(L,params,Nreps)

rng(1); % set seed for phase randomization

% parameters for A and B matrices
if(nargin<=2)
    G = 0.14;       % Viscous Constant: Ns/m; originally 0.14
    I = 0.1;        % Inertia Kgm2; originally 0.1
    tau = 0.066;    % Muscle time constant, s; originally 0.066
else
    G = params(1); I = params(2); tau = params(3);
end

delt = 0.005; % time step length in secs
sigma = 150*sqrt(delt); % execution noise
t = 0:delt:20-delt; % x axis for graphs

% sum-of-sines target movement
nstep = ceil(23/delt); % number of time steps

% create state space model in discrete time
A = [0 1 0
    0 -G/I 1/I
    0 0 -1/tau];
B = [0 0 1/tau]';

Ad = expm(A*delt);
order = size(Ad,1); % order of the system
Ad = [Ad zeros(order)
    zeros(order) Ad];

Bd = delt*B;
Bd = [Bd zeros(order,1)
    zeros(order,1) Bd];

% generate target movement
freqs = 0.05*primes(45); % target frequencies
freqX = freqs(1:2:end)'; 
freqY = freqs(2:2:end)';
amp = 0.015*ones(1,length(freqX)); % amplitudes
phase = 2*pi*rand(length(freqX),2)-pi; % phases of sum of sines
target(1,:,:) = amp*sin(freqX*2*pi*(0:delt:nstep*delt-delt) + repmat(phase(:,1),1,nstep));
phase = 2*pi*rand(length(freqX),2)-pi;
target(2,:,:) = amp*sin(freqY*2*pi*(0:delt:nstep*delt-delt) + repmat(phase(:,2),1,nstep));
target = squeeze(sum(target,2));

%state vector
x = zeros(2*order,Nreps,nstep);
u = zeros(size(Bd,2),Nreps,nstep); % movement commands
hand = zeros(2,Nreps,nstep);
x([1 4],:,1) = hand(:,:,1) - repmat(target(:,1),[1 Nreps]); % x position

% rng(i_seed); % set seed for noise during execution
% simulate trajectory
for i = 2:nstep
    u(:,:,i) = -L*x(:,:,i-1);
    x(:,:,i) = Ad*x(:,:,i-1) + Bd*(u(:,:,i)+sigma*randn(1,Nreps));
    
    % calculate absolute hand position
    hand([1 2],:,i) = hand([1 2],:,i-1) + (x([1 4],:,i) - x([1 4],:,i-1));
    x([1 4],:,i) = hand(:,:,i) - repmat(target(:,i),[1 Nreps]);
end

% compute fourier transforms
e = 2/delt; % remove first two seconds of tracking data for cleaner spectra
delay = 0.4/delt; % modeling sensorimotor delay
start = e+delay+1;
stop = (20/delt);

traj = hand(:,:,(e+1):(e+1)+stop);
traj(:,Nreps+1,:) = target(:,start:start+stop);
traj_avg = mean(traj,3);
traj_fft = fft(traj - repmat(traj_avg,[1 1 length(traj)]),[],3);

% find indices of target frequencies
idx = find(abs(traj_fft(1,Nreps+1,:))>1);
idx = idx(idx<2000);
idy = find(abs(traj_fft(2,Nreps+1,:))>1);
idy = idy(idy<2000);

% complex ratios
ratio = traj_fft(1,1:Nreps,idx)./repmat(traj_fft(1,Nreps+1,idx),[1 Nreps]); % XX
ratio(2,:,:) = traj_fft(2,1:Nreps,idy)./repmat(traj_fft(2,Nreps+1,idy),[1 Nreps]); % YY
ratio(3,:,:) = traj_fft(2,1:Nreps,idx)./repmat(traj_fft(1,Nreps+1,idx),[1 Nreps]); % XY
ratio(4,:,:) = traj_fft(1,1:Nreps,idy)./repmat(traj_fft(2,Nreps+1,idy),[1 Nreps]); % YX

sim.t = t;
sim.hand = hand;
sim.target = target;
sim.idx = [idx idy];
sim.ratio = ratio;
sim.freq = [freqX freqY];
sim.L = L;