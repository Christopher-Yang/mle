function e2 = sim_error(L,params)

% load single_subj
load sim_data
persistent graph % used for graphing pre-fitting simulation data
global ratio_opt % saves optimal ratio for graphing

% simulate multiple times with noise to get a noisy phaser plot
if(nargin<2)
    params = [0.15 0.11 0.06];
else
    params = params;
end

Nreps = 50; % number of simulations to run
sim = freq_sim_noisy_arm_3state(L,params,Nreps); % perform simulations

Nfreq = length(sim.freq);
ratio = permute(sim.ratio,[3 2 1]);
ratio_opt = ratio;

% estimate Gaussian distribution for real and complex parts
names = {'x','y','xy','yx'};
z = NaN([2 Nreps Nfreq length(names)]);
z_mean = NaN([2 Nfreq length(names)]);
zN = NaN([2 Nreps Nfreq length(names)]);
Sigma = NaN([2 2 Nfreq length(names)]);
p = NaN(Nfreq,2);

for j = 1:length(names)
    z_data(1,:,:,j) = real(dat.(names{j}));
    z_data(2,:,:,j) = imag(dat.(names{j}));
    for i=1:Nfreq
        z(1,:,i,j) = real(ratio(i,:,j));
        z(2,:,i,j) = imag(ratio(i,:,j));
        
        % get mean
        z_mean(:,i,j) = mean(z(:,:,i,j),2);
        
        % get covariance matrix
        zN(:,:,i,j) = z(:,:,i,j)-repmat(z_mean(:,i,j),[1 Nreps]);
        Sigma(:,:,i,j) = zN(:,:,i,j)*zN(:,:,i,j)'/Nreps;
        for k = 1:size(z_data,2)
            p(k,i,j) = exp(-0.5*(z_data(:,k,i,j)-z_mean(:,i,j))'*inv(Sigma(:,:,i,j))*(z_data(:,k,i,j)-z_mean(:,i,j)))/sqrt(((2*pi)^2)*det(Sigma(:,:,i,j)));
        end
    end
    % calculate mean of empirical data
    a = mean(dat.(names{j}));
    mu(1,:,j) = real(a);
    mu(2,:,j) = imag(a);
end

% negative log likelihood
e2 = sum(sum(sum(-log(p))));

% only plot the pre-fitting phasors
if isempty(graph)
    col1 = [1 0.83 0.33];
    col2 = [1 0 0];
    colors = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

    figure(1); clf
    subplot(2,4,1); hold on
    plot([-1.5 1.5],[0 0],'k')
    plot([0 0],[-1.5 1.5],'k')
    axis square
    title('XX')
    ylabel('Before')
    
    subplot(2,4,2); hold on
    plot([-1.5 1.5],[0 0],'k')
    plot([0 0],[-1.5 1.5],'k')
    axis square
    title('YY')
    
    subplot(2,4,3); hold on
    plot([-1.5 1.5],[0 0],'k')
    plot([0 0],[-1.5 1.5],'k')
    axis square
    title('XY')
    
    subplot(2,4,4); hold on
    plot([-1.5 1.5],[0 0],'k')
    plot([0 0],[-1.5 1.5],'k')
    axis square
    title('YX')
    
    for i = 1:Nfreq
        subplot(2,4,1)
        plot(ratio(i,:,1),'.','Color',colors(i,:),'MarkerSize',5)
        error_ellipse(cov(real(dat.x(:,i)),imag(dat.x(:,i))),mu(:,i,1),'conf',0.95,'style','k');
        
        subplot(2,4,2)
        plot(ratio(i,:,2),'.','Color',colors(i,:),'MarkerSize',5)
        error_ellipse(cov(real(dat.y(:,i)),imag(dat.y(:,i))),mu(:,i,2),'conf',0.95,'style','k');
        
        subplot(2,4,3)
        plot(ratio(i,:,3),'.','Color',colors(i,:),'MarkerSize',5)
        error_ellipse(cov(real(dat.xy(:,i)),imag(dat.xy(:,i))),mu(:,i,3),'conf',0.95,'style','k');
        
        subplot(2,4,4)
        plot(ratio(i,:,4),'.','Color',colors(i,:),'MarkerSize',5)
        error_ellipse(cov(real(dat.yx(:,i)),imag(dat.yx(:,i))),mu(:,i,4),'conf',0.95,'style','k');
    end
    graph=0;
end