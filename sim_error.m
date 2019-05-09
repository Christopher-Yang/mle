function e2 = sim_error(L,params)

% load real_data;

% load 3state;
% load single_subj
load sim_data
persistent graph
global ratio_opt
global z_opt
global Sigma_opt

% simulate multiple times with noise to get a noisy phaser plot
if(nargin<2)
    params = [0.15 0.11 0.06];
else
    params = params;
end

REAL_DATA = 0;
Nreps = 10;
block = 1;
%L = [80 9 5];

% for i=1:Nreps
sim = freq_sim_noisy_arm_3state(L,params,Nreps);
% end
Nfreq = length(sim.freq);
ratio = permute(sim.ratio,[3 2 1]);
ratio_opt = ratio;

% for i=1:Nreps
%     for j=1:Nfreq
%         ratio(i,j) = sim(i).ratio(j);
%     end
% end

% for i=1:length(real_data.fx)
%     f_ind(i) = find(abs(sim(1).freq - real_data.fx(i))<.001); % find corresponding frequency for data in the simulation
% end

% estimate Gaussian distribution for real and complex parts
names = {'x','y','xy','yx'};
z = NaN([2 Nreps Nfreq length(names)]);
z_mean = NaN([2 Nfreq length(names)]);
zN = NaN([2 Nreps Nfreq length(names)]);
Sigma = NaN([2 2 Nfreq length(names)]);
p = NaN(Nfreq,2);
% z_data = cat(3,[real(real_data.x_fft(1,:)); imag(real_data.x_fft(1,:))],[real(real_data.x_fft(1,:)); imag(real_data.x_fft(1,:))]);
% z_data = cat(3,[real(mean(dat.x,2)'); imag(mean(dat.x,2)')], [real(mean(dat.y,2)'); imag(mean(dat.y,2)')]);

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
        
%         subplot(1,2,j)
%         plot(z_mean(1,i,j),z_mean(2,i,j),'r.')
%         plot_ellipse(z_mean(:,i,j),Sigma(:,:,i,j));
    end
    z_opt = z_mean;
    Sigma_opt = Sigma;
    a = mean(dat.(names{j}));
    mu(1,:,j) = real(a);
    mu(2,:,j) = imag(a);
end

% if REAL_DATA
%     z_data(1,:,1) = real(real_data.x_fft(block,:));
%     z_data(2,:,1) = imag(real_data.x_fft(block,:));
%     z_data(1,:,2) = real(real_data.y_fft(block,:));
%     z_data(2,:,2) = imag(real_data.y_fft(block,:));
% else
%     z_data(1,:) = real(d.ratio);
%     z_data(2,:) = imag(d.ratio);
% end
% e2 = sum(sum(sum((z_data-z_mean).^2)));
e2 = sum(sum(sum(-log(p))));

% plot data
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
%     if REAL_DATA
%         subplot(2,2,1)
%         plot(dat.x','k.','LineWidth',1.5)
%         subplot(2,2,2)
%         plot(dat.y','k.','LineWidth',1.5)
%     else
%         subplot(2,2,1)
%         for i = 1:Nfreq
%             error_ellipse(cov(real(dat.x(:,i)),imag(dat.x(:,i))),mu(:,i,1),'conf',0.95,'style','k');
%         end
%         subplot(2,2,2)
%         for i = 1:Nfreq
%             error_ellipse(cov(real(dat.y(:,i)),imag(dat.y(:,i))),mu(:,i,2),'conf',0.95,'style','k');
%         end
%     end
    graph=0;
end

%error = z_mean - 

%polarplot(data.rot.avg.x_x.fft(1,:),'-o','LineWidth',1.5);
% polarplot(data.rot.avg.x_x.fft(2,:),'-o','LineWidth',1.5);
% polarplot(data.rot.avg.x_x.fft(5,:),'-o','LineWidth',1.5);

% legend('Inf Horizon','X -> X data','Y -> Y data');
% legend('.0001x effort','.01x','Default','100x','10000x','X -> X data','Y -> Y data');
% legend('10 ms delay','100 ms','200 ms (default)','400 ms','600 ms','X -> X data','Y -> Y data');