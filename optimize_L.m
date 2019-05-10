% optimize gain
clear all

load sim_data
Linit = dat.Linit;
params = [0.14 0.1 0.066];
f_targ = @(L) sim_error(L,params);
Nfreq = 7;
options = optimoptions('fmincon','Display','iter','MaxIterations',10);
Lopt = fmincon(f_targ,Linit,[],[],[],[],[],[],[],options);

% options = optimoptions('fminunc','Display','iter','MaxFunctionEvaluations',50000,'MaxIterations',3000,'Algorithm','quasi-newton');
% 
% Lopt = fminunc(f_targ,Linit,options);

global ratio_opt

% generate colormap for plotting
col1 = [1 0.83 0.33];
col2 = [1 0 0];
colors = [linspace(col1(1),col2(1),Nfreq)', linspace(col1(2),col2(2),Nfreq)', linspace(col1(3),col2(3),Nfreq)'];

% calculate mean of empirical data
names = {'x','y','xy','yx'};
for i = 1:length(names)
    a = mean(dat.(names{i}));
    mu(1,:,i) = real(a);
    mu(2,:,i) = imag(a);
end

% plot post-fitting simulation data
subplot(2,4,5); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
axis square
ylabel('After')

subplot(2,4,6); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
axis square

subplot(2,4,7); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
axis square

subplot(2,4,8); hold on
plot([-1.5 1.5],[0 0],'k')
plot([0 0],[-1.5 1.5],'k')
axis square

for i = 1:Nfreq
    subplot(2,4,5)
    plot(ratio_opt(i,:,1),'.','Color',colors(i,:),'MarkerSize',5)
    error_ellipse(cov(real(dat.x(:,i)),imag(dat.x(:,i))),mu(:,i,1),'conf',0.95,'style','k');
    
    subplot(2,4,6)
    plot(ratio_opt(i,:,2),'.','Color',colors(i,:),'MarkerSize',5)
    error_ellipse(cov(real(dat.y(:,i)),imag(dat.y(:,i))),mu(:,i,2),'conf',0.95,'style','k');
    
    subplot(2,4,7)
    plot(ratio_opt(i,:,3),'.','Color',colors(i,:),'MarkerSize',5)
    error_ellipse(cov(real(dat.xy(:,i)),imag(dat.xy(:,i))),mu(:,i,3),'conf',0.95,'style','k');
    
    subplot(2,4,8)
    plot(ratio_opt(i,:,4),'.','Color',colors(i,:),'MarkerSize',5)
    error_ellipse(cov(real(dat.yx(:,i)),imag(dat.yx(:,i))),mu(:,i,4),'conf',0.95,'style','k');
end