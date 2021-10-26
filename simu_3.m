
% Validation of normality

rng(2608);

switch(TrueBetaSetting)
	
	case 1
		gamma = -1/3;  alpha = 1/5;  s_prop = 1/5;
		lambda = 0.1;  epsilon_0 = 1e-7;  IterMax = 500;
	
	case 2
		gamma = -1/2;  alpha = 1/5;  s_prop = 1/5;
		lambda = 0.1;  epsilon_0 = 1e-7;  IterMax = 500;
		
	case 3
		gamma = -2/3;  alpha = 1/3;  s_prop = 1/5;
		lambda = 0.1;  epsilon_0 = 1e-7;  IterMax = 500;
		
	case 4
		gamma = -2/3;  alpha = 1/3;  s_prop = 1/5;
		lambda = 10;  epsilon_0 = 1e-7;  IterMax = 500;
	
	case 5
		gamma = -2/3;  alpha = 1/3;  s_prop = 1/5;
		lambda = 200;  epsilon_0 = 1e-7;  IterMax = 500;
		
	case 6
		gamma = -1/3;  alpha = gamma+0.05;  s_prop = 1/5;
		lambda = 200;  epsilon_0 = 1e-7;  IterMax = 500;
		
		
end
		

n_vec = [50,100,200,400,800,1600,3200,6400,12800];

for n = n_vec

fprintf(1,'\ncase=%d, n=%d\n',TrueBetaSetting,n);

beta_list_Gradient = [];
beta_list_Newton   = [];

for(exp_index = 1:Nexp)
	
	fprintf(1,'%d ',exp_index);
	
	s = floor(s_prop * n);
	
	%%% Set up beta parameters
	beta_true = zeros(n,1);
		beta_true(1:s) = alpha*log(n);
		beta_true((s+1):n) = gamma*log(n);

	%%% Generate network
	W = beta_true*ones(1,n);  W = W + W';  W = 1./(1+exp(-W));  W = W-diag(diag(W));
	A = generate_A(W);
	% imagesc(A);

	%%% Gradient method

	% initial guess
	d_vec = sum(A);
	beta_0 = log(d_vec+1)-mean(log(d_vec+1));
	rho_n = sum(A(:))/(n*(n-1));
	stepsize = 1/(2*(n-1)+lambda);
	
	FileNameString1 = sprintf('Setting_%d',  TrueBetaSetting);
	FileNameString2 = sprintf('n_%d',  n);
	
	[beta_est] = ...
		Fast_gradient_deg(sum(A), beta_0, lambda, stepsize, epsilon_0, IterMax,1);
	beta_list_Gradient(:, exp_index) = beta_est;
	
	
	
	[beta_est] = ...
		Fast_Newton_method_deg(sum(A), beta_0, lambda, epsilon_0, IterMax);
	
	beta_list_Newton(:, exp_index) = beta_est;
	
end

fprintf(1,'\n');

writematrix(beta_list_Gradient,...
	sprintf('./numerical_results/normality_gradient_case_%d_n_%d.txt', TrueBetaSetting,n));

writematrix(beta_list_Newton,...
	sprintf('./numerical_results/normality_Newton_case_%d_n_%d.txt', TrueBetaSetting,n));


true_location_vec = ones(1,n)*((n-2)/n*eye(n)+1/n*ones(n,n))*beta_true/(2*(n-1));



figure('visible','off');

font_size = 20;

data = beta_list_Gradient(1,:);  data = data(:);
pd = fitdist(data,'normal');
x_pdf = linspace(min(data),max(data));
y_pdf = pdf(pd,x_pdf);
fig = histogram(data,'Normalization','pdf');  hold on;
x1 = xline(beta_true(1),'g:',{'True $\beta_1$ ($\lambda=0$)'},'linewidth',2,'LabelVerticalAlignment','middle', 'interpreter','latex');
x2 = xline(true_location_vec(1),'m:',{'Mean-predict (large $\lambda$)'},'linewidth',2,'LabelVerticalAlignment','middle', 'interpreter','latex');
x1.FontSize = font_size;
x2.FontSize = font_size;
line(x_pdf,y_pdf,'LineWidth',2,'color','red','linestyle','--');   hold off;
title(sprintf('Histogram of $\\widehat\\beta_{\\lambda,1}$ (gradient)\nSetting %d, $\\lambda=%1.1f$, $n=%d$', TrueBetaSetting,lambda,n),'interpreter','latex');
set(gca,'FontSize',font_size);
set(gca,'color',[192,192,192]/255);
set(gcf, 'InvertHardCopy', 'off');
saveas(fig, sprintf('./figures/normality_Gradient_beta_1_case_%d_n_%d.png', TrueBetaSetting,n));


data = beta_list_Gradient(1:2,:);  data = data';
fig = scatter(data(:,1),data(:,2),'k');
title(sprintf('Scatter plot of $(\\widehat\\beta_{\\lambda,1}, \\widehat\\beta_{\\lambda,2})$ (gradient)\nSetting %d, $\\lambda=%1.1f$, $n=%d$', TrueBetaSetting,lambda,n),'interpreter','latex');
set(gca,'FontSize',font_size);
set(gca,'color',[192,192,192]/255);
set(gcf, 'InvertHardCopy', 'off');
saveas(fig, sprintf('./figures/normality_Gradient_beta_12_case_%d_n_%d.png', TrueBetaSetting,n));


data = beta_list_Newton(1,:);  data = data(:);
pd = fitdist(data,'normal');
x_pdf = linspace(min(data),max(data));
y_pdf = pdf(pd,x_pdf);
fig = histogram(data,'Normalization','pdf');  hold on;
x1 = xline(beta_true(1),'g:',{'True $\beta_1$ ($\lambda=0$)'},'linewidth',2,'LabelVerticalAlignment','middle', 'interpreter','latex');
x2 = xline(true_location_vec(1),'m:',{'Mean-predict (large $\lambda$)'},'linewidth',2,'LabelVerticalAlignment','middle', 'interpreter','latex');
x1.FontSize = font_size;
x2.FontSize = font_size;
line(x_pdf,y_pdf,'LineWidth',2,'color','red','linestyle','--');   hold off;
title(sprintf('Histogram of $\\widehat\\beta_{\\lambda,1}$ (Newton)\nSetting %d, $\\lambda=%1.1f$, $n=%d$', TrueBetaSetting,lambda,n),'interpreter','latex');
set(gca,'FontSize',font_size);
set(gca,'color',[192,192,192]/255);
set(gcf, 'InvertHardCopy', 'off');
saveas(fig, sprintf('./figures/normality_Newton_beta_1_case_%d_n_%d.png', TrueBetaSetting,n));

data = beta_list_Newton(1:2,:);  data = data';
fig = scatter(data(:,1),data(:,2),'k');
title(sprintf('Scatter plot of $(\\widehat\\beta_{\\lambda,1}, \\widehat\\beta_{\\lambda,2})$ (Newton)\nSetting %d, $\\lambda=%1.1f$, $n=%d$', TrueBetaSetting,lambda,n),'interpreter','latex');
set(gca,'FontSize',font_size);
set(gca,'color',[192,192,192]/255);
set(gcf, 'InvertHardCopy', 'off');
saveas(fig, sprintf('./figures/normality_Newton_beta_12_case_%d_n_%d.png', TrueBetaSetting,n));













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data = beta_list_Gradient(n,:);  data = data(:);  data = outlier_elim(data);
pd = fitdist(data,'normal');
x_pdf = linspace(min(data),max(data));
y_pdf = pdf(pd,x_pdf);
fig = histogram(data,'Normalization','pdf');  hold on;
x1 = xline(beta_true(n),'g:',{'True $\beta_n$ ($\lambda=0$)'},'linewidth',2,'LabelVerticalAlignment','middle', 'interpreter','latex');
x2 = xline(true_location_vec(1),'m:',{'Mean-predict (large $\lambda$)'},'linewidth',2,'LabelVerticalAlignment','middle', 'interpreter','latex');
x1.FontSize = font_size;
x2.FontSize = font_size;
line(x_pdf,y_pdf,'LineWidth',2,'color','red','linestyle','--');   hold off;
title(sprintf('Histogram of $\\widehat\\beta_{\\lambda,%d}$ (gradient)\nSetting %d, $\\lambda=%1.1f$, $n=%d$', n,TrueBetaSetting,lambda,n),'interpreter','latex');
set(gca,'FontSize',font_size);
set(gca,'color',[192,192,192]/255);
set(gcf, 'InvertHardCopy', 'off');
saveas(fig, sprintf('./figures/normality_Gradient_beta_n_case_%d_n_%d.png', TrueBetaSetting,n));


data = beta_list_Gradient((n-1):n,:);  data = data';  data = outlier_elim(data);
fig = scatter(data(:,1),data(:,2),'k');
title(sprintf('Scatter plot of $(\\widehat\\beta_{\\lambda,%d}, \\widehat\\beta_{\\lambda,%d})$ (gradient)\nSetting %d, $\\lambda=%1.1f$, $n=%d$', n-1,n,TrueBetaSetting,lambda,n),'interpreter','latex');
set(gca,'FontSize',font_size);
set(gca,'color',[192,192,192]/255);
set(gcf, 'InvertHardCopy', 'off');
saveas(fig, sprintf('./figures/normality_Gradient_beta_n-1n_case_%d_n_%d.png', TrueBetaSetting,n));


data = beta_list_Newton(n,:);  data = data(:);  data = outlier_elim(data);
pd = fitdist(data,'normal');
x_pdf = linspace(min(data),max(data));
y_pdf = pdf(pd,x_pdf);
fig = histogram(data,'Normalization','pdf');  hold on;
x1 = xline(beta_true(n),'g:',{'True $\beta_n$ ($\lambda=0$)'},'linewidth',2,'LabelVerticalAlignment','middle', 'interpreter','latex');
x2 = xline(true_location_vec(1),'m:',{'Mean-predict (large $\lambda$)'},'linewidth',2,'LabelVerticalAlignment','middle', 'interpreter','latex');
x1.FontSize = font_size;
x2.FontSize = font_size;
line(x_pdf,y_pdf,'LineWidth',2,'color','red','linestyle','--');   hold off;
title(sprintf('Histogram of $\\widehat\\beta_{\\lambda,%d}$ (Newton)\nSetting %d, $\\lambda=%1.1f$, $n=%d$', n,TrueBetaSetting,lambda,n),'interpreter','latex');
set(gca,'FontSize',font_size);
set(gca,'color',[192,192,192]/255);
set(gcf, 'InvertHardCopy', 'off');
saveas(fig, sprintf('./figures/normality_Newton_beta_n_case_%d_n_%d.png', TrueBetaSetting,n));

data = beta_list_Newton((n-1):n,:);  data = data';  data = outlier_elim(data);
fig = scatter(data(:,1),data(:,2),'k');
title(sprintf('Scatter plot of $(\\widehat\\beta_{\\lambda,%d}, \\widehat\\beta_{\\lambda,%d})$ (Newton)\nSetting %d, $\\lambda=%1.1f$, $n=%d$', n-1,n,TrueBetaSetting,lambda,n),'interpreter','latex');
set(gca,'FontSize',font_size);
set(gca,'color',[192,192,192]/255);
set(gcf, 'InvertHardCopy', 'off');
saveas(fig, sprintf('./figures/normality_Newton_beta_n-1n_case_%d_n_%d.png', TrueBetaSetting,n));


end
