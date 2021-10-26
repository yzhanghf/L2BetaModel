
% Automatic tuning of \lambda

rng(2608);

LambdaRange = exp([0:0.5:4])-1;
nlambda = length(LambdaRange);

switch TrueBetaSetting
	
	case 1
		gamma = -1/3;  alpha = 1/2;  s_prop = 1/10;
		n = 500;  epsilon_0 = 1e-5;  IterMax = 100;
	
	case 4
		gamma = -1/3;  alpha = -1/3+0.05;  s_prop = 1/3;
		n = 500;  epsilon_0 = 1e-5;  IterMax = 300;
	
end

s = floor(s_prop * n);
beta_true = zeros(n,1);
	beta_true(1:s) = alpha*log(n);
	beta_true((s+1):n) = gamma*log(n);
W = beta_true*ones(1,n);  W = W + W';  W = 1./(1+exp(-W));  W = W-diag(diag(W));


A = generate_A(W);
deg_obs = sum(A);  deg_obs = deg_obs(:);


gradient_TuningCriterion = zeros(nlambda,1);
gradient_TrueL2Loss = zeros(nlambda,1);

for(lambda_index = 1:nlambda)
	
	fprintf(1,'%d ',lambda_index);
	lambda = LambdaRange(lambda_index);
	stepsize = 1/((n-2)/2+lambda);
	
	% Run the algorithm once and compute the true estimation error
	beta_0 = zeros(n,1);
	[beta_est, beta_history, epsilon, iter_actual] = ...
			Gradient_descend(A, beta_0, lambda, stepsize, epsilon_0, IterMax);
	gradient_TrueL2Loss(lambda_index) = norm(beta_est-beta_true)/norm(beta_true);
	
	
	% Tuning method: thresholding ridge estimation + AIC
	eb = beta_est*ones(1,n);  eb = eb + eb';  eb = exp(eb);  eb = eb - diag(diag(eb));
	NegLogLike = sum(sum(log(1+eb)))/2 - sum(beta_est .* deg_obs);
	AIC = (n-2)*max(deg_obs)/(max(deg_obs)*(n-2)/(n-1)+lambda) + 2*max(deg_obs)/(2*max(deg_obs)+lambda);
	
	gradient_TuningCriterion(lambda_index) = AIC + NegLogLike;
	
end

fprintf(1,'\n');

print_gradient_TuningCriterion = (gradient_TuningCriterion - min(gradient_TuningCriterion)) / range(gradient_TuningCriterion) * max(gradient_TrueL2Loss);


font_size = 18;  MarkerSize = 15;  LineWidth = 2;

fig = figure('visible','off');
plot1 = plot(log(LambdaRange+1), print_gradient_TuningCriterion, 'k--o');
hold on;
plot2 = plot(log(LambdaRange+1), gradient_TrueL2Loss, 'r-*');
xlim([1 log(max(LambdaRange)+1)]);  xticks(log(LambdaRange+1));  xticklabels(round(LambdaRange,2));
xlabel('Penalty scale $\lambda$', 'interpreter','latex', 'FontSize', font_size);
ylabel('Estimation error; AIC($\lambda$)','interpreter','latex', 'FontSize', font_size);
plot1.MarkerSize = MarkerSize;  plot2.MarkerSize = MarkerSize;
plot1.LineWidth = LineWidth;    plot2.LineWidth = LineWidth;
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',font_size);
saveas(fig, sprintf('./figures/tuning_lambda_Setting_%d.png',TrueBetaSetting));



