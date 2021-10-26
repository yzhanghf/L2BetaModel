
% L2 error rate as n changes

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



gradient_time_mat = [];
gradient_error_mat = [];
gradient_relative_error_mat = [];
Newton_time_mat = [];
Newton_error_mat = [];
Newton_relative_error_mat = [];


font_size = 20;

Nexp = 100;


for(n_index = 1:length(n_vec))
	
	n = n_vec(n_index);
	fprintf(1,'n=%d:\n',n);
	
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
		% beta_0 = zeros(n,1);
		d_vec = sum(A);
		beta_0 = log(d_vec+1)-mean(log(d_vec+1));
		rho_n = sum(A(:))/(n*(n-1));
		stepsize = 1/(2*(n-1)+lambda);
		
		FileNameString1 = sprintf('Setting_%d',  TrueBetaSetting);
		FileNameString2 = sprintf('n_%d',  n);
		
		tic;
		beta_est = Fast_gradient_deg(sum(A), beta_0, lambda, stepsize, epsilon_0, IterMax,1);
		gradient_time_mat(n_index, exp_index) = toc;
		
		gradient_error_mat(n_index, exp_index) = ...
			norm(beta_est-beta_true)/sqrt(n);
		gradient_relative_error_mat(n_index, exp_index) = ...
			norm(beta_est-beta_true) / norm(beta_true);
			
		
		tic;
		beta_0 = zeros(n,1);
		beta_est = Fast_Newton_method_deg(sum(A), beta_0, lambda, epsilon_0, IterMax);
		Newton_time_mat(n_index, exp_index) = toc;
		
		Newton_error_mat(n_index, exp_index) = ...
			norm(beta_est-beta_true)/sqrt(n);
		Newton_relative_error_mat(n_index, exp_index) = ...
			norm(beta_est-beta_true) / norm(beta_true);
		
		
		
		
	end
	
	fprintf(1,'\n');
	
	dlmwrite(...
		sprintf('./numerical_results/%s_gradient_time_%s.txt',FileNameString1,FileNameString2),...
		gradient_time_mat(n_index, :)'...
	);
	dlmwrite(...
		sprintf('./numerical_results/%s_gradient_error_%s.txt',FileNameString1,FileNameString2),...
		gradient_error_mat(n_index, :)'...
	);
	dlmwrite(...
		sprintf('./numerical_results/%s_gradient_relative_error_%s.txt',FileNameString1,FileNameString2),...
		gradient_relative_error_mat(n_index, :)'...
	);
	
	dlmwrite(...
		sprintf('./numerical_results/%s_Newton_time_%s.txt',FileNameString1,FileNameString2),...
		Newton_time_mat(n_index, :)'...
	);
	dlmwrite(...
		sprintf('./numerical_results/%s_Newton_error_%s.txt',FileNameString1,FileNameString2),...
		Newton_error_mat(n_index, :)'...
	);
	dlmwrite(...
		sprintf('./numerical_results/%s_Newton_relative_error_%s.txt',FileNameString1,FileNameString2),...
		Newton_relative_error_mat(n_index, :)'...
	);

end









