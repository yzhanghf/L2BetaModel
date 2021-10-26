function [beta_est, beta_history, epsilon, iter_actual] = Fast_Newton_method_deg(d_vec, beta_0, lambda, epsilon_0, IterMax)
	
	d_tab = tabulate(d_vec);
	d_tab(d_tab(:,2)==0,:) = [];
	[~,ord] = sort(d_tab(:,1));
	d_tab = d_tab(ord,:);
	
	unique_degs = d_tab(:,1);
	deg_freq = d_tab(:,2); % vector of n_k's
	deg_prop = d_tab(:,3)/100.0;
	
	if(unique_degs(1)==0)
		unique_degs(1) = [];
		deg_freq(1) = [];
		deg_prop(1) = [];
	end
	
	m = length(unique_degs);
	
	dictionary = {};
	for(ii = 1:m)
		dictionary{ii} = find(d_vec==unique_degs(ii));
	end
	dictionary_0 = find(d_vec==0);
	
	n = length(d_vec);
	
	delta_0 = zeros(m,1);
	for(ii = 1:m)
		delta_0(ii) = mean(beta_0(dictionary{ii}));
	end
	delta_old = delta_0(:);
	
	epsilon = 2*epsilon_0;
	count = 1;
	
	delta_history = delta_old;
	
	% start iteration
	while( (count <= IterMax) & (epsilon > epsilon_0) )
		
		avg_delta_old = sum(delta_old .* deg_prop);
		
		FMatHat = delta_old * ones(1,m);
		FMatHat = FMatHat + FMatHat';  FMatHat = 1./(1+exp(-FMatHat));
		FMatHat = FMatHat - diag(diag(FMatHat));
		
		% new gradient vector
		G = deg_freq .* (FMatHat*deg_freq) + deg_freq .* (deg_freq-1) .* 1./(1+exp(-2*delta_old)) - deg_freq .* unique_degs + deg_freq .* (delta_old - avg_delta_old) * lambda;
		
		% new Jacobian matrix
		HMatHat = delta_old * ones(1,m);
		HMatHat = HMatHat + HMatHat';  HMatHat = 1./(1+exp(-HMatHat));  HMatHat = HMatHat .* (1-HMatHat);
		
		HMatHat_diag = diag(diag(HMatHat));
		HMatHat_offdiag = HMatHat - HMatHat_diag;
		
		Jacobian = diag(deg_freq) * HMatHat_offdiag * diag(deg_freq);
		Jacobian = Jacobian + diag(sum(Jacobian));
		temp2    = lambda * deg_freq * deg_freq' / n;
		temp2    = temp2 - diag(diag(temp2));
		Jacobian = Jacobian - temp2 + diag(deg_freq .*(1-deg_freq/n) * lambda) + 2*diag(deg_freq.*(deg_freq-1)) .* HMatHat_diag;
		% regularize
		Jacobian = Jacobian + eye(size(Jacobian,1));
		
		% update
		delta_new = delta_old - Jacobian\G;
		
		% compute epsilon
		epsilon = norm(sqrt(deg_prop).*(delta_new - delta_old)) / norm(sqrt(deg_prop).*delta_old);
		
		
		delta_old = delta_new;
		count = count + 1;
		delta_history(:,count) = delta_old;
		
	end
	
	delta_est = delta_new;
	beta_est = zeros(n,1);
	
	beta_history = [];
	
	for(ii = 1:m)
		beta_est(dictionary{ii}) = delta_est(ii);
		beta_history(dictionary{ii},:) = ones(length(dictionary{ii}),1)*delta_history(ii,:);
	end
	% directly assign values to zero degree nodes for numerical stability, this part can be tweaked:
	beta_est(dictionary_0) = 2*delta_est(1);
	% - 0.9*log(n);
	beta_history(dictionary_0,:) = 2*delta_est(1);
	% - 0.9*log(n);
	
	iter_actual = count;

end % end function