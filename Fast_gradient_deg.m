function [beta_est, count, stepsize, delta_old, delta_new] = Fast_gradient_deg(d_vec, beta_0, lambda, stepsize, epsilon_0, IterMax, adaptive_stepsize)
	
	if(~exist('adaptive_stepsize','var'))
		adaptive_stepsize = 0;
	end
	
	d_tab = tabulate(d_vec);
	d_tab(d_tab(:,2)==0,:) = [];
	[~,ord] = sort(d_tab(:,1));
	d_tab = d_tab(ord,:);
	
	unique_degs = d_tab(:,1);
	deg_freq = d_tab(:,2); % vector of n_k's
	deg_prop = d_tab(:,3)/100.0;
	m = length(unique_degs);
	
	dictionary = {};
	for(ii = 1:m)
		dictionary{ii} = find(d_vec==unique_degs(ii));
	end
	
	n = length(d_vec);
	
	delta_0 = zeros(m,1);
	for(ii = 1:m)
		delta_0(ii) = mean(beta_0(dictionary{ii}));
	end
	delta_old = delta_0(:);
	
	
	epsilon = 2*epsilon_0;
	count = 1;
	
	% start iteration
	while( (count <= IterMax) & (epsilon > epsilon_0) )
		
		avg_delta_old = sum(delta_old .* deg_prop);
		
		FMatHat = delta_old * ones(1,m);
		FMatHat = FMatHat + FMatHat';  FMatHat = 1./(1+exp(-FMatHat));
		FMatHat = FMatHat - diag(diag(FMatHat));
		
		G = deg_freq .* (FMatHat*deg_freq) + deg_freq .* (deg_freq-1) .* 1./(1+exp(-2*delta_old)) - deg_freq .* unique_degs + deg_freq .* (delta_old - avg_delta_old) * lambda;
		
		% gradient descent
		delta_new = delta_old - stepsize * G;
		
		
		if(adaptive_stepsize == 1)
			
			% compute objective function value under old and new delta estimates
			
			avg_delta_old = sum(delta_old .* deg_prop);
			
			FM_old = delta_old*ones(1,m);
			FM_old = FM_old + FM_old';  FM_old = log(1+exp(FM_old));
			FM_old = FM_old - diag(diag(FM_old));
			like_old = deg_freq' * FM_old * deg_freq/2 + sum( deg_freq .* (deg_freq-1)/2 .* log(1+exp(2*delta_old)) ) - sum( deg_freq .* unique_degs .* delta_old ) + lambda/2 * sum( deg_freq .* (delta_old - avg_delta_old).^2 );
			
			
			avg_delta_new = sum(delta_new .* deg_prop);
			
			FM_new = delta_new*ones(1,m);
			FM_new = FM_new + FM_new';  FM_new = log(1+exp(FM_new));
			FM_new = FM_new - diag(diag(FM_new));
			like_new = deg_freq' * FM_new * deg_freq/2 + sum( deg_freq .* (deg_freq-1)/2 .* log(1+exp(2*delta_new)) ) - sum( deg_freq .* unique_degs .* delta_new ) + lambda/2 * sum( deg_freq .* (delta_new - avg_delta_new).^2 );
			
			% decide on stepsize
			if(like_new<like_old)
				% only update epsilon when delta_new is recognized
				epsilon = norm(sqrt(deg_prop).*(delta_new - delta_old)) / norm(sqrt(deg_prop).*delta_old);
				
				stepsize = stepsize*1.2;
				delta_old = delta_new;
			else
				stepsize = stepsize*0.8;
			end;
		else
			delta_old = delta_new;
		end
		
		count = count + 1;
		
	end
	
	delta_est = delta_new;
	beta_est = zeros(n,1);
	
	for(ii = 1:m)
		beta_est(dictionary{ii}) = delta_est(ii);
	end
	
end