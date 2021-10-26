

%%%%%%%% data 2: Kaggle Kg_nt network
% http://snap.stanford.edu/data/facebook-large-page-page-network.html

DataName = "Kaggle_kg_nt";


df1 = readtable('./data/kg_Table.txt');
in_cord19 = readmatrix('./data/kg_in_cord19.txt');
d_vec = df1.degree;
d_dict = df1.DOI;
nation_label = [df1.Other,df1.China,df1.USA,df1.UK,df1.EU,df1.JpKr,df1.India];

n = length(d_vec);

% use a hot start
beta_0 = log(d_vec+1)-mean(log(d_vec+1));


% plot beta_lambda track over lambda values
beta_track = [];
lambda_vec = [25,50,100,200,400,800];
for(lambda_idx = 1:length(lambda_vec))
	beta_track(:,lambda_idx) = ...
		Fast_Newton_method_deg(d_vec, beta_0, lambda_vec(lambda_idx), 1e-7, 500);
end

delta_track = [];
uni1 = unique(beta_track(:,length(lambda_vec)));
for(idx = 1:length(uni1))
	delta_track(idx,:) = mean(  beta_track( abs(beta_track(:,length(lambda_vec))-uni1(idx))<1e-9, :), 1  );
end

fig = figure;
plot(1:length(lambda_vec),...
	sign(delta_track).*log(abs(delta_track)),...
	'k-');
xticks(1:length(lambda_vec)); xticklabels(lambda_vec);
xlabel('$\lambda$', 'interpreter','latex');
ylabel('$\mathrm{Sign}(\widehat\beta_\lambda)\circ\log(|\widehat\beta_\lambda|)$', 'interpreter','latex');
set(gca,'fontsize',15);

saveas(fig, './data/kg_beta_sensitivity.png')


% beta_est = Fast_gradient_deg(d_vec, beta_0, 10, 1/n, 1e-7, 500, 1);
beta_est = Fast_Newton_method_deg(d_vec, beta_0, 100, 1e-7, 500);


% ONLY KEEP THOSE IN CORD-19!
beta_est = beta_est(in_cord19>0.5);
nation_label = nation_label(in_cord19>0.5,:);


writematrix(beta_est, './data/kg_beta_est.txt', 'Delimiter',',');

[sort_beta_est, order_beta_est] = sort(beta_est,'descend');
d_dict(order_beta_est(1:30))


nation_label = nation_label > 0;


mean_list = [];
std_list  = [];
entry_nos = [];
for(ii = 1:7)
	temp = beta_est(nation_label(:,ii)==1);
	entry_nos(ii) = length(temp);
	mean_list(ii) = mean(temp);
	std_list(ii)  = std(temp);
end


xaxislimit=[-7,-4];


NameList = {'Other','China','US','UK','EU','JpKr','India','beta_est'};

output_ptr = 1;
output_ptr = fopen('./data/output_kg_beta_est_by_region.txt','w')


fprintf(output_ptr, 'Region & %s & %s & %s & %s & %s & %s & %s \\\\\\hline\n', ...
		NameList{1:7});
fprintf(output_ptr, 'Entry count & %d & %d & %d & %d & %d & %d & %d\\\\\n', entry_nos);
fprintf(output_ptr, 'mean($\\hat\\beta_\\lambda$) & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f \\\\\n', mean_list);
fprintf(output_ptr, 'std($\\hat\\beta_\\lambda$) & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) \\\\\\hline\n', std_list);



fclose(output_ptr);





fig = figure;


h = {};
plot_order = [3,5,2,4,6,7];
color_list = {'k','r', 'c', 'm', 'b', 'g', 'y'};
alpha_list = [0.7, 0.7, 0.9, 0.5, 0.3, 0.7, 0.9];
for (h_idx = plot_order)
	plot_vec = beta_est(nation_label(:,h_idx)==1);
	plot_vec = plot_vec(xaxislimit(1)<plot_vec & xaxislimit(2)>plot_vec);
	h{h_idx} = histogram(plot_vec,'FaceAlpha',alpha_list(h_idx));
	xlim(xaxislimit);
	h{h_idx}.FaceColor = color_list{h_idx};
	hold on;
end
hold off;


legend(NameList{plot_order});

title({'Overlapping histograms of $\hat\beta_\lambda$','Y-axis is count'},'Interpreter','LaTeX')
xlabel('Estimated $\hat\beta_\lambda$','Interpreter','LaTeX')
ylabel('Count')
set(gca, 'FontSize', 15)

saveas(fig, './data/kg_beta_overlap_hist.png')


fig = figure;
plot_vec = beta_est(nation_label(:,1)==1);
plot_vec = plot_vec(xaxislimit(1)<plot_vec & xaxislimit(2)>plot_vec);
h00 = histogram(plot_vec);
h00.FaceColor = color_list{1};
legend('Other regions');
title({'Overlapping histograms of $\hat\beta_\lambda$','Y-axis is count'},'Interpreter','LaTeX')
xlabel('Estimated $\hat\beta_\lambda$','Interpreter','LaTeX')
ylabel('Count')
set(gca, 'FontSize', 15)

saveas(fig, './data/kg_beta_overlap_hist_others.png')






output_ptr = fopen('./data/output_kg_beta_est_by_region.txt','w')


fprintf(output_ptr, 'Region & %s & %s & %s & %s & %s & %s & %s \\\\\\hline\n', ...
		NameList{1:7});
fprintf(output_ptr, 'Entry count & %d & %d & %d & %d & %d & %d & %d\\\\\n', entry_nos);
fprintf(output_ptr, 'mean($\\hat\\beta_\\lambda$) & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f & %1.3f \\\\\n', mean_list);
fprintf(output_ptr, 'std($\\hat\\beta_\\lambda$) & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) & (%1.3f) \\\\\\hline\n', std_list);

fclose(output_ptr);

