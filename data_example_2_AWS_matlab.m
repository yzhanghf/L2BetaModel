

%%%%%%%% data 2: Amazon COVID-19 Knowledge graph
% Blog: https://aws.amazon.com/cn/blogs/database/building-and-querying-the-aws-covid-19-knowledge-graph/
% Data: https://dj2taa9i652rf.cloudfront.net/

DataName = 'AWS-COVID-19-KG';

deg_vec = readmatrix('data/AWS-csv/deg_vec.txt');

n = length(deg_vec);

% beta_0 = log(deg_vec+1)-mean(log(deg_vec+1));
beta_0 = zeros(size(deg_vec));
% beta_est = Fast_gradient_deg(deg_vec, beta_0, 10, 1/n, 1e-10, 2000, 1);

% plot beta_lambda track over lambda values
beta_track = [];
lambda_vec = [1,2.5,5,10,20,40,80,160,320];
beta_history = {};
iter_actual = [];
for(lambda_idx = 1:length(lambda_vec))
	[beta_track(:,lambda_idx), beta_history{lambda_idx}, ~, iter_actual(lambda_idx)]= ...
		Fast_Newton_method_deg(deg_vec, beta_0, lambda_vec(lambda_idx), 1e-10, 2000);
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
hold on;
yline( log(n),'b--','LineWidth',2);
yline(-log(n),'b--','LineWidth',2);
hold off;
xticks(1:length(lambda_vec)); xticklabels(lambda_vec);
xlabel('$\lambda$', 'interpreter','latex');
xlim([1,length(lambda_vec)])
ylabel('$\mathrm{Sign}(\widehat\beta_\lambda)\circ\log(|\widehat\beta_\lambda|)$', 'interpreter','latex');
set(gca,'fontsize',15);


saveas(fig, './data/AWS-csv/AWS_beta_sensitivity.png')



beta_est = Fast_Newton_method_deg(deg_vec, beta_0, 40, 1e-10, 2000);

writematrix(beta_est, './data/AWS-csv/AWS_beta_est.txt', 'Delimiter',',');

hist(beta_est)

% country_tab = readmatrix('data/AWS-csv/nation.csv');





concept_tab = readtable('data/AWS-csv/AWS_concept_rank.csv');

writetable(head(concept_tab,50), 'data/AWS-csv/AWS_concept_score_rank_top_50.csv');
writetable(flip(tail(concept_tab,50)), 'data/AWS-csv/AWS_concept_score_rank_bottom_50.csv');

take = 30;

all_concepts = unique(concept_tab.entity_String);
all_concepts2 = {'brand', ...
				'diag.',...
				'generic',...
				'procedure',...
				'organ',...
				'test',...
				'treatment'};
for idx = 1:length(all_concepts)
	xxx = all_concepts{idx};
	yyy = all_concepts2{idx};
	subtab = concept_tab(find(strcmp(concept_tab.entity_String,xxx)),:);
	subtab_head = head(subtab,take);
	subtab_tail = flip(tail(subtab,take));
	writetable(subtab_head, sprintf('data/AWS-csv/AWS_%s_concept_score_rank_top_%d.csv',yyy,take));
	writetable(subtab_tail, sprintf('data/AWS-csv/AWS_%s_concept_score_rank_bottom_%d.csv',yyy,take));
	
	fp = fopen(sprintf('data/AWS-csv/AWS_%s_concept_score_rank_top_bottom_%d.tex',yyy,take), 'w');
	% fprintf(fp, '\\multicolumn{4}{c|}{Top %d}           & \\multicolumn{4}{c}{Bottom %d} \\\\\\hline\n	\\multicolumn{1}{c}{Concept}  & \\multicolumn{1}{c}{Paper \\#} & \\multicolumn{1}{c}{Avg($d_i$)} & \\multicolumn{1}{c|}{WABS} & \\multicolumn{1}{c}{Concept}  & \\multicolumn{1}{c}{Paper \\#} & \\multicolumn{1}{c}{Avg($d_i$)} & \\multicolumn{1}{c}{WABS} \\\\\\hline\n',take,take);
	for(jjj = 1:take)
		fprintf(fp, '%s & %d & %1.2f & %1.2f & %s & %d & %1.2f & %1.2f\\\\\n', ...
				subtab_head.concept_String{jjj}, subtab_head.paper_number(jjj), subtab_head.paper_avg_deg(jjj), subtab_head.score(jjj),...
				subtab_tail.concept_String{jjj}, subtab_tail.paper_number(jjj), subtab_tail.paper_avg_deg(jjj), subtab_tail.score(jjj));
	end
	fprintf(fp, '\\\\\\hline\n');
	fclose(fp);
	
end





