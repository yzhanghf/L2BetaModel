rng(2608)

%%%%%%%% data 1: COVID-19 student social network
% paper: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0236337#sec025
% data:  https://osf.io/p8yrg/?show=view

DataName = "COVID-19-student-network";

var_names = ["studentindex",...
			"gender","depression","anxiety","stress","loneliness",...
			"friend","p.interaction","e.support","inf.support","co-study"];

Table_marginal_2019_09 = readmatrix('data/student_loneliness/marginal-2019-09.csv');
Table_marginal_2020_04 = readmatrix('data/student_loneliness/marginal-2020-04.csv');

Combined_table = [];
count = 1;
Common_index = [];
for(nn = 1:length(Table_marginal_2019_09))
	pos=find(Table_marginal_2020_04(:,1)==Table_marginal_2019_09(nn,1));
	if length(pos)>=1
		Combined_table(count,:,1) = Table_marginal_2019_09(nn,:);
		Combined_table(count,:,2) = Table_marginal_2020_04(pos(1),:);
		Common_index(count)       = Table_marginal_2019_09(nn,1);
		count = count + 1;
	end
end


% model estimation
lambda = 10;
epsilon_0 = 1e-3;
IterMax = 1000;

n = length(Table_marginal_2019_09);
beta_0 = zeros(n,1);
% stepsize = 1/((n-2)/2+lambda);
output_marginal_2019_09 = [];
output_marginal_2019_09(:,1:6) = Table_marginal_2019_09(:,1:6);
for(icol = 7:11)
	output_marginal_2019_09(:,icol) = ...
		Newton_method_deg(Table_marginal_2019_09(:,icol), beta_0, lambda, epsilon_0, IterMax);
		% Gradient_descend_deg(Table_2019_04(:,icol), beta_0, lambda, stepsize, epsilon_0, IterMax);
end


n = length(Table_marginal_2020_04);
beta_0 = zeros(n,1);
% stepsize = 1/((n-2)/2+lambda);
output_marginal_2020_04 = [];
output_marginal_2020_04(:,1:6) = Table_marginal_2020_04(:,1:6);
for(icol = 7:11)
	output_marginal_2020_04(:,icol) = ...
		Newton_method_deg(Table_marginal_2020_04(:,icol), beta_0, lambda, epsilon_0, IterMax);
		% Gradient_descend_deg(Table_2019_04(:,icol), beta_0, lambda, stepsize, epsilon_0, IterMax);
end

output_diff = [];
for(nn = 1:length(Common_index))
	x1 = find( (output_marginal_2019_09(:,1)==Common_index(nn)) );
	x2 = find( (output_marginal_2020_04(:,1)==Common_index(nn)) );
	output_diff(nn,3:11) = output_marginal_2020_04(x2,3:11) - output_marginal_2019_09(x1,3:11);
	output_diff(nn,1:2) = output_marginal_2020_04(nn,1:2);
end



%%% Interpretation

[~,ax] = plotmatrix(output_marginal_2019_09(:,3:11));
for(ii=1:9)
	ax(9,ii).XLabel.String = var_names(ii+2);
	ax(9,ii).XTick = [];
	ax(9,ii).XLabel.Rotation = 30;
	ax(9,ii).XLabel.Position(1) = ax(9,ii).XLabel.Position(1);
	ax(9,ii).XLabel.Position(2) = ax(9,ii).XLabel.Position(2)+1;
	ax(9,ii).XLabel.FontSize = 15;
	ax(9,ii).XLabel.HorizontalAlignment = 'right';
	ax(9,ii).XLabel.VerticalAlignment = 'top';
	
	
	ax(ii,1).YLabel.String = var_names(ii+2);
	ax(ii,1).YTick = [];
	ax(ii,1).YLabel.Rotation = 30;
	ax(ii,1).YLabel.Position(1) = ax(ii,1).YLabel.Position(1)-0;
	ax(ii,1).YLabel.Position(2) = ax(ii,1).YLabel.Position(2)+1;
	ax(ii,1).YLabel.FontSize = 15;
	ax(ii,1).YLabel.HorizontalAlignment = 'right';
	ax(ii,1).YLabel.VerticalAlignment = 'middle';
end
set(gcf,'papersize',get(gcf,'papersize')+[10,10]);
set(gcf, 'paperposition',get(gcf, 'paperposition')+[10,10,0,0]);
saveas(gcf, './data/student_loneliness/output_2019_09_var_scatter.png')


[~,ax] = plotmatrix(output_marginal_2020_04(:,3:11));
for(ii=1:9)
	ax(9,ii).XLabel.String = var_names(ii+2);
	ax(9,ii).XTick = [];
	ax(9,ii).XLabel.Rotation = 30;
	ax(9,ii).XLabel.Position(1) = ax(9,ii).XLabel.Position(1);
	ax(9,ii).XLabel.Position(2) = ax(9,ii).XLabel.Position(2)+1;
	ax(9,ii).XLabel.FontSize = 15;
	ax(9,ii).XLabel.HorizontalAlignment = 'right';
	ax(9,ii).XLabel.VerticalAlignment = 'top';
	
	
	ax(ii,1).YLabel.String = var_names(ii+2);
	ax(ii,1).YTick = [];
	ax(ii,1).YLabel.Rotation = 30;
	ax(ii,1).YLabel.Position(1) = ax(ii,1).YLabel.Position(1)-0;
	ax(ii,1).YLabel.Position(2) = ax(ii,1).YLabel.Position(2)+1;
	ax(ii,1).YLabel.FontSize = 15;
	ax(ii,1).YLabel.HorizontalAlignment = 'right';
	ax(ii,1).YLabel.VerticalAlignment = 'middle';
end
set(gcf,'papersize',get(gcf,'papersize')+[0,10]);
set(gcf, 'paperposition',get(gcf, 'paperposition')+[1,10,0,0]);
saveas(gcf, './data/student_loneliness/output_2020_04_var_scatter.png')


[~,ax] = plotmatrix(output_diff(:,3:11));
for(ii=1:9)
	ax(9,ii).XLabel.String = var_names(ii+2);
	ax(9,ii).XTick = [];
	ax(9,ii).XLabel.Rotation = 30;
	ax(9,ii).XLabel.Position(1) = ax(9,ii).XLabel.Position(1);
	ax(9,ii).XLabel.Position(2) = ax(9,ii).XLabel.Position(2)+1;
	ax(9,ii).XLabel.FontSize = 15;
	ax(9,ii).XLabel.HorizontalAlignment = 'right';
	ax(9,ii).XLabel.VerticalAlignment = 'top';
	
	
	ax(ii,1).YLabel.String = var_names(ii+2);
	ax(ii,1).YTick = [];
	ax(ii,1).YLabel.Rotation = 30;
	ax(ii,1).YLabel.Position(1) = ax(ii,1).YLabel.Position(1)-0;
	ax(ii,1).YLabel.Position(2) = ax(ii,1).YLabel.Position(2)+1;
	ax(ii,1).YLabel.FontSize = 15;
	ax(ii,1).YLabel.HorizontalAlignment = 'right';
	ax(ii,1).YLabel.VerticalAlignment = 'middle';
end
set(gcf,'papersize',get(gcf,'papersize')+[0,10]);
set(gcf, 'paperposition',get(gcf, 'paperposition')+[1,10,0,0]);
saveas(gcf, './data/student_loneliness/output_diff_var_scatter.png')


















output_marginal_2019_09(:,3:11) = normalize(output_marginal_2019_09(:,3:11),1);
output_marginal_2020_04(:,3:11) = normalize(output_marginal_2020_04(:,3:11),1);
output_diff(:,3:11)    = normalize(output_diff(:,3:11),1);


% CCA
mental_subset = 3:6;
[marginal_A_2019_09,marginal_B_2019_09] = canoncorr(output_marginal_2019_09(:,mental_subset),output_marginal_2019_09(:,7:11));
[marginal_A_2020_04,marginal_B_2020_04] = canoncorr(output_marginal_2020_04(:,mental_subset),output_marginal_2020_04(:,7:11));
[A_diff, B_diff] = canoncorr(output_diff(:,mental_subset),output_diff(:,7:11));


CCA_output = ...
	table(var_names(3:11)',...
	[marginal_A_2019_09(:,1);marginal_B_2019_09(:,1)]*sign(marginal_A_2019_09(length(mental_subset),1)),...
	[marginal_A_2020_04(:,1);marginal_B_2020_04(:,1)]*sign(marginal_A_2020_04(length(mental_subset),1)),...
	[A_diff(:,1);B_diff(:,1)]*sign(A_diff(length(mental_subset),1)));
CCA_output.Properties.VariableNames = ...
	{'Variable','2019-09','2020-04','Diff'};
CCA_output








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% permutation inference (bootstrap)
N_boot = 100000;
CCA_A_coef_boot_2019_09 = [];  CCA_B_coef_boot_2019_09 = [];
CCA_A_coef_boot_2020_04 = [];  CCA_B_coef_boot_2020_04 = [];
CCA_A_coef_boot_diff = [];  CCA_B_coef_boot_diff = [];
nnn1 = length(output_marginal_2019_09);
nnn2 = length(output_marginal_2020_04);
nnn = length(output_diff);
for(nbt = 1:N_boot)
	sfidx_2019_09 = randsample(nnn1,nnn1);
	sfidx_2020_04 = randsample(nnn2,nnn2);
	sfidx = randsample(nnn,nnn);
	[temp1,temp2] = canoncorr(output_marginal_2019_09(sfidx_2019_09,mental_subset),output_marginal_2019_09(:,7:11));
		CCA_A_coef_boot_2019_09(:,nbt) = temp1(:,1)*sign(temp1(length(mental_subset),1));
		CCA_B_coef_boot_2019_09(:,nbt) = temp2(:,1)*sign(temp1(length(mental_subset),1));
	[temp1,temp2] = canoncorr(output_marginal_2020_04(sfidx_2020_04,mental_subset),output_marginal_2020_04(:,7:11));
		CCA_A_coef_boot_2020_04(:,nbt) = temp1(:,1)*sign(temp1(length(mental_subset),1));
		CCA_B_coef_boot_2020_04(:,nbt) = temp2(:,1)*sign(temp1(length(mental_subset),1));
	[temp1,temp2] = canoncorr(output_diff(sfidx,mental_subset),output_diff(:,7:11));
		CCA_A_coef_boot_diff(:,nbt) = temp1(:,1)*sign(temp1(length(mental_subset),1));
		CCA_B_coef_boot_diff(:,nbt) = temp2(:,1)*sign(temp1(length(mental_subset),1));
end

p_val_A_2019_09 = [];  p_val_B_2019_09 = [];
p_val_A_2020_04 = [];  p_val_B_2020_04 = [];
p_val_A_diff = [];  p_val_B_diff = [];
for(ii = 1:length(mental_subset))
	p_val_A_2019_09(ii) = mean(abs(CCA_A_coef_boot_2019_09(ii,:)) > abs(marginal_A_2019_09(ii)));
	p_val_A_2020_04(ii) = mean(abs(CCA_A_coef_boot_2020_04(ii,:)) > abs(marginal_A_2020_04(ii)));
	p_val_A_diff(ii) = mean(abs(CCA_A_coef_boot_diff(ii,:)) > abs(A_diff(ii)));
end
for(jj = 1:5)
	p_val_B_2019_09(jj) = mean(abs(CCA_B_coef_boot_2019_09(jj,:)) > abs(marginal_B_2019_09(jj)));
	p_val_B_2020_04(jj) = mean(abs(CCA_B_coef_boot_2020_04(jj,:)) > abs(marginal_B_2020_04(jj)));
	p_val_B_diff(jj) = mean(abs(CCA_B_coef_boot_diff(jj,:)) > abs(B_diff(jj)));
end

fw = fopen('./data/student_loneliness/output_CCA_table.txt', 'w');
fprintf(fw, 'Variable & 2019-09 & p-val. & 2020-04 & p-val. & Diff & p-val.\\\\\\hline\n');
p_val_2019_09 = [p_val_A_2019_09,p_val_B_2019_09];
p_val_2020_04 = [p_val_A_2020_04,p_val_B_2020_04];
p_val_diff = [p_val_A_diff,p_val_B_diff];
for(irow = 1:size(CCA_output,1))
	fprintf(fw, '%s & %1.3f & (%1.3f) & %1.3f & (%1.3f) & %1.3f & (%1.3f)\\\\\n', ...
		CCA_output{irow,1},...
		CCA_output{irow,2},p_val_2019_09(irow),...
		CCA_output{irow,3},p_val_2020_04(irow),...
		CCA_output{irow,4},p_val_diff(irow)...
	);
	if(irow==1 | irow==length(mental_subset))
		fprintf(fw, '\\hline\n');
	end
end
fprintf(fw,'\\hline');
fclose(fw);


fig = figure;

extraversion_new_2019_09 = readmatrix('data/student_loneliness/extraversion_new-2019-09.csv');
extraversion_new_2020_04 = readmatrix('data/student_loneliness/extraversion_new-2020-04.csv');

s1 = scatter(extraversion_new_2019_09,output_marginal_2019_09(:,7:11)*marginal_B_2019_09(:,1),'g*');
hold on;
s2 = scatter(extraversion_new_2020_04,output_marginal_2020_04(:,7:11)*marginal_B_2020_04(:,1),'k^');
reg_line = lsline;
l1 = reg_line(1);  l2 = reg_line(2);
l1.Color = 'blue';  l1.LineStyle = '--';  l1.LineWidth = 2;
l2.Color = 'red';  l2.LineStyle = '--';  l2.LineWidth = 2;
hold off;

ylabel('Loneliness score','interpreter','latex', 'FontSize', 15);
xlabel('CCA: $\hat\beta$','interpreter','latex', 'FontSize', 15);
xlim([1.5,5]);  ylim([-2.5,2.5]);
set(gca,'fontsize',15);

corr(extraversion_new_2019_09,output_marginal_2019_09(:,7:11)*marginal_B_2019_09(:,1),'rows','complete')
corr(extraversion_new_2020_04,output_marginal_2020_04(:,7:11)*marginal_B_2020_04(:,1),'rows','complete')

saveas(fig, './data/student_loneliness/output_loneliness_vs_connection.png')







%%% SENSITIVITY ANALYSIS
% model estimation
lambda_vec = 10.^([-1:0.5:2]);

epsilon_0 = 1e-5;
IterMax = 1000;
sens_output_marginal_2019_09 = [];
sens_output_marginal_2020_04 = [];
for(lambda_index = 1:length(lambda_vec))

	lambda = lambda_vec(lambda_index);
	
	n = length(Table_marginal_2019_09);
	beta_0 = zeros(n,1);
	% stepsize = 1/((n-2)/2+lambda);
	for(icol = 7:11)
		sens_output_marginal_2019_09(:,lambda_index,icol-6) = ...
			Newton_method_deg(Table_marginal_2019_09(:,icol), beta_0, lambda, epsilon_0, IterMax);
			% Gradient_descend_deg(Table_2019_04(:,icol), beta_0, lambda, stepsize, epsilon_0, IterMax);
	end


	n = length(Table_marginal_2020_04);
	beta_0 = zeros(n,1);
	% stepsize = 1/((n-2)/2+lambda);
	for(icol = 7:11)
		sens_output_marginal_2020_04(:,lambda_index,icol-6) = ...
			Newton_method_deg(Table_marginal_2020_04(:,icol), beta_0, lambda, epsilon_0, IterMax);
			% Gradient_descend_deg(Table_2019_04(:,icol), beta_0, lambda, stepsize, epsilon_0, IterMax);
	end
	
end

colors={'k','r','b','m','g'};

fig = figure;
ppp = {};
for(kk = 1:5)
	for(ii = 1:size(sens_output_marginal_2019_09(:,:,kk),1))
		if(ii == 1)
			ppp{kk} = plot(1:length(lambda_vec), sens_output_marginal_2019_09(ii,:,kk), colors{kk}, 'DisplayName',var_names{kk});
			hold on;
		else
			plot(1:length(lambda_vec), sens_output_marginal_2019_09(ii,:,kk), colors{kk});
		end
	end
	xticks(1:length(lambda_vec)); xticklabels(round(lambda_vec,2));
	xlabel('$\lambda$', 'interpreter','latex');
	ylabel('$\widehat\beta_\lambda$', 'interpreter','latex');
	set(gca,'fontsize',15);
end
hold off;
legend([ppp{1:5}],var_names{7:11},'location','southeast');
saveas(fig, './data/student_loneliness/output_sensitivity_beta_2019_09.png')



fig = figure;
ppp = {};
for(kk = 1:5)
	for(ii = 1:size(sens_output_marginal_2020_04(:,:,kk),1))
		if(ii == 1)
			ppp{kk} = plot(1:length(lambda_vec), sens_output_marginal_2020_04(ii,:,kk), colors{kk}, 'DisplayName',var_names{kk});
			hold on;
		else
			plot(1:length(lambda_vec), sens_output_marginal_2020_04(ii,:,kk), colors{kk});
		end
	end
	xticks(1:length(lambda_vec)); xticklabels(round(lambda_vec,2));
	xlabel('$\lambda$', 'interpreter','latex');
	ylabel('$\widehat\beta_\lambda$', 'interpreter','latex');
	set(gca,'fontsize',15);
end
hold off;
legend([ppp{1:5}],var_names{7:11},'location','southeast');
saveas(fig, './data/student_loneliness/output_sensitivity_beta_2020_04.png')






