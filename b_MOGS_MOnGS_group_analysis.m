clear; clc;
%% Analysis 1: Effects of each group on micro-offline and online gains for Mogs_Hipp
%
% Requirements: Access to files generated from "a_MOGS_MOnGS_calc.m".
%
% The paths to these folders must be specified here:

out_dir = 'G:\My Drive\DataScience\UT_Austin\MOGS\processed_data'; % where your data output from a_MOGS_MOnGS_calc.m went
bp = 11; % Set break point (Bonstrup et al. 2019 used 11).

% Notes: 
% 1) The first analysis determined whether the group effect depended on
%    break point. Twenty-four one-way ANOVAs were conducted (one for each
%    break point). The graph shows the p-value for all 24 ANOVAs.
% 2) The second analysis is a one-way ANOVA using MOGS on all trials. A
%    violin plot is generated showing the average MOGS for all three groups.
% 3) The last analysis direclty compares MOGS between the REST and ENC
%    groups.
%
%    Output:
%    bp_off.mat --> bp_off (av. MOGS calculated with all inter-trial periods)
%

%% Load in data
load([out_dir, '/Microoffline_data.mat'])
load([out_dir, '/Microonline_data.mat'])
load([out_dir, '/ENC_SEM_Memory_Performance.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('         Started MOGS group Analysis         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');

%% Contrasting MOGS and MOnGs against zero

% MOGS
disp('Testing MOGS against zero...')
zrs = zeros(15,1);
rest_mogs = mean(microoffline(1:11,:,1),1,'omitmissing')';
disp(['Mean MOGS = ', num2str(mean(rest_mogs(1:bp))), ', SEM = ', num2str(std(rest_mogs(1:bp))/sqrt(length(rest_mogs))), '.'])
[~,p,~,stats] = ttest(rest_mogs, zrs,'Alpha',0.05,'Tail','right');

if p < 0.05
    disp('Reject the null hypothesis that MOGS is equal to zero in the REST group');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
else
    disp('Fail to reject the null hypothesis that MOGS is equal to zero in the REST group');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
end

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

% MOnGS
disp('Testing MOnGS against zero...')
zrs = zeros(15,1);
rest_mons = mean(microonline(1:11,:,1),1,'omitmissing')';
disp(['Mean MOnGS = ', num2str(mean(rest_mons(1:bp))), ', SEM = ', num2str(std(rest_mons(1:bp))/sqrt(length(rest_mons))), '.'])
[~,p,~,stats] = ttest(rest_mons, zrs,'Alpha',0.05,'Tail','left');

if p < 0.05
    disp('Reject the null hypothesis that MOnGS is equal to zero in the REST group');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
else
    disp('Fail to reject the null hypothesis that MOnGS is equal to zero in the REST group');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
end

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

% MOGS vs. MOnGS
disp('Testing MOGS against MOnGS...')
[~,p,~,stats] = ttest(rest_mogs, rest_mons,'Alpha',0.05,'Tail','right');

if p < 0.05
    disp('Reject the null hypothesis that MOGS and MOnGS are equal in the REST group');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
else
    disp('Fail to reject the null hypothesis that MOGS and MOnGS are equal to zero in the REST group');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
end

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

%% Contrasting associative memory performance against zero for the ENC group
disp('Testing associative memory performance against zero for the ENC group...')
ENC_mem = mean(mem_perf_sorted(:,:,1),1,'omitmissing');
zrs = zeros(15,1)';
disp(['Mean encoding success = ', num2str(mean(ENC_mem(1:bp))), ', SEM = ', num2str(std(ENC_mem(1:bp))/sqrt(length(ENC_mem))), '.'])
[~,p,~,stats] = ttest(ENC_mem, zrs,'Alpha',0.05,'Tail','right');

if p < 0.05
    disp('Reject the null hypothesis that memory performance equals zero');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
else
    disp('Fail to reject the null hypothesis that memory performance equals zero');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
end

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

%% Contrasting episodic memory performance against zero for the SEM group
disp('Testing Episodic memory performance against zero for the SEM group...')
SEM_mem = mean(mem_perf_sorted(:,:,2),1,'omitmissing')';
zrs = zeros(15,1);
disp(['Mean encoding success = ', num2str(mean(SEM_mem(1:bp))), ', SEM = ', num2str(std(SEM_mem(1:bp))/sqrt(length(SEM_mem))), '.'])
[~,p,~,stats] = ttest(SEM_mem, zrs,'Alpha',0.05,'Tail','right');

if p < 0.05
    disp('Reject the null hypothesis that memory performance equals zero');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
else
    disp('Fail to reject the null hypothesis that memory performance equals zero');
    disp(['t(', num2str(stats.df), ') = ' num2str(stats.tstat), ', p = ' num2str(p) '.']);
end

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));
%% Contrasting av. MOGS between groups using all blocks

% Organize data
bp_off = zeros(15,3);
for grp=1:3
    for subject=1:15    
        bp_off(subject, grp) = mean(microoffline(1:bp,subject,grp),'omitmissing');
    end
end

% Plot average MOGS for all groups
a = violinplot(bp_off); 
clr = {'black','blue','red'};

for i=1:3
  a(i).ViolinColor=clr{i}; a(i).ScatterPlot.SizeData=104; a(i).ScatterPlot.MarkerFaceAlpha=0.9; 
  a(i).ScatterPlot.MarkerFaceColor=clr{i}; a(i).ScatterPlot.MarkerEdgeColor='black'; a(i).EdgeColor='black'; a(i).BoxColor='black';
  a(i).ShowMean=1; a(i).MedianPlot.SizeData=150; a(i).MedianPlot.MarkerEdgeColor=[0,0,0]; a(i).MedianPlot.LineWidth=2.5;
  a(i).MeanPlot.LineWidth=2.5; a(i).MeanPlot.Color=[0,0,0];
end

ylim([-8 8]);
yline(0, 'k--', 'LineWidth', 3)
set(gca,'FontSize',24);
xticklabels({'REST','ENC','SEM'}); 
title('Average MOGS by group','Fontsize',54);
xlabel('Group','FontSize',36)
ylabel('MOGS','FontSize',36)
set(gca,'FontSize',42,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
  
% Contrast MOGS by group using ANOVA
aov_tmp_off = anova(bp_off); 
pval_tmp_off = aov_tmp_off.stats{1,5}; 
disp(['One-way ANOVA comparing MOGS between groups : F(' num2str(aov_tmp_off.stats{1,2}) ',' num2str(aov_tmp_off.stats{2,2}) ') = ' num2str(aov_tmp_off.stats{1,4}) ', p = ' num2str(pval_tmp_off) '.']); 
[~,p,~,stats] = ttest2(bp_off(:,1), bp_off(:,2), "Tail","both","Vartype","unequal"); tp_val = p; disp(['T-test between rest and ENC: t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(tp_val) '.']);
[~,p,~,stats] = ttest2(bp_off(:,3), bp_off(:,2), "Tail","both","Vartype","unequal"); tp_val = p; disp(['T-test between SEM and ENC: t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(tp_val) '.']);
disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

clearvars -except bp_off microonline microoffline aov_tmp_off out_dir bp 

%% Contrasting av. MOnGS between groups using all blocks

% Organize data
bp_on = zeros(15,3);
for grp=1:3
    for subject=1:15    
        bp_on(subject, grp) = mean(microonline(1:bp,subject,grp),'omitmissing');
    end
end

% Plot average MOGS for all groups
a = violinplot(bp_on); 
clr = {'black','blue','red'};

for i=1:3
  a(i).ViolinColor=clr{i}; a(i).ScatterPlot.SizeData=104; a(i).ScatterPlot.MarkerFaceAlpha=0.9; 
  a(i).ScatterPlot.MarkerFaceColor=clr{i}; a(i).ScatterPlot.MarkerEdgeColor='black'; a(i).EdgeColor='black'; a(i).BoxColor='black';
  a(i).ShowMean=1; a(i).MedianPlot.SizeData=150; a(i).MedianPlot.MarkerEdgeColor=[0,0,0]; a(i).MedianPlot.LineWidth=2.5;
  a(i).MeanPlot.LineWidth=2.5; a(i).MeanPlot.Color=[0,0,0];
end

ylim([-8 8]);
yline(0, 'k--', 'LineWidth', 3)
set(gca,'FontSize',24);
xticklabels({'REST','ENC','SEM'}); 
title('Average MOnGS by group','Fontsize',54);
xlabel('Group','FontSize',36)
ylabel('MOnGS','FontSize',36)
set(gca,'FontSize',42,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
  
% Contrast MOnGS by group using ANOVA
aov_tmp_on = anova(bp_on); 
pval_tmp_on = aov_tmp_on.stats{1,5}; 
disp(['One-way ANOVA comparing MOnGS between groups : F(' num2str(aov_tmp_on.stats{1,2}) ',' num2str(aov_tmp_on.stats{2,2}) ') = ' num2str(aov_tmp_on.stats{1,4}) ', p = ' num2str(pval_tmp_on) '.']); 
[~,p,~,stats] = ttest2(bp_on(:,1), bp_on(:,2), "Tail","both","Vartype","unequal"); tp_val = p; disp(['T-test between rest and ENC: t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(tp_val) '.']);
[~,p,~,stats] = ttest2(bp_on(:,3), bp_on(:,2), "Tail","both","Vartype","unequal"); tp_val = p; disp(['T-test between SEM and ENC: t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(tp_val) '.']);
disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

clearvars -except bp_off microonline microoffline aov_tmp_off aov_tmp_on out_dir bp bp_on

%% Breakpoint analysis

% Breakpoint analysis (MOGS/MONS)
bp_groups = zeros(15, 24, 3);
for grp=1:3
    for subject = 1:15
        for bp = 1:24 % all possible breakpoints
            bp_groups(subject,bp,grp)= mean(microoffline(1:bp-1,subject,grp),1, 'omitmissing');
        end
    end
end
bp_groups(:,1,:) = [];

% Make data structure with MOGS for all break points.
bp_groups(10,1,:) = NaN;
bp_groups_m = zeros(23,3);
bp_groups_s = zeros(23,3);
for i=1:23 % Block
    for j=1:3 % Group
        bp_groups_m(i,j) = mean(bp_groups(:,i,j), 'omitmissing');
        bp_groups_s(i,j) = std(bp_groups(:,i,j), 'omitmissing')/(sqrt(length(bp_groups(:,i,j))));
    end
end

% Get the group effect p-value for all 24 break point ANOVAs 
ps = zeros(23,1);
disp(' --- ANOVA on each breakpoint (MOGS) --- ');
for bp = 1:23
    tmp_bp_mogs = [bp_groups(:,bp,1) bp_groups(:,bp,2) bp_groups(:,bp,3)]; 
    aov_bp_mogs = anova(tmp_bp_mogs); pval_bpmogs = aov_bp_mogs.stats{1,5}; 
    disp(['One-way ANOVA comparing MOGS bp = ' num2str(bp+1) ': F(' num2str(aov_bp_mogs.stats{1,2}) ',' num2str(aov_bp_mogs.stats{2,2}) ') = ' num2str(aov_bp_mogs.stats{1,4}) ', p = ' num2str(pval_bpmogs) '.']);
    clear tmp_bp_mogs
    ps(bp,1) = pval_bpmogs; clear pval_bpmogs;
end

% Plot the results
yyaxis left;
pplot = plot(ps); 
pplot.LineWidth = 5; 
pplot.Color = 'k'; 
xlim([1 23]); 
set(gca,'FontSize',24);
ylim([0 1]); 
xlabel('Break point','Fontsize',36)
ylabel(['p-value' newline],'Fontsize',36)
yline(0.05, 'k--', 'LineWidth', 3);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
sgt = sgtitle(['P-values for Group effect by break point' newline],'Color','black');
yticks([0 0.2 0.4 0.6 0.8 1])
sgt.FontSize = 54;
set(gca, 'box', 'off')
ax=gca;
ax.YAxis(1).Color = 'k';
hold on

% MOGS across groups
yyaxis right;
ax.YAxis(2).Color = 'k';
x = 1:23;
y = bp_groups_m';
err = bp_groups_s';
shadedErrorBar(x(1,:),y(1,:),err(1,:),'lineprops',{'-ko','MarkerFaceColor','k','MarkerSize',12});
shadedErrorBar(x(1,:),y(2,:),err(2,:),'lineprops',{'-bo','MarkerFaceColor','b','MarkerSize',12});
shadedErrorBar(x(1,:),y(3,:),err(3,:),'lineprops',{'-ro','MarkerFaceColor','r','MarkerSize',12});
ylim([-4 3]);
xlim([1 23]);
xticks([5, 10, 15, 20]);
xticklabels({'','','',''})
ylabel('MOGS','FontSize',36)
set(get(gca,'YLabel'),'Rotation',270); % Rotate it
hold off
disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

clearvars -except bp_off microonline microoffline aov_bp_mogs aov_tmp out_dir bp bp_on bp_microoffline bp_groups


%% Method 2 for determining whether break point matters.
bp_data = pagetranspose(bp_groups(:,1:23,:));
bp_data_org = [bp_data(:,:,1)';bp_data(:,:,2)';bp_data(:,:,3)'];

% Create single vector out of bp_data_org
bp_data_org = bp_data_org(:);
disp(' --- Trial by Group ANOVA ---')

% Create independent variable for ANOVAn for Block
trl_covar = 2:24;
blk_iv = repmat(trl_covar, 1, 45)';

% Create independent variable for ANOVAn for GROUP
Groups = [1,2,3];
Group_iv = repelem(Groups, 1, 345)'; % Repeats each element 3 times

% Create subject variable for ANOVAn for GROUP
subs = 1:45;
subs_iv = repelem(subs,1, 23)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl = [Group_iv, subs_iv, blk_iv, bp_data_org]; 

% Prepare the table for LME analysis
tbl = array2table(tbl, 'VariableNames', {'Group', 'Subject', 'Block', 'Microoffline'});

% Run LME Contrasting BLOCK and GROUP with a BLOCK COVARIATE
lme = fitlme(tbl,'Microoffline~Group+Block+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME:');
disp(lme);
disp(['p-value for interaction: ' num2str(round(lme.Coefficients{2,6},3))]);

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));
clearvars -except bp_off bp_on microonline microoffline aov_bp_mogs aov_tmp out_dir

%% Linear mixed effects model at trial level (MOGS)

% Only keep first 10 MOGS
microoffline10 = microoffline(1:10,:, :); % Keep only the first 10 MOGS

% take microoffline and make into into a single column of data
trial_mogs = microoffline10(:); % Reshape microoffline into a single column

% Create independent variable for ANOVAn for Block
trl_covar = 1:10;
blk_iv = repmat(trl_covar, 1, 45)';

% Create independent variable for ANOVAn for GROUP
Groups = [1,2,3];
Group_iv = repelem(Groups, 1, 150)'; % Repeats each element 3 times

% Create subject variable for ANOVAn for GROUP
subs = 1:45;
subs_iv = repelem(subs,1, 10)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl = [Group_iv subs_iv blk_iv trial_mogs]; 

% Remove subject 40 from tbl
tbl(subs_iv == 40,:) = [];

% Prepare the table for LME analysis
tbl = array2table(tbl, 'VariableNames', {'Group', 'Subject', 'Block', 'Microoffline'});

% Run LME Contrasting BLOCK and GROUP 
lme = fitlme(tbl,'Microoffline~Group*Block+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME:');
disp(lme);
disp(['p-value for interaction: ' num2str(round(lme.Coefficients{4,6},3))]);


disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

% Plot results
plot_mogs_m = zeros(10,3);
plot_mogs_s = zeros(10,3);
for i=1:3
    if i == 3
        % Remove column 10
        tmp_sem = [microoffline10(:,1:9,3), microoffline10(:,11:15,3)];
        plot_mogs_m(:,3) = mean(tmp_sem,2,'omitmissing');
        plot_mogs_s(:,3) = std(tmp_sem,0,2,'omitmissing')/sqrt(14);
    else
        plot_mogs_m(:,i) = mean(microoffline10(:,:,i),2,'omitmissing');
        plot_mogs_s(:,i) = std(microoffline10(:,:,i),0,2,'omitmissing')/sqrt(15);
    end
end

x = 1:10;
y = plot_mogs_m';
err = plot_mogs_s';
shadedErrorBar(x(1,:),y(1,:),err(1,:),'lineprops',{'-ko','MarkerFaceColor','k','MarkerSize',12});
hold on
shadedErrorBar(x(1,:),y(2,:),err(2,:),'lineprops',{'-bo','MarkerFaceColor','b','MarkerSize',12});
shadedErrorBar(x(1,:),y(3,:),err(3,:),'lineprops',{'-ro','MarkerFaceColor','r','MarkerSize',12});
ylim([-4 8]);
xlim([1 10]);
xticks([2,4,6,8,10]);
xticklabels({'','','','','',''})
yticks([-4, -2,0,2,4,6,8])
yticklabels({'','','','','','',''})
ylabel('','FontSize',36)
box off
hold off
disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

%% Linear mixed effects model at trial level (MOnGS)
% Only keep first 11 MOnGS
microonline10 = microonline(1:11,:, :); % Keep only the first 11 MOnGS

% take microonline and make into into a single column of data
trial_mongs = microonline10(:); % Reshape microonline into a single column

% Create independent variable for ANOVAn for Block
trl_covar = 1:11;
blk_iv = repmat(trl_covar, 1, 45)';

% Create independent variable for ANOVAn for GROUP
Groups = [1,2,3];
Group_iv = repelem(Groups, 1, 165)'; % Repeats each element 3 times

% Create subject variable for ANOVAn for GROUP
subs = 1:45;
subs_iv = repelem(subs,1, 11)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl = [Group_iv subs_iv blk_iv trial_mongs]; 

% Remove subject 40 from tbl
tbl(subs_iv == 40,:) = [];

% Prepare the table for LME analysis
tbl = array2table(tbl, 'VariableNames', {'Group', 'Subject', 'Block', 'Microonline'});

% Run LME Contrasting BLOCK and GROUP 
lme = fitlme(tbl,'Microonline~Group*Block+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME:');
disp(lme);
disp(['p-value for interaction: ' num2str(round(lme.Coefficients{4,6},3))]);

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

% Break up interaction by Group

% --- REST ---
microonline10_rest = microonline(1:11,:, 1); % Keep only the first 11 MOnGS

% take microonline and make into into a single column of data
trial_mongs_rest = microonline10_rest(:); % Reshape microonline into a single column

% Create independent variable for ANOVAn for Block
blk_iv_rest = repmat(trl_covar, 1, 15)';

% Create subject variable for ANOVAn for GROUP
subs = 1:15;
subs_iv_rest = repelem(subs,1, 11)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl_rest = [subs_iv_rest blk_iv_rest trial_mongs_rest]; 

% Prepare the table for LME analysis
tbl_rest = array2table(tbl_rest, 'VariableNames', {'Subject', 'Block', 'Microonline'});

% Run LME Contrasting BLOCK and GROUP 
lme_rest = fitlme(tbl_rest,'Microonline~Block+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME:');
disp(lme);
disp(['p-value for trial effect: ' num2str(round(lme_rest.Coefficients{2,6},3))]);
disp([' --- No trial effect for REST group --- '])

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

% --- ENC ---
microonline10_enc = microonline(1:11,:, 2); % Keep only the first 11 MOnGS

% take microonline and make into into a single column of data
trial_mongs_enc = microonline10_enc(:); % Reshape microonline into a single column

% Create independent variable for ANOVAn for Block
blk_iv_enc = repmat(trl_covar, 1, 15)';

% Create subject variable for ANOVAn for GROUP
subs = 1:15;
subs_iv_enc = repelem(subs,1, 11)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl_enc = [subs_iv_enc blk_iv_enc trial_mongs_enc]; 

% Prepare the table for LME analysis
tbl_enc = array2table(tbl_enc, 'VariableNames', {'Subject', 'Block', 'Microonline'});

% Run LME Contrasting BLOCK and GROUP 
lme_enc = fitlme(tbl_enc,'Microonline~Block+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME:');
disp(lme);
disp(['p-value for trial effect: ' num2str(round(lme_enc.Coefficients{2,6},3))]);
disp([' --- P < 0.05 for trial effect ENC group --- '])

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

% --- SEM (PICK UP HERE!) ---
microonline10_sem = microonline(1:11,:, 3); % Keep only the first 11 MOnGS

% Remove subject 10
microonline10_sem(:,10) = []; % Remove subject 10

% take microonline and make into into a single column of data
trial_mongs_sem = microonline10_sem(:); % Reshape microonline into a single column

% Create independent variable for ANOVAn for Block
blk_iv_sem = repmat(trl_covar, 1, 14)';

% Create subject variable for ANOVAn for GROUP
subs = 1:14;
subs_iv_sem = repelem(subs,1, 11)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl_sem = [subs_iv_sem blk_iv_sem trial_mongs_sem]; 

% Prepare the table for LME analysis
tbl_sem = array2table(tbl_sem, 'VariableNames', {'Subject', 'Block', 'Microonline'});

% Run LME Contrasting BLOCK and GROUP 
lme_sem = fitlme(tbl_sem,'Microonline~Block+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME:');
disp(lme);
disp(['p-value for trial effect: ' num2str(round(lme_sem.Coefficients{2,6},3))]);
disp([' --- P < 0.05 for trial effect ENC group --- '])

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

% Plot results
plot_mongs_m = zeros(11,3);
plot_mongs_s = zeros(11,3);
for i=1:3
    if i == 3
        % Remove column 10
        tmp_sem = [microonline10(:,1:9,3), microonline10(:,11:15,3)];
        plot_mongs_m(:,3) = mean(tmp_sem,2,'omitmissing');
        plot_mongs_s(:,3) = std(tmp_sem,0,2,'omitmissing')/sqrt(14);
    else
        plot_mongs_m(:,i) = mean(microonline10(:,:,i),2,'omitmissing');
        plot_mongs_s(:,i) = std(microonline10(:,:,i),0,2,'omitmissing')/sqrt(15);
    end
end

x = 1:11;
y = plot_mongs_m';
err = plot_mongs_s';
shadedErrorBar(x(1,:),y(1,:),err(1,:),'lineprops',{'-ko','MarkerFaceColor','k','MarkerSize',12});
hold on
shadedErrorBar(x(1,:),y(2,:),err(2,:),'lineprops',{'-bo','MarkerFaceColor','b','MarkerSize',12});
shadedErrorBar(x(1,:),y(3,:),err(3,:),'lineprops',{'-ro','MarkerFaceColor','r','MarkerSize',12});
ylim([-8 8]);
xlim([1 10]);
xticks([2,4,6,8,10]);
xticklabels({'','','','','',''})
yticks([-8,-6,-4, -2,0,2,4,6,8])
yticklabels({'','','','','','','','',''})
ylabel('','FontSize',36)
box off
hold off
disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));
%% Finish up 
save([out_dir, '/bp_off.mat'],'bp_off','-v7.3');
save([out_dir, '/bp_on.mat'],'bp_on','-v7.3');

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('        Finished MOGS group Analysis         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');

