clear; clc;
%% Analysis 4: Effects of each group on micro-offline and online gains for Mogs_Hipp
%
% Requirements: Access to "_processed_data" files generated from "a_MOGS_calc.m".
%
% The paths to these folders must be specified here:

home_dir = '/Users/mvf327/Library/CloudStorage/Box-Box/UT_Austin/Personal_Lab_Files/Papers/MOGS_Hipp/_analysis'; % root folder
load_dir = [home_dir, '/_processed_data'];
out_dir = [home_dir '/_processed_data'];

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
%    bp_24_off.mat --> bp_24_off (av. MOGS calculated with all inter-trial periods)
%

%% Load in data
load([load_dir, '/Microoffline_data.mat'])
load([load_dir, '/Microonline_data.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('         Started MOGS group Analysis         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');


%% Breakpoint analysis

% Breakpoint analysis (MOGS/MONS)
bp_microoffline = zeros(15, 24, 3);
for grp=1:3
    for subject = 1:15
        for bp = 1:24 % all possible breakpoints
            bp_microoffline(subject,bp,grp)= mean(microoffline(1:bp-1,subject,grp),1, 'omitmissing');
        end
    end
end

% Get the group effect p-value for all 24 break point ANOVAs 
ps = zeros(24,1);
disp(' --- ANOVA on each breakpoint (MOGS) --- ');
for bp = 1:23
    tmp_bp_mogs = [bp_microoffline(:,bp+1,1) bp_microoffline(:,bp+1,2) bp_microoffline(:,bp+1,3)]; 
    aov_bp_mogs = anova(tmp_bp_mogs); pval_bpmogs = aov_bp_mogs.stats{1,5}; 
    disp(['One-way ANOVA comparing MOGS bp = ' num2str(bp+1) ': F(' num2str(aov_bp_mogs.stats{1,2}) ',' num2str(aov_bp_mogs.stats{2,2}) ') = ' num2str(aov_bp_mogs.stats{1,4}) ', p = ' num2str(pval_bpmogs) '.']);
    clear tmp_bp_mogs
    ps(bp,1) = pval_bpmogs; clear pval_bpmogs;
end

% Plot the results
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
sgt.FontSize = 54;
set(gca, 'box', 'off')
disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

%% Contrasting av. MOGS between groups using all blocks

% Organize data
bp_24_off = zeros(15,3);
for grp=1:3
    for subject=1:15    
        bp_24_off(subject, grp) = mean(microoffline(1:24,subject,grp),'omitmissing');
    end
end

% Plot average MOGS for all groups
a = violinplot(bp_24_off); 
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
aov_tmp = anova(bp_24_off); 
pval_tmp = aov_tmp.stats{1,5}; 
disp(['One-way ANOVA comparing MOGS between groups : F(' num2str(aov_tmp.stats{1,2}) ',' num2str(aov_tmp.stats{2,2}) ') = ' num2str(aov_tmp.stats{1,4}) ', p = ' num2str(pval_tmp) '.']); 
[h,p,ci,stats] = ttest2(bp_24_off(:,1), bp_24_off(:,2)); tp_val = p; disp(['T-test between rest and ENC: t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(tp_val) '.']);
disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

clearvars -except bp_24_off microonline microoffline aov_bp_mogs aov_tmp out_dir


%% Finish up 
save([out_dir, '/bp_24_off.mat'],'bp_24_off','-v7.3');

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('        Finished MOGS group Analysis         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');

