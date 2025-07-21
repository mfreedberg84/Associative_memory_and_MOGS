clear; clc;
%% Analysis 7: Effects of group on performance and retention 
%
% Requirements: Access files generated from "a_MOGS_MOnGS_calc.m".
%
% indicate where your data output from a_MOGS_MOnGS_calc.m went 
% and where the figures should go.
out_dir =  ''; % where your data output from a_MOGS_MOnGs_calc.m went
fig_dir = '';

% Notes:
% 1) The first figure generated is a heatmap of the number of correctly
%    performed sequences (seq_score) performed for all participants. This
%    figure is not used in the report but is helpful for identifying the
%    problematic subject in the SEM group.
% 2) The second figure generated is a learning curve plot for all three
%    groups using the number of correctly performed sequences (seq_score) as
%    the dependent variable.
% 3) The program automatically invokes a few sub-functions that do not need 
%    to be run by the user in addition to b1_acquisition_retention. They must 
%    be in the same folder as the b1_acquisition_retent.m script. They are 
%    "g2_curvefit_mogs.m," "g3_exponentialfit.m," "g4_exponential_ls.m," 
%    and "g5_exponential.ms." These functions will determine the
%    acquisition asymptote, magnitude, and rate for all participants using
%    their learning curve described in 2.
% 4) One-way ANOVAs contrasting groups are performed for acquisition
%    asymptote, magnitude, and rate, and retention.
%
%    Output:
%    Correct_sequence_data.txt          --> Input to curvefit program
%    Correct_sequence_data_output.txt   --> Output of curvefit program
%    a_data.mat                         --> asymptote data
%    b_data.mat                         --> magnitude data
%    c_data.mat                         --> rate data
%    diff_perf.mat                      --> macro-online gains (From trials 24-25)
%
% Plots:
% 1) Figure 3A
% 2) Figure 3B

%% Load in data
load([out_dir, '/Microoffline_data.mat'])
load([out_dir, '/Microonline_data.mat'])
load([out_dir, '/ITI.mat'])
load([out_dir, '/Number_of_correct_sequences.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('         Started Acquisition and Retention Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');

%% Quality check of participant data
% quality check of data (Subject 10 in the SEM group has many trials with
% no correct completed sequences

for grp = 1:3
    subplot(1,3,grp); 
    if grp == 1; ttl = 'REST';elseif grp==2;ttl='ENC'; else;ttl='SEM';end
    imagesc(seq_score(:,1:24,grp));
    title(ttl,'FontSize',36);
    xlabel('trial','FontSize',36)
    if grp == 1; ylabel('Participant','FontSize',16);end
    xtk = get(gca,'XTickLabel'); set(gca,'XTickLabel',xtk,'fontsize',24); 
    ytk = get(gca,'YTickLabel'); set(gca,'XTickLabel',ytk,'fontsize',24); 
end

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
sgt = sgtitle('Number of correct sequences by gorup and participant','Color','black');
sgt.FontSize = 36;
disp('press any button to continue...')
pause;close;
disp(repmat(char(8), 1, 32));

% Remove subject with poor performance
seq_score(10,:,3) = NaN; 

%% Plot the mean number of correctly performed sequences for each group

% Isolate acquisition data
numb_corr_performance = seq_score(:,1:24,:);

% Isolate retention data
numb_corr_retention = seq_score(:,25:33,:);

% Calculate mean performance (1-24) for all three groups
mean_perf = zeros(3,24);
sem_perf = zeros(3,24);
for grp=1:3 %group
    for trl=1:24 % trial
        mean_perf(grp,trl) = mean(numb_corr_performance(:,trl,grp), 'omitnan');
        sem_perf(grp,trl) = std(numb_corr_performance(:,trl,grp),'omitmissing')/sqrt(size(numb_corr_performance,1));
    end
end

% Plot performance during acquisition trials
x = 1:24;
y = mean_perf;
err = sem_perf;
shadedErrorBar(x(1,:),y(1,:),err(1,:),'lineprops',{'-ko','MarkerFaceColor','k','MarkerSize',12});
hold on
shadedErrorBar(x(1,:),y(2,:),err(2,:),'lineprops',{'-bo','MarkerFaceColor','b','MarkerSize',12});
shadedErrorBar(x(1,:),y(3,:),err(3,:),'lineprops',{'-ro','MarkerFaceColor','r','MarkerSize',12});
ylim([0 10]);
xlim([1 24]);
xlabel('trial','FontSize',36);
ylabel(['Mean Number of correctly ' newline 'performed sequences'],'FontSize',36)
set(gca,'FontSize',42,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
title([' Effect of group on number of correct sequences (acquisition) ' newline]); 
legend('Rest','ENC','SEM','Location','southeast');
hold off

% Save figure
savefig([fig_dir, '\Figure_3A.fig']);
exportgraphics(gca, [fig_dir, '\Figure_3A.tiff'], 'Resolution', 300)

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

%% Fit a three parameter exponential growth curve to participant's data

% Organize data
clc;
disp('');
disp(' -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-');
disp(' Number of correct sequences analysis')
disp(' -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-');
disp('');
disp(' --- ANOVAs comparing learning parameters between groups (# corr seq ) --- ');

t=1;
cf_temp = zeros(1056,4);
for i=1:3 % Groups
    temp_cfit_data = seq_score(:,1:24,i);
    for j=1:15 % participants
        for k=1:length(temp_cfit_data)
            if and(j==10,i==3); continue;end
            cf_temp(t,1) = j;
            cf_temp(t,2) = i;
            cf_temp(t,3) = k;
            cf_temp(t,4) = temp_cfit_data(j,k);
            t=t+1;
        end
    end
    clear temp_cfit_data
end

% Use curvefit progam to get acquisition asymptotes, magnitudes, and rates.
writematrix(cf_temp,[out_dir, '/Correct_sequence_data.txt'],'Delimiter','tab');
load([out_dir, '/Correct_sequence_data.txt'])
h2_curvefit_mogs;  
writematrix(output,[out_dir, '/' outfile])
cf_output = readmatrix([out_dir, '/Correct_sequence_data_output.txt']);

% Put nans in the place of SEM subject 10's data
cf_output =[cf_output(1:39,:); [nan nan nan nan nan nan nan]; cf_output(40:end,:)];

% Perform ANOVAs for asymptote, magnitude, and rate.
a_data = [cf_output(1:15,3) cf_output(16:30,3) cf_output(31:45,3)];
b_data = [cf_output(1:15,4) cf_output(16:30,4) cf_output(31:45,4)];
c_data = [cf_output(1:15,5) cf_output(16:30,5) cf_output(31:45,5)];
aov_a = anova(a_data); pval_a = aov_a.stats{1,5}; disp(['One-way ANOVA comparing asymptote between groups: F(' num2str(aov_a.stats{1,2}) ',' num2str(aov_a.stats{2,2}) ') = ' num2str(aov_a.stats{1,4}) ', p = ' num2str(pval_a) '.']);
aov_b = anova(b_data); pval_b = aov_b.stats{1,5}; disp(['One-way ANOVA comparing magnitude between groups: F(' num2str(aov_b.stats{1,2}) ',' num2str(aov_b.stats{2,2}) ') = ' num2str(aov_b.stats{1,4}) ', p = ' num2str(pval_b) '.']);
aov_c = anova(c_data); pval_c = aov_c.stats{1,5}; disp(['One-way ANOVA comparing rate between groups: F(' num2str(aov_c.stats{1,2}) ',' num2str(aov_c.stats{2,2}) ') = ' num2str(aov_c.stats{1,4}) ', p = ' num2str(pval_c) '.']);
disp('press any button to continue...');
pause;
disp(repmat(char(8), 1, 32));

%% Analyze retention data

disp(' --- ANOVAs comparing retention between groups (# corr seq ) --- ');
% Calculate mean retention across groups
pre_perf = zeros(3,15);
post_perf = zeros(3,15);
for grp=1:3 %group
    for subject=1:15 % subject
        pre_perf(grp,subject) = mean(seq_score(subject,19:24,grp),'omitnan');
        post_perf(grp,subject) = mean(seq_score(subject,25:33,grp),'omitnan');
    end
end

% Organize data for violinplot
pre_perf = pre_perf';
post_perf = post_perf';
diff_perf = post_perf-pre_perf;

% retention plot
a = violinplot(diff_perf);
a(1).ViolinColor='black';
a(1).ScatterPlot.SizeData=100;
a(1).ScatterPlot.MarkerFaceAlpha=0.9;
a(1).ScatterPlot.MarkerFaceColor='black';
a(1).ScatterPlot.MarkerEdgeColor='black';
a(1).EdgeColor='black';
a(1).BoxColor='black';
a(1).ShowMean=1;
a(1).MedianPlot.SizeData=150;
a(1).MedianPlot.MarkerEdgeColor=[0,0,0];
a(1).MedianPlot.LineWidth=2.5;
a(1).MeanPlot.LineWidth=2.5;
a(1).MeanPlot.Color=[0,0,0];
a(2).ViolinColor='blue';
a(2).ScatterPlot.SizeData=100;
a(2).ScatterPlot.MarkerFaceAlpha=0.9;
a(2).ScatterPlot.MarkerFaceColor='blue';
a(2).ScatterPlot.MarkerEdgeColor='black';
a(2).EdgeColor='black';
a(2).BoxColor='black';
a(2).ShowMean=1;
a(2).MedianPlot.SizeData=150;
a(2).MedianPlot.MarkerEdgeColor=[0,0,0];
a(2).MedianPlot.LineWidth=2.5;
a(2).MeanPlot.LineWidth=2.5;
a(2).MeanPlot.Color=[0,0,0];
a(3).ViolinColor='red';
a(3).ScatterPlot.SizeData=100;
a(3).ScatterPlot.MarkerFaceAlpha=0.9;
a(3).ScatterPlot.MarkerFaceColor='red';
a(3).ScatterPlot.MarkerEdgeColor='black';
a(3).EdgeColor='black';
a(3).BoxColor='black';
a(3).ShowMean=1;
a(3).MedianPlot.SizeData=150;
a(3).MedianPlot.MarkerEdgeColor=[0,0,0];
a(3).MedianPlot.LineWidth=2.5;
a(3).MeanPlot.LineWidth=2.5;
a(3).MeanPlot.Color=[0,0,0];
ylim([-1 5]);
yline(0, 'k--', 'LineWidth', 1)
set(gca,'FontSize',24);
ylabel(['Retention: ∆ number of correctly ' newline 'performed sequences'])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
xticklabels({'REST','ENC','SEM'}); title(' Effect of group on retention (# correct sequences) '); 

aov_retention = anova(diff_perf); pval_retention = aov_retention.stats{1,5}; disp(['One-way ANOVA comparing retention between groups: F(' num2str(aov_retention.stats{1,2}) ',' num2str(aov_retention.stats{2,2}) ') = ' num2str(aov_retention.stats{1,4}) ', p = ' num2str(pval_retention) '.']);
clearvars -except out_dir fig_dir a_data b_data back_dir c_data data diff_perf home_dir ITI mean_perf mem_perf mem_perf_sorted microonline microoffline seq_score numb_corr_performance numb_corr_retention seq_score  

% Save figure
savefig([fig_dir, '\Figure_3B.fig']);
exportgraphics(gca, [fig_dir, '\Figure_3B.tiff'], 'Resolution', 300)

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));



%% Save data and Finish up
save([out_dir, '/a_data.mat'],'a_data','-v7.3');
save([out_dir, '/b_data.mat'],'b_data','-v7.3');
save([out_dir, '/c_data.mat'],'c_data','-v7.3');
save([out_dir, '/diff_perf.mat'],"diff_perf",'-v7.3')

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('         Finished Acquisition and Retention Analysis       ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');