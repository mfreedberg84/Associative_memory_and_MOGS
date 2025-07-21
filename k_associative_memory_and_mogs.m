clear; clc;
%% Analysis 4: Determine whether associative memory performance for the ENC group is related to micro-offline gains.
%
% Requirements: Access to "_processed_data" files generated from
% "a_MOGS_MOnGS_calc.m", "c_group_differences_MOGS.m", and
% "e_group_differences_MOnGS.m"
%
% indicate where your data output from a_MOGS_MOnGS_calc.m went  
% and where the figures should go
out_dir = '';
fig_dir = '';

% Notes:
%   1) This script analyzes the association between associative memory
%   performance (encoding success) and MOGS for the ENC and SEM groups.
%   2) Associations are tested using Pearson's correlations
%            
%   Plots:  
%   1) Figure 3E

%% Load data
load([out_dir,'/Mogs_by_part.mat'])
load([out_dir,'/Mongs_by_part.mat'])
load([out_dir,'/ENC_SEM_Memory_Performance'])


disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
disp('        Started associative memory and MOGS Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');

%% --- Examining correlation between associative memory performance and MOGS for ENC and SEM ---

% Calculate average mem performance for ENC and SEM groups
avgMemPerformance_ENC = mean(mem_perf_sorted(:,:, 1),1, 'omitnan')';
avgMemPerformance_SEM = mean(mem_perf_sorted(:,:, 2),1, 'omitnan')';


tmp_ep = [avgMemPerformance_ENC(:,1) Mogs_by_part(:,2)]; 
tmp_ep(any(isnan(tmp_ep), 2), :) = [];
[ENC_MOGS_r, ENC_MOGS_p] = corrcoef(tmp_ep(:,1), tmp_ep(:,2));
disp(['ENC - associative memory performance and MOGS: r(' num2str(length(tmp_ep(:,1))-2) ') = ' num2str(ENC_MOGS_r(1,2)) ', p = ' num2str(ENC_MOGS_p(1,2))]); 
scatter(tmp_ep(:,1), tmp_ep(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); 
hold on; tlp = polyfit(tmp_ep(:,1), tmp_ep(:,2), 1); 
tlpx = [min(tmp_ep(:,1)) max(tmp_ep(:,1))]; 
tlpy = polyval(tlp, tlpx); 
plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 

tmp_sem = [avgMemPerformance_SEM(:,1) Mogs_by_part(:,3)];
tmp_sem(any(isnan(tmp_sem), 2), :) = [];
[SEM_MOGS_r, SEM_MOGS_p] = corrcoef(tmp_sem(:,1), tmp_sem(:,2));
disp(['SEM - associative memory performance and MOGS: r(' num2str(length(tmp_sem(:,1))-2) ') = ' num2str(SEM_MOGS_r(1,2)) ', p = ' num2str(SEM_MOGS_p(1,2))]); 
scatter(tmp_sem(:,1), tmp_sem(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
hold on; tlp = polyfit(tmp_sem(:,1), tmp_sem(:,2), 1); 
tlpx = [min(tmp_sem(:,1)) max(tmp_sem(:,1))]; 
tlpy = polyval(tlp, tlpx); 
plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
hold off; 

ylim([-6 8]);
xlim([0 1]);
xticks([0 0.2 0.4 0.6 0.8 1.0])
%xticklabels({'','','','','','','','','',''})
yticks([-6 -4 -2 0 2 4 6 8])
%yticklabels({'','','','','','','',''})
yline(0, 'k--', 'LineWidth', 3)
set(gca,'FontSize',26);
title(['Association between associative memory' newline 'performance and MOGS by group'],'Fontsize',16);
xlabel('associative memory performance','FontSize',16)
ylabel('MOGS','FontSize',16)
set(gca,'FontSize',16,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
set(gcf,'Color','white')

% Save figure
savefig([fig_dir, '\Figure_3E.fig']);
exportgraphics(gca, [fig_dir, '\Figure_3E.tiff'], 'Resolution', 300)

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

%% Finish up

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
disp('       Finished associative memory and MOGS Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');

clearvars -except back_dir ep_corr_a sm_corr_a home_dir bp_24_off 