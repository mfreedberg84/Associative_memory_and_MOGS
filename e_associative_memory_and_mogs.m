clear; clc;
%% Analysis 5: Determine whether associative memory performance for the ENC group is related to micro-offline gains.
%
% Requirements: Access to "_processed_data" files generated from "a_MOGS_calc.m" and "e_MOGS_group_analysis.m.
%
% The paths to these folders must be specified here:

home_dir =  ''; %root_dir
load_dir = '';
out_dir = '';

% Notes:
%   1) This script analyzes the association between associative memory
%   performance (encoding success) and MOGS for the ENC and SEM groups.
%   2) Associations are tested using Pearson's correlations
%            
%   Output:  
%   None

%% Load data
load([load_dir,'/bp_off.mat'])
load([load_dir,'/Associative_Memory_ENC.mat'])
load([load_dir,'/Associative_Memory_SEM.mat'])


disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('        Started associative memory and MOGS Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');

% --- Examining correlation between associative memory performance and MOGS for ENC and SEM ---

tmp_ep = [ep_corr_a(:,1) bp_off(:,2)]; 
tmp_ep(any(isnan(tmp_ep), 2), :) = [];
[ENC_MOGS_r, ENC_MOGS_p] = corrcoef(tmp_ep(:,1), tmp_ep(:,2));
disp(['ENC - associative memory performance and MOGS: r(' num2str(length(tmp_ep(:,1))-2) ') = ' num2str(ENC_MOGS_r(1,2)) ', p = ' num2str(ENC_MOGS_p(1,2))]); 
scatter(tmp_ep(:,1), tmp_ep(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); 
hold on; tlp = polyfit(tmp_ep(:,1), tmp_ep(:,2), 1); 
tlpx = [min(tmp_ep(:,1)) max(tmp_ep(:,1))]; 
tlpy = polyval(tlp, tlpx); 
plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 

sm_correction = bp_off(:,3);
sm_correction(10) = [];
tmp_sem = [sm_corr_a(:,1) sm_correction(:,1)]; 
[SEM_MOGS_r, SEM_MOGS_p] = corrcoef(tmp_sem(:,1), tmp_sem(:,2));
disp(['SEM - associative memory performance and MOGS: r(' num2str(length(tmp_sem(:,1))-2) ') = ' num2str(SEM_MOGS_r(1,2)) ', p = ' num2str(SEM_MOGS_p(1,2))]); 
scatter(tmp_sem(:,1), tmp_sem(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
hold on; 
tlp = polyfit(tmp_sem(:,1), tmp_sem(:,2), 1); 
tlpx = [min(tmp_sem(:,1)) max(tmp_sem(:,1))]; 
tlpy = polyval(tlp, tlpx); 
plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
hold off; 

ylim([-6 8]);
xlim([6 24]);
xticks([6 8 10 12 14 16 18 20 22 24])
%xticklabels({'','','','','','','','','',''})
yticks([-6 -4 -2 0 2 4 6 8])
%yticklabels({'','','','','','','',''})
yline(0, 'k--', 'LineWidth', 3)
set(gca,'FontSize',24);
title(['Association between associative memory' newline 'performance and MOGS by group'],'Fontsize',54);
xlabel('associative memory performance','FontSize',36)
ylabel('MOGS','FontSize',36)
set(gca,'FontSize',42,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
disp('press any button to continue...');
set(gcf,'Color','white')
frame = getframe(gcf); % Capture the current figure
img = frame.cdata;
imwrite(img, 'MOGS.tiff'); % Save as TIFF
pause;close;
disp(repmat(char(8), 1, 32));


%% Finish up

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('       Finished associative memory and MOGS Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');

clearvars -except back_dir ep_corr_a sm_corr_a home_dir bp_24_off 