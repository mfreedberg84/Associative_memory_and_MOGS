clear; clc;
%% Analysis 9: Acquisition and Retention correlations with memory performance
%
% Requirements: Access to "_processed_data" files generated from "h1_acquisition_retention.m".
%
% indicate where your data output from a_MOGS_MOnGS_calc.m went  
% and where the figures should go
out_dir = '';
fig_dir = '';

% Notes:
%   1) This script uses the curvefit and retention data from the h1 script. 
%   2) Each of these are regressed against associative memory performance.
%   3) associative memory performance is the average of recalled associative memories
%      and ranges from 1-24.
%   4) On all plots correlations for ENC (blue) and SEM (red) are shown,
%      along with trend lines
%   5) Statistics are output on command line.
%
%   Plots
%   1) Figure 3A
%   2) Figure 3B
%   3) Figure 3C
%   4) Figure 3D

%% Load in data
load([out_dir '\ENC_SEM_Memory_Performance.mat'])
load([out_dir '\a_data.mat'])
load([out_dir '\b_data.mat'])
load([out_dir '\c_data.mat'])
load([out_dir '\diff_perf.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
disp('         Started associative Memory and Acquisition Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');

%% Plotting memory performance for each group

% Calculate average mem performance for ENC and SEM groups
avgMemPerformance_ENC = mean(mem_perf_sorted(:,:, 1),1, 'omitnan')';
avgMemPerformance_SEM = mean(mem_perf_sorted(:,:, 2),1, 'omitnan')';

% Remove subject 10 for SEM
avgMemPerformance_SEM(10,:) =nan; % Remove subject 10 from the performance data

% Make a violin plot of average memory performance for ENC and SEM
figure;
violinplot([avgMemPerformance_ENC, avgMemPerformance_SEM], {'ENC', 'SEM'});
ylabel('Average Memory Performance');
title('Violin Plot of Memory Performance for ENC and SEM Groups');
% Annotate average values
text(1, 0.6, ['Mean is ' num2str(mean(avgMemPerformance_ENC), '%.3f')], 'Color', 'b', 'FontSize', 16);
text(1, 0.55, ['SEM is ' num2str(std(avgMemPerformance_ENC)/sqrt(15), '%.3f')], 'Color', 'b', 'FontSize', 16);

text(2, 0.6, ['Mean is ' num2str(mean(avgMemPerformance_SEM,'omitmissing'), '%.3f')], 'Color', 'b', 'FontSize', 16);
text(2, 0.55, ['SEM is ' num2str(mean(avgMemPerformance_SEM,'omitmissing')/sqrt(14), '%.3f')], 'Color', 'b', 'FontSize', 16);

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));
%% Examining correlation between associative memory performance 

% Prepping data
enc_corr_a(:,1:2) = [avgMemPerformance_ENC a_data(:,2)];
enc_corr_b(:,1:2) = [avgMemPerformance_ENC b_data(:,2)];
enc_corr_c(:,1:2) = [avgMemPerformance_ENC c_data(:,2)];
enc_corr_ret(:,1:2) = [avgMemPerformance_ENC diff_perf(:,2)];

sem_corr_a(:,1:2) = [avgMemPerformance_SEM a_data(:,3)];
sem_corr_b(:,1:2) = [avgMemPerformance_SEM b_data(:,3)];
sem_corr_c(:,1:2) = [avgMemPerformance_SEM c_data(:,3)];
sem_corr_ret(:,1:2) = [avgMemPerformance_SEM diff_perf(:,3)];


%% Plot associative memory performance by asymptote
sem_corr_a(10,:) = [];
disp(' *** Examining correlation between associative memory performance and ABR/Macro-offline gain *** ');
scatter(enc_corr_a(:,1), enc_corr_a(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); hold on; 
tlp = polyfit(enc_corr_a(:,1), enc_corr_a(:,2), 1); tlpx = [min(enc_corr_a(:,1)) max(enc_corr_a(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 
scatter(sem_corr_a(:,1), sem_corr_a(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
tlp = polyfit(sem_corr_a(:,1), sem_corr_a(:,2), 1); tlpx = [min(sem_corr_a(:,1)) max(sem_corr_a(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
ylim([0 10]);
xlim([0 1]);
set(gca,'FontSize',16);
title(['Association between associative memory' newline 'performance and asymptote'],'Fontsize',54);
xlabel('associative memory performance','FontSize',16)
ylabel('Learning Asymptote','FontSize',16)
set(gca,'FontSize',16,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('','ENC', "",'SEM','Location','southeast');
hold off

% Statistical testing on correlations(Pearson's)
[ENC_ep_a_r, ENC_ep_a_p] = corrcoef(enc_corr_a(:,1), enc_corr_a(:,2));
disp(['ENC: associative memory performance and asymptote: r(' num2str(length(enc_corr_a(:,1))-2) ') = ' num2str(ENC_ep_a_r(1,2)) ', p = ' num2str(ENC_ep_a_p(1,2))]); 
[SEM_ep_a_r, SEM_ep_a_p] = corrcoef(sem_corr_a(:,1), sem_corr_a(:,2));
disp(['SEM: associative memory performance and asymptote: r(' num2str(length(sem_corr_a(:,1))-2) ') = ' num2str(SEM_ep_a_r(1,2)) ', p = ' num2str(SEM_ep_a_p(1,2))]); 

% Save figure
savefig([fig_dir, '\Figure_3A.fig']);
exportgraphics(gca, [fig_dir, '\Figure_3A.tiff'], 'Resolution', 300)

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

%% Plot associative memory performance by magnitude
sem_corr_b(10,:) = [];
scatter(enc_corr_b(:,1), enc_corr_b(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); hold on; 
tlp = polyfit(enc_corr_b(:,1), enc_corr_b(:,2), 1); tlpx = [min(enc_corr_b(:,1)) max(enc_corr_b(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 
scatter(sem_corr_b(:,1), sem_corr_b(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
tlp = polyfit(sem_corr_b(:,1), sem_corr_b(:,2), 1); tlpx = [min(sem_corr_b(:,1)) max(sem_corr_b(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
ylim([0 10]);
xlim([0 1]);
set(gca,'FontSize',16);
title(['Association between associative memory' newline 'performance and magnitude'],'Fontsize',54);
xlabel('associative memory performance','FontSize',16)
ylabel('Learning Magnitude','FontSize',16)
set(gca,'FontSize',16,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('','ENC', "",'SEM','Location','southeast');

% Statistical testing on correlations(Pearson's)
[ENC_ep_b_r, ENC_ep_b_p] = corrcoef(enc_corr_b(:,1), enc_corr_b(:,2));
disp(['ENC: associative memory performance and magnitude: r(' num2str(length(enc_corr_b(:,1))-2) ') = ' num2str(ENC_ep_b_r(1,2)) ', p = ' num2str(ENC_ep_b_p(1,2))]); 
[SEM_ep_b_r, SEM_ep_b_p] = corrcoef(sem_corr_b(:,1), sem_corr_b(:,2));
disp(['SEM: associative memory performance and magnitude: r(' num2str(length(sem_corr_b(:,1))-2) ') = ' num2str(SEM_ep_b_r(1,2)) ', p = ' num2str(SEM_ep_b_p(1,2))]); 

% Save figure
savefig([fig_dir, '\Figure_3B.fig']);
exportgraphics(gca, [fig_dir, '\Figure_3B.tiff'], 'Resolution', 300)

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

%% Plot associative memory performance by rate
sem_corr_c(10,:) = [];
scatter(enc_corr_c(:,1), enc_corr_c(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); hold on; 
tlp = polyfit(enc_corr_b(:,1), enc_corr_b(:,2), 1); tlpx = [min(enc_corr_b(:,1)) max(enc_corr_b(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 
scatter(sem_corr_c(:,1), sem_corr_c(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
tlp = polyfit(sem_corr_c(:,1), sem_corr_c(:,2), 1); tlpx = [min(sem_corr_c(:,1)) max(sem_corr_c(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
ylim([0 10]);
xlim([0 1]);
set(gca,'FontSize',24);
title(['Association between associative memory' newline 'performance and rate'],'Fontsize',54);
xlabel('associative memory performance','FontSize',16)
ylabel('Learning Rate','FontSize',16)
set(gca,'FontSize',16,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('','ENC', "",'SEM','Location','southeast');

% Statistical testing on correlations(Pearson's)
[ENC_ep_c_r, ENC_ep_c_p] = corrcoef(enc_corr_c(:,1), enc_corr_c(:,2));
disp(['ENC: associative memory performance and rate: r(' num2str(length(enc_corr_c(:,1))-2) ') = ' num2str(ENC_ep_c_r(1,2)) ', p = ' num2str(ENC_ep_c_p(1,2))]); 
[SEM_ep_c_r, SEM_ep_c_p] = corrcoef(sem_corr_c(:,1), sem_corr_c(:,2));
disp(['SEM: associative memory performance and rate: r(' num2str(length(sem_corr_c(:,1))-2) ') = ' num2str(SEM_ep_c_r(1,2)) ', p = ' num2str(SEM_ep_c_p(1,2))]); 

% Save figure
savefig([fig_dir, '\Figure_3C.fig']);
exportgraphics(gca, [fig_dir, '\Figure_3C.tiff'], 'Resolution', 300)

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

%% Plot associative memory performance by rate
sem_corr_ret(10,:) = [];
scatter(enc_corr_ret(:,1), enc_corr_ret(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); hold on; 
tlp = polyfit(enc_corr_ret(:,1), enc_corr_ret(:,2), 1); tlpx = [min(enc_corr_ret(:,1)) max(enc_corr_ret(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 
scatter(sem_corr_ret(:,1), sem_corr_ret(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
tlp = polyfit(sem_corr_ret(:,1), sem_corr_ret(:,2), 1); tlpx = [min(sem_corr_ret(:,1)) max(sem_corr_ret(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
ylim([0 6]);
xlim([0 1]);
set(gca,'FontSize',16);
title(['Association between associative memory' newline 'performance and retention'],'Fontsize',54);
xlabel('associative memory performance','FontSize',16)
ylabel(['Retention magnitude' newline],'FontSize',16)
set(gca,'FontSize',16,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('','ENC', "",'SEM','Location','southeast');

% Statistical testing on correlations(Pearson's)
[ENC_ep_ret_r, ENC_ep_ret_p] = corrcoef(enc_corr_ret(:,1), enc_corr_ret(:,2));
disp(['ENC: associative memory performance and retention: r(' num2str(length(enc_corr_ret(:,1))-2) ') = ' num2str(ENC_ep_ret_r(1,2)) ', p = ' num2str(ENC_ep_ret_p(1,2))]); 
[SEM_ep_ret_r, SEM_ep_ret_p] = corrcoef(sem_corr_ret(:,1), sem_corr_ret(:,2));
disp(['SEM: associative memory performance and retention: r(' num2str(length(sem_corr_ret(:,1))-2) ') = ' num2str(SEM_ep_ret_r(1,2)) ', p = ' num2str(SEM_ep_ret_p(1,2))]); 

% Save figure
savefig([fig_dir, '\Figure_3D.fig']);
exportgraphics(gca, [fig_dir, '\Figure_3D.tiff'], 'Resolution', 300)

disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

clearvars -except enc_corr_a enc_corr_b enc_corr_c enc_corr_ret sm_corr_a sem_corr_b sem_corr_c sem_corr_ret out_dir load_dir


%% Finish up
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
disp('        Finished associative Memory and Acquisition Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
