clear; clc;
%% Analysis 3: Determine whether associative memory performance is associated with ESL and MOGS
%
% Requirements: Access to "_processed_data" files generated from "c1_acquisition_retention.m".
%
% The paths to these folders must be specified here:
home_dir = ''; %root_dir
load_dir = '';
out_dir = '';

% Notes:
%   1) This script uses the curvefit and retention data from the b1 script. 
%   2) Each of these are regressed against associative memory performance.
%   3) associative memory performance is the sum of recalled associative memories
%      and ranges from 1-24.
%   4) On all plots correlations for ENC (blue) and SEM (red) are shown,
%      along with trend lines
%   5) Statistics are output on command line.
%
%   Output
%   associative_Memory_ENC.mat         --> ep_corr_a (number of associative
%   memories recalled and asymptote data for ENC group)
%   associative_Memory_SEM.mat         --> sm_corr_a (number of associative
%   memories recalled and asymptote data for SEM group)

%% Load in data
load([load_dir, '/ENC_SEM_Memory_Performance.mat'])
load([load_dir, '/a_data.mat'])
load([load_dir, '/b_data.mat'])
load([load_dir, '/c_data.mat'])
load([load_dir, '/diff_perf.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('         Started associative Memory and Acquisition Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');

%% Examining correlation between associative memory performance 
% and ABR/Macro-offline gains for ENC

% associative Group

disp(' *** Examining correlation between associative memory performance and ABR/Macro-offline gains for ENC *** ');

ep_corr_a = zeros(15,2);
ep_corr_b = zeros(15,2);
ep_corr_c = zeros(15,2);
ep_corr_ret = zeros(15,2);
for subject=1:15 % associative group
    ep_corr_a(subject,1:2) = [sum(mem_perf_sorted(:,subject,1)) a_data(subject,2)];
    ep_corr_b(subject,1:2) = [sum(mem_perf_sorted(:,subject,1)) b_data(subject,2)];
    ep_corr_c(subject,1:2) = [sum(mem_perf_sorted(:,subject,1)) c_data(subject,2)];
    ep_corr_ret(subject,1:2) = [sum(mem_perf_sorted(:,subject,1)) diff_perf(subject,2)];
end

tmp_sm = zeros(15,1);
sm_corr_a = zeros(15,2);
sm_corr_b = zeros(15,2);
sm_corr_c = zeros(15,2);
sm_corr_ret = zeros(15,2);
for subject=1:15 % Semantic group
    if subject == 10
        continue
    else
        tmp_sm(subject,1) = sum(mem_perf_sorted(:,subject,2));
        sm_corr_a(subject,1:2) = [sum(mem_perf_sorted(:,subject,2)) a_data(subject,3)];
        sm_corr_b(subject,1:2) = [sum(mem_perf_sorted(:,subject,2)) b_data(subject,3)];
        sm_corr_c(subject,1:2) = [sum(mem_perf_sorted(:,subject,2)) c_data(subject,3)];
        sm_corr_ret(subject,1:2) = [sum(mem_perf_sorted(:,subject,2)) diff_perf(subject,3)];
    end
end

%% Plot associative memory performance by asymptote
sm_corr_a(10,:)=[];
scatter(ep_corr_a(:,1), ep_corr_a(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); hold on; 
tlp = polyfit(ep_corr_a(:,1), ep_corr_a(:,2), 1); tlpx = [min(ep_corr_a(:,1)) max(ep_corr_a(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 
scatter(sm_corr_a(:,1), sm_corr_a(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
tlp = polyfit(sm_corr_a(:,1), sm_corr_a(:,2), 1); tlpx = [min(sm_corr_a(:,1)) max(sm_corr_a(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
ylim([0 10]);
xlim([0 24]);
set(gca,'FontSize',24);
title(['Association between associative memory' newline 'performance and asymptote'],'Fontsize',54);
xlabel('associative memory performance','FontSize',36)
ylabel('Learning Asymptote','FontSize',36)
set(gca,'FontSize',42,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('','ENC', "",'SEM','Location','southeast');

% Statistical testing on correlations(Pearson's)
[ENC_ep_a_r, ENC_ep_a_p] = corrcoef(ep_corr_a(:,1), ep_corr_a(:,2));
disp(['ENC: associative memory performance and asymptote: r(' num2str(length(ep_corr_a(:,1))-2) ') = ' num2str(ENC_ep_a_r(1,2)) ', p = ' num2str(ENC_ep_a_p(1,2))]); 
[SEM_ep_a_r, SEM_ep_a_p] = corrcoef(sm_corr_a(:,1), sm_corr_a(:,2));
disp(['SEM: associative memory performance and asymptote: r(' num2str(length(sm_corr_a(:,1))-2) ') = ' num2str(SEM_ep_a_r(1,2)) ', p = ' num2str(SEM_ep_a_p(1,2))]); 
disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

%% Plot associative memory performance by magnitude
sm_corr_b(10,:)=[];
scatter(ep_corr_b(:,1), ep_corr_b(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); hold on; 
tlp = polyfit(ep_corr_b(:,1), ep_corr_b(:,2), 1); tlpx = [min(ep_corr_b(:,1)) max(ep_corr_b(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 
scatter(sm_corr_b(:,1), sm_corr_b(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
tlp = polyfit(sm_corr_b(:,1), sm_corr_b(:,2), 1); tlpx = [min(sm_corr_b(:,1)) max(sm_corr_b(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
ylim([0 10]);
xlim([0 24]);
set(gca,'FontSize',24);
title(['Association between associative memory' newline 'performance and magnitude'],'Fontsize',54);
xlabel('associative memory performance','FontSize',36)
ylabel('Learning Magnitude','FontSize',36)
set(gca,'FontSize',42,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('','ENC', "",'SEM','Location','southeast');

% Statistical testing on correlations(Pearson's)
[ENC_ep_b_r, ENC_ep_b_p] = corrcoef(ep_corr_b(:,1), ep_corr_b(:,2));
disp(['ENC: associative memory performance and magnitude: r(' num2str(length(ep_corr_b(:,1))-2) ') = ' num2str(ENC_ep_b_r(1,2)) ', p = ' num2str(ENC_ep_b_p(1,2))]); 
[SEM_ep_b_r, SEM_ep_b_p] = corrcoef(sm_corr_b(:,1), sm_corr_b(:,2));
disp(['SEM: associative memory performance and magnitude: r(' num2str(length(sm_corr_b(:,1))-2) ') = ' num2str(SEM_ep_b_r(1,2)) ', p = ' num2str(SEM_ep_b_p(1,2))]); 
disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

%% Plot associative memory performance by rate
sm_corr_c(10,:)=[];
scatter(ep_corr_c(:,1), ep_corr_c(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); hold on; 
tlp = polyfit(ep_corr_b(:,1), ep_corr_b(:,2), 1); tlpx = [min(ep_corr_b(:,1)) max(ep_corr_b(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 
scatter(sm_corr_c(:,1), sm_corr_c(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
tlp = polyfit(sm_corr_c(:,1), sm_corr_c(:,2), 1); tlpx = [min(sm_corr_c(:,1)) max(sm_corr_c(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
ylim([0 10]);
xlim([0 24]);
set(gca,'FontSize',24);
title(['Association between associative memory' newline 'performance and rate'],'Fontsize',54);
xlabel('associative memory performance','FontSize',36)
ylabel('Learning Rate','FontSize',36)
set(gca,'FontSize',42,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('','ENC', "",'SEM','Location','southeast');

% Statistical testing on correlations(Pearson's)
[ENC_ep_c_r, ENC_ep_c_p] = corrcoef(ep_corr_c(:,1), ep_corr_c(:,2));
disp(['ENC: associative memory performance and rate: r(' num2str(length(ep_corr_c(:,1))-2) ') = ' num2str(ENC_ep_c_r(1,2)) ', p = ' num2str(ENC_ep_c_p(1,2))]); 
[SEM_ep_c_r, SEM_ep_c_p] = corrcoef(sm_corr_c(:,1), sm_corr_c(:,2));
disp(['SEM: associative memory performance and rate: r(' num2str(length(sm_corr_c(:,1))-2) ') = ' num2str(SEM_ep_c_r(1,2)) ', p = ' num2str(SEM_ep_c_p(1,2))]); 
disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

%% Plot associative memory performance by rate
sm_corr_ret(10,:)=[];
scatter(ep_corr_ret(:,1), ep_corr_ret(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[0 0 1], 'LineWidth',1.5); hold on; 
tlp = polyfit(ep_corr_ret(:,1), ep_corr_ret(:,2), 1); tlpx = [min(ep_corr_ret(:,1)) max(ep_corr_ret(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"b"); 
scatter(sm_corr_ret(:,1), sm_corr_ret(:,2),150,'MarkerEdgeColor',[0 0 0], 'MarkerFaceColor',[1 0 0], 'LineWidth',1.5); 
tlp = polyfit(sm_corr_ret(:,1), sm_corr_ret(:,2), 1); tlpx = [min(sm_corr_ret(:,1)) max(sm_corr_ret(:,1))]; tlpy = polyval(tlp, tlpx); plot(tlpx, tlpy, 'LineWidth', 6,'Color',"r"); 
ylim([0 6]);
xlim([0 24]);
set(gca,'FontSize',24);
title(['Association between associative memory' newline 'performance and retention'],'Fontsize',54);
xlabel('associative memory performance','FontSize',36)
ylabel(['Retention magnitude' newline],'FontSize',36)
set(gca,'FontSize',42,'XColor','k','YColor','k');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
legend('','ENC', "",'SEM','Location','southeast');

% Statistical testing on correlations(Pearson's)
[ENC_ep_ret_r, ENC_ep_ret_p] = corrcoef(ep_corr_ret(:,1), ep_corr_ret(:,2));
disp(['ENC: associative memory performance and retention: r(' num2str(length(ep_corr_ret(:,1))-2) ') = ' num2str(ENC_ep_ret_r(1,2)) ', p = ' num2str(ENC_ep_ret_p(1,2))]); 
[SEM_ep_ret_r, SEM_ep_ret_p] = corrcoef(sm_corr_ret(:,1), sm_corr_ret(:,2));
disp(['SEM: associative memory performance and retention: r(' num2str(length(sm_corr_ret(:,1))-2) ') = ' num2str(SEM_ep_ret_r(1,2)) ', p = ' num2str(SEM_ep_ret_p(1,2))]); 
disp('press any button to continue...');
pause; close;
disp(repmat(char(8), 1, 32));

clearvars -except ep_corr_a ep_corr_b ep_corr_c ep_corr_ret sm_corr_a sm_corr_b sm_corr_c sm_corr_ret out_dir load_dir


%% Finish up
save([out_dir, '/associative_Memory_ENC.mat'],'ep_corr_a','-v7.3')
save([out_dir, '/associative_Memory_SEM.mat'],'sm_corr_a','-v7.3')

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('        Finished associative Memory and Acquisition Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
