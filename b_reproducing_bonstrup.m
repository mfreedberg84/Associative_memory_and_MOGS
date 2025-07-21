clear; clc;
%% Analysis 1:Reproducing Bonstrup et al. (2019)
%
% Requirements: Access to files generated from "a_MOGS_MOnGS_calc.m".
%
% indicate where your data output from a_MOGS_MOnGS_calc.m went and
% the break point
out_dir = '';
bp = 11;

% Tests: 
% 1) Average MOGS for first 11 trials against zero.
% 2) Average MOnGS for the first 11 trials against zero.
% 3) Average MOGS vs. MOnGS for first 11 trials
% 4) Observed Power analysis of MOGS against zero
%
% Plots:
% 1) Power analysis results

%% Load in data
load([out_dir, '/Microoffline_data.mat'])
load([out_dir, '/Microonline_data.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('         Analysis 1: Reproducing results of Bonstrup et al. (2019)         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');

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

%% Post-hoc power analysis

% Define parameters
testtype = 't';         % One-sample t-test
mu0 = 0;                % Mean under the null hypothesis (testing against zero)
sigma0 = 1;             % Estimated standard deviation under the null hypothesis
p0 = [mu0, sigma0];     % Null hypothesis parameters
mu1 = mean(rest_mogs);  % Mean under the alternative hypothesis (expected mean difference)
p1 = mu1;               % Alternative hypothesis parameter
pwr = [];              % Desired power
nn=1:100;

% Calculate required sample size
pwrout = sampsizepwr('t',[0 0.4117],0.3469,[],nn,'Alpha',0.05);

figure;
plot(nn,pwrout,'b-')
% Make vertical red line at n=15
line([15 15], ylim, 'Color', 'r', 'LineStyle', '--');
title('Power versus Sample Size')
xlabel('Sample Size')
ylabel('Power')
disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

%% Finishing up
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - -');
disp('            Finished reproducibility Analysis           ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - -');

