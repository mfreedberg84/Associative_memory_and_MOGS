clear; clc;
%% Analysis 8: post hoc sample size power analysis
%
% Requirements: Access to files generated from "a_MOGS_MOnGS_calc.m".
%
% indicate where your data output from a_MOGS_MOnGS_calc.m went 
out_dir = ''; %root_dir;

% Notes: 
% 1) The effect size generated from the group effect comes from 10000
% permuted simulations.
%
% Plots:
% None

%% Load in data

load([out_dir, '/Microoffline_data.mat'])
load([out_dir, '/Microonline_data.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('        Started Sample Size Analysis         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');
%% Linear mixed effects model at trial level (MOGS)

% Only keep first 10 MOGS
microoffline10 = microoffline(1:10,:, :); % Keep only the first 10 MOGS

% take microoffline and make into into a single column of data
trial_mogs = microoffline10(:); % Reshape data

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
tbl = array2table(tbl, 'VariableNames', ...
    {'Group', 'Subject', 'Block', 'Microoffline'} ...
    );

% Run LME Contrasting BLOCK and GROUP 
lme = fitlme(tbl,'Microoffline~Group*Block+(1|Subject)');

significance = anova(lme, ...
        'DFMethod', ...
        'Satterthwaite'); 

disp('Observed results')
disp(significance)

% calculate effect size from lme results
% Extract the coefficients and calculate the effect size
coefficients = lme.Coefficients.Estimate;
effectSize = coefficients(2) / std(tbl.Microoffline,'omitmissing');
disp(['Effect Size: ', num2str(effectSize)]);

% Calculate required sample size from effectSize
requiredSampleSize = ceil((norminv(0.95) + norminv(0.85))^2 * (std(tbl.Microoffline, 'omitmissing')^2) / effectSize^2);
disp(['Required Sample Size: ', num2str(requiredSampleSize)]);

%% Finishing up
disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('        Finished Sample Size Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');