clear; clc;
%% Analysis 8: post hoc sample size power analysis
%
% Requirements: Access to files generated from "a_MOGS_MOnGS_calc.m".
%
% indicate where your data output from a_MOGS_MOnGS_calc.m went 
out_dir = '';

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

% simulations to get power
% Set the number of simulations for power analysis
numSimulations = 1000;

% Initialize an array to store the results of each simulation
powerResults = zeros(numSimulations, 1);

rng(42);
hw = waitbar(0,'Running simulations...');
for sim = 1:numSimulations

    % Progress bar
    if mod(sim,floor(numSimulations/10))<1e-2
         waitbar(sim/numSimulations,hw);
    end

    % Perform the simulation and store the results
    perm_order = randperm(440);
    perm_mogs = tbl.Microoffline(perm_order);
    tmp_tbl = tbl;
    tmp_tbl.Microoffline = perm_mogs;
    lmeSim = fitlme(tmp_tbl, ...
        'Microoffline ~ Group*Block + (1|Subject)' ...
        );

    powerResults(sim) = anova(lmeSim, ...
        'DFMethod', ...
        'Satterthwaite').pValue(1) < 0.05; 
    len_disp = floor(log10(abs(sim))) + 1;
end
close(findall(groot,'type','figure'))

% Calculate the power as the proportion of significant results
powerEstimate = mean(powerResults);
disp(['Estimated Power: ', num2str(powerEstimate)]);

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