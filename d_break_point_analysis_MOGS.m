clear; clc;
%% Analysis 3: Break point analysis (MOGS)
%
% Requirements: Access to files generated from "a_MOGS_MOnGS_calc.m".
%
% indicate where your data output from a_MOGS_MOnGS_calc.m went 
out_dir = '';

% Tests: 
% 1) LME contrasting MOGS across groups with trial number as covariate and
% subject as random factor
%
% Plots:
% 1) Figure 2B

%% Load in data
load([out_dir, '/Microoffline_data.mat'])
load([out_dir, '/Microonline_data.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('              Break Point Analysis           ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');

% Breakpoint analysis (MOGS/MONS)
bp_groups = zeros(15, 24, 3);
for grp=1:3
    for subject = 1:15
        for bp = 1:24 % all possible breakpoints
            bp_groups(subject,bp,grp)= mean( ...
                microoffline(1:bp-1,subject,grp),1, 'omitmissing' ...
                );
        end
    end
end
bp_groups(:,1,:) = [];

%% Testing whether results are dependent on break point
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
tbl = array2table(tbl, 'VariableNames', ...
    {'Group', 'Subject', 'Trial_numb', 'Microoffline' ...
    });

% Run LME Contrasting BLOCK and GROUP with a BLOCK COVARIATE
lme = fitlme(tbl,'Microoffline~Group+Trial_numb+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME:');

% Group effect
disp(['Group effect: t(', ...
    num2str(lme.Coefficients{2,5}), ...
    ') = ', ...
    num2str(round(lme.Coefficients{2,4},2)), ...
    ', p = ', ...
    num2str(round(lme.Coefficients{2,6},3)) ...
    ]);

% Trial number
disp(['Trial effect: t(', ...
    num2str(lme.Coefficients{3,5}), ...
    ') = ', ...
    num2str(round(lme.Coefficients{3,4},2)), ...
    ', p = ', ...
    num2str(round(lme.Coefficients{3,6},3)) ...
    ]);

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

clearvars -except bp_off microonline microoffline aov_bp_mogs aov_tmp out_dir bp bp_on bp_microoffline bp_groups

%% Finish up 

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('        Finished Break Point Analysis        ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');

