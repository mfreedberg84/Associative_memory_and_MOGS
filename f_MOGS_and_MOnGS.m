clear; clc;
%% Analysis 5: Correlation between MOGS and MOnGS
%
% Requirements: Access to "_processed_data" files generated 
% from "a_MOGS_MOnGS_calc.m" and "b_MOGS_MOnGS_group_analysis.m.
%
% The paths to these folders must be specified here:

home_dir =  ''; %root_dir
load_dir = '';
out_dir = '';

% Notes:
%  1) This script performs a simple Pearson's R correlation between 
%  MOGS and MOnGS for the first 11 trials.
%  2) Associations are tested using Pearson's correlations
%            
%  Output:None

%% Load data
load([load_dir,'/Microoffline_data.mat'])
load([load_dir,'/Microonline_data.mat'])


disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('            Correlation between MOGS and MOnGS           ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');

%% All participants

% Remove trials 12-32 for MOGS
MOGS_11 = microoffline(1:11,:,:);

% Remove trials 12-33 for MOnGS
MOnGS_11 =  microonline(1:11,:,:);

% Calculate average MOGS and MOnGS for the first 11 trials
avgMOGS = mean(MOGS_11, 2,'omitmissing');
avgMOnGS = mean(MOnGS_11, 2,'omitmissing');

% Further means across all subjects
% Calculate further means across all subjects for MOGS and MOnGS
avgMOGS_trl = mean(avgMOGS, 3, 'omitmissing');
avgMOnGS_trl = mean(avgMOnGS, 3, 'omitmissing');

% Regress MOGS against MOnGS
% Calculate the correlation coefficients between MOGS and MOnGS
[correlationR, correlationP] = corrcoef(avgMOGS_trl, avgMOnGS_trl);
disp(['All: Correlation between MOGS and MOnGS: r(' ...
    num2str(length(avgMOGS_trl)-2) ...
    ') = ' ...
    num2str(correlationR(1,2)) ...
    ', p = ' ...
    num2str(correlationP(1,2)) ...
    ]);

%% REST
% Regress MOGS against MOnGS (REST)
% Calculate the correlation coefficients between MOGS and MOnGS
[correlationR, correlationP] = corrcoef(avgMOGS(:,1,1), avgMOnGS(:,1,1));
disp(['REST: Correlation between MOGS and MOnGS: r(' ...
    num2str(length(avgMOGS(:,1,1))-2) ...
    ') = ' ...
    num2str(correlationR(1,2)) ...
    ', p = ' ...
    num2str(correlationP(1,2)) ...
    ]);

%% ENC

% Regress MOGS against MOnGS (ENC)
% Calculate the correlation coefficients between MOGS and MOnGS
[correlationR, correlationP] = corrcoef(avgMOGS(:,1,2), avgMOnGS(:,1,2));
disp(['ENC: Correlation between MOGS and MOnGS: r(' ...
    num2str(length(avgMOGS(:,1,1))-2) ...
    ') = ' ...
    num2str(correlationR(1,2)) ...
    ', p = ' ...
    num2str(correlationP(1,2)) ...
    ]);

%% SEM

% Regress MOGS against MOnGS (SEM)
% Calculate the correlation coefficients between MOGS and MOnGS
[correlationR, correlationP] = corrcoef(avgMOGS(:,1,3), avgMOnGS(:,1,3));
disp(['SEM: Correlation between MOGS and MOnGS: r(' ...
    num2str(length(avgMOGS(:,1,1))-2) ...
    ') = ' ...
    num2str(correlationR(1,2)) ...
    ', p = ' ...
    num2str(correlationP(1,2)) ...
    ]);


%% Finish up

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('              Finished MOGS/MOnGS Correlation            ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
