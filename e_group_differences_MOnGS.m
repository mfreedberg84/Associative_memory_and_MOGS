clear; clc;
%% Analysis 4: Group differences in Mongs
%
% Requirements: Access to files generated from "a_MOGS_MOnGS_calc.m".
%
% indicate where your data output from a_MOGS_MOnGS_calc.m went, 
% where the figures should go, and the break point
out_dir = '';
fig_dir = '';
bp = 11;

% Tests: 
% 1) LME contrasting MOGS across groups
%
% Output: 
% Mongs_by_part.m
%
% Plots:
% 1) Figure 2C

%% Load in data
load([out_dir, '/Microonline_data.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('           Group differences in Mongs         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');

%% Linear mixed effects model at trial level (MOnGS)

% Only keep first 11 MOnGS
microonline10 = microonline(1:11,:, :); % Keep only the first 11 MOnGS

% take microonline and make into into a single column of data
trial_mongs = microonline10(:); % Reshape microonline into a single column

% Create independent variable for ANOVAn for Block
trl_covar = 1:11;
blk_iv = repmat(trl_covar, 1, 45)';

% Create independent variable for ANOVAn for GROUP
Groups = [1,2,3];
Group_iv = repelem(Groups, 1, 165)'; % Repeats each element 3 times

% Create subject variable for ANOVAn for GROUP
subs = 1:45;
subs_iv = repelem(subs,1, 11)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl = [Group_iv subs_iv blk_iv trial_mongs]; 

% Remove subject 40 from tbl
tbl(subs_iv == 40,:) = [];

% Prepare the table for LME analysis
tbl = array2table(tbl, 'VariableNames', ...
    {'Group', 'Subject', 'Block', 'Microonline' ...
    });

% Run LME Contrasting BLOCK and GROUP 
lme = fitlme(tbl,'Microonline~Group*Block+(1|Subject)');

% Group effect
disp(['Group effect: t(', ...
    num2str(lme.Coefficients{2,5}), ...
    ') = ', ...
    num2str(round(lme.Coefficients{2,4},2)), ...
    ', p = ', ...
    num2str(round(lme.Coefficients{2,6},3)) ...
    ]);

% Trial effect
disp(['Trial effect: t(', ...
    num2str(lme.Coefficients{3,5}), ...
    ') = ', ...
    num2str(round(lme.Coefficients{3,4},2)), ...
    ', p = ', ...
    num2str(round(lme.Coefficients{3,6},3)) ...
    ]);

% Interaction
disp(['Group x Trial: t(', ...
    num2str(lme.Coefficients{4,5}), ...
    ') = ', ...
    num2str(round(lme.Coefficients{4,4},2)), ...
    ', p = ', ...
    num2str(round(lme.Coefficients{4,6},3)) ...
    ]);

disp(' ');
disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

%% Follow up significant interaction by splitting groups

% --- REST ---
microonline10_rest = microonline(1:11,:, 1); % Keep only the first 11 MOnGS

% take microonline and make into into a single column of data
trial_mongs_rest = microonline10_rest(:); % Reshape microonline into a single column

% Create independent variable for ANOVAn for Block
blk_iv_rest = repmat(trl_covar, 1, 15)';

% Create subject variable for ANOVAn for GROUP
subs = 1:15;
subs_iv_rest = repelem(subs,1, 11)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl_rest = [subs_iv_rest blk_iv_rest trial_mongs_rest]; 

% Prepare the table for LME analysis
tbl_rest = array2table(tbl_rest, 'VariableNames', {'Subject', 'Block', 'Microonline'});

% Run LME Contrasting BLOCK and GROUP 
lme_rest = fitlme(tbl_rest,'Microonline~Block+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME (rest):');
% Trial effect
disp(['Trial effect: t(', ...
    num2str(lme_rest.Coefficients{2,5}), ...
    ') = ', ...
    num2str(round(lme_rest.Coefficients{2,4},2)), ...
    ', p = ', ...
    num2str(round(lme_rest.Coefficients{2,6},3)) ...
    ]);

disp(' ');
disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

% --- ENC ---
microonline10_enc = microonline(1:11,:, 2); % Keep only the first 11 MOnGS

% take microonline and make into into a single column of data
trial_mongs_enc = microonline10_enc(:); % Reshape microonline into a single column

% Create independent variable for ANOVAn for Block
blk_iv_enc = repmat(trl_covar, 1, 15)';

% Create subject variable for ANOVAn for GROUP
subs = 1:15;
subs_iv_enc = repelem(subs,1, 11)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl_enc = [subs_iv_enc blk_iv_enc trial_mongs_enc]; 

% Prepare the table for LME analysis
tbl_enc = array2table(tbl_enc, 'VariableNames', {'Subject', 'Block', 'Microonline'});

% Run LME Contrasting BLOCK and GROUP 
lme_enc = fitlme(tbl_enc,'Microonline~Block+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME (enc):');
% Trial effect
disp(['Trial effect: t(', ...
    num2str(lme_enc.Coefficients{2,5}), ...
    ') = ', ...
    num2str(round(lme_enc.Coefficients{2,4},2)), ...
    ', p = ', ...
    num2str(round(lme_enc.Coefficients{2,6},3)) ...
    ]);

disp(' ');
disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

% --- SEM (PICK UP HERE!) ---
microonline10_sem = microonline(1:11,:, 3); % Keep only the first 11 MOnGS

% Remove subject 10
microonline10_sem(:,10) = []; % Remove subject 10

% take microonline and make into into a single column of data
trial_mongs_sem = microonline10_sem(:); % Reshape microonline into a single column

% Create independent variable for ANOVAn for Block
blk_iv_sem = repmat(trl_covar, 1, 14)';

% Create subject variable for ANOVAn for GROUP
subs = 1:14;
subs_iv_sem = repelem(subs,1, 11)';

% Concatenate subs, Groups, blk, and microoffline vector
tbl_sem = [subs_iv_sem blk_iv_sem trial_mongs_sem]; 

% Prepare the table for LME analysis
tbl_sem = array2table(tbl_sem, 'VariableNames', {'Subject', 'Block', 'Microonline'});

% Run LME Contrasting BLOCK and GROUP 
lme_sem = fitlme(tbl_sem,'Microonline~Block+(1|Subject)');

% Display the results of the ANOVA
disp('Results of the LME (sem):');
% Trial effect
disp(['Trial effect: t(', ...
    num2str(lme_sem.Coefficients{2,5}), ...
    ') = ', ...
    num2str(round(lme_sem.Coefficients{2,4},2)), ...
    ', p = ', ...
    num2str(round(lme_sem.Coefficients{2,6},3)) ...
    ]);

disp(' ');
disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

%% Plot results for MOnGS across trials
plot_mongs_m = zeros(11,3);
plot_mongs_s = zeros(11,3);
for i=1:3
    if i == 3
        % Remove column 10
        tmp_sem = [microonline10(:,1:9,3), microonline10(:,11:15,3)];
        plot_mongs_m(:,3) = mean(tmp_sem,2,'omitmissing');
        plot_mongs_s(:,3) = std(tmp_sem,0,2,'omitmissing')/sqrt(14);
    else
        plot_mongs_m(:,i) = mean(microonline10(:,:,i),2,'omitmissing');
        plot_mongs_s(:,i) = std(microonline10(:,:,i),0,2,'omitmissing')/sqrt(15);
    end
end

x = 1:11;
y = plot_mongs_m';
err = plot_mongs_s';
shadedErrorBar(x(1,:),y(1,:),err(1,:),'lineprops',{'-ko','MarkerFaceColor','k','MarkerSize',12});
hold on
shadedErrorBar(x(1,:),y(2,:),err(2,:),'lineprops',{'-bo','MarkerFaceColor','b','MarkerSize',12});
shadedErrorBar(x(1,:),y(3,:),err(3,:),'lineprops',{'-ro','MarkerFaceColor','r','MarkerSize',12});
ylim([-8 8]);
xlim([1 10]);
xticks([0,2,4,6,8,10]);
xticklabels({'0','2','4','6','8','10'})
yticks([-8,-6,-4, -2,0,2,4,6,8])
yticklabels({'-8','-6','-4','-2','0','2','4','6','8'})
xlabel('Trial Number','FontSize',16);
title(['MOnGS across early learning trials'], fontsize = 24);
ylabel('MOnGS','FontSize',16)
%legend()
% Rename legend items
%legend({'REST', 'ENC', 'SEM'}, 'Location', 'Best','FontSize',16);
box off
hold off

% Save figure
savefig([fig_dir, '\Figure_2C.fig']);
exportgraphics(gca, [fig_dir, '\Figure_2C.tiff'], 'Resolution', 300)

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));
%% Finish up 

% Calculate average MOGS for first 10 breaks
avgMOnGS = mean(microonline10,1,'omitmissing'); % Calculate average MOGS for the first 10 breaks
Mongs_by_part = [avgMOnGS(1,:,1)' avgMOnGS(1,:,2)' avgMOnGS(1,:,3)'];
% Remove SEM participant 10
Mongs_by_part(10,3) = nan; % Remove the average for SEM participant 10
save([out_dir '/MOnGS_by_part.mat']',"Mongs_by_part",'-v7.3')

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('        Finished MOnGS group Analysis         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');

