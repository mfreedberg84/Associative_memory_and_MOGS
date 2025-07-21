clear; clc;
%% Analysis 2: Group differences in Mogs
%
% Requirements: Access to files generated from "a_MOGS_MOnGS_calc.m".
%
% indicate where your data output from a_MOGS_MOnGS_calc.m went, 
% where the figures should go, and the break point
out_dir = '';
fig_dir = '';
mkdir(fig_dir);
bp = 11;

% Output
% 1) Mogs_by_part.m
%
% Tests: 
% 1) LME contrasting MOGS across groups
%
% Plots:
% 1) Figure 2A

%% Load in data
load([out_dir, '/Microoffline_data.mat'])
load([out_dir, '/Microonline_data.mat'])

disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('           Group differences in Mogs         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');

%% Linear mixed effects model at trial level (MOGS)

% Only keep first 10 MOGS
microoffline10 = microoffline(1:10,:, :); % Keep only the first 10 MOGS

% take microoffline and make into into a single column of data
trial_mogs = microoffline10(:); % Reshape microoffline into a single column

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
tbl = array2table(tbl, 'VariableNames', {'Group', 'Subject', 'Block', 'Microoffline'});

% Run LME Contrasting BLOCK and GROUP 
lme = fitlme(tbl,'Microoffline~Group*Block+(1|Subject)');

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

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));

%% Plot results for all blocks

% Calculate MOGs for all intertrial periods
microoffline23 = microoffline(1:23,:, :);

plot_mogs_m = zeros(23,3);
plot_mogs_s = zeros(23,3);
for i=1:3
    if i == 3
        % Remove column 10
        tmp_sem = [microoffline23(:,1:9,3), microoffline23(:,11:15,3)];
        plot_mogs_m(:,3) = mean(tmp_sem,2,'omitmissing');
        plot_mogs_s(:,3) = std(tmp_sem,0,2,'omitmissing')/sqrt(14);
    else
        plot_mogs_m(:,i) = mean(microoffline23(:,:,i),2,'omitmissing');
        plot_mogs_s(:,i) = std(microoffline23(:,:,i),0,2,'omitmissing')/sqrt(15);
    end
end

x = 1:23;
y = plot_mogs_m';
err = plot_mogs_s';
shadedErrorBar(x(1,:),y(1,:),err(1,:),'lineprops',{'-ko','MarkerFaceColor','k','MarkerSize',12});
hold on
shadedErrorBar(x(1,:),y(2,:),err(2,:),'lineprops',{'-bo','MarkerFaceColor','b','MarkerSize',12});
shadedErrorBar(x(1,:),y(3,:),err(3,:),'lineprops',{'-ro','MarkerFaceColor','r','MarkerSize',12});
ylim([-8 8]);
xlim([1 23]);
xticks([1:23])
xticklabels({'1','','','','5','','','','','10','','','','','15','','','','','20','','','',''});
yticks([-8,-6,-4, -2,0,2,4,6,8])
yticklabels({'-8','-6','-4','-2','0','2','4','6','8'})
line([bp-0.5 bp-0.5], ylim, 'Color', 'k', 'LineStyle', '--');
% Create horizontal line at y=0
% Create a horizontal line at y=0
yline(0, 'k--');
xlabel('Trial Number','FontSize',16);
title(['Group differences in MOGS'], fontsize = 24);
ylabel('MOGS','FontSize',16)
legend()
% Rename legend items
legend({'REST', 'ENC', 'SEM'}, 'Location', 'Best','FontSize',24);
legend('boxoff');
box off
hold off

% Save figure
savefig([fig_dir, '\Figure_2A.fig']);
exportgraphics(gca, [fig_dir, '\Figure_2A.tiff'], 'Resolution', 300)

disp('press any button to continue...');
pause;close;
disp(repmat(char(8), 1, 32));
%% Finish up 

% Calculate average MOGS for first 10 breaks
avgMOGS = mean(microoffline10,1,'omitmissing'); % Calculate average MOGS for the first 10 breaks
Mogs_by_part = [avgMOGS(1,:,1)' avgMOGS(1,:,2)' avgMOGS(1,:,3)'];
save([out_dir '/MOGS_by_part.mat']',"Mogs_by_part",'-v7.3')


disp(' - - - - - - - - - - - - - - - - - - - - - - ');
disp('        Finished MOGS group Analysis         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - ');

