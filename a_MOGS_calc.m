clear; clc;
%% Analysis for Mogs_Hipp. 
% This script organizes the data for all analyses. Each row in MOGS_data.txt 
% represents a % single keypress during the ESL task or each trial during the
% associative memory recall test (at the end). 
%
% Requirements: Access to "MOGS_data.txt" File.
%
% The paths to these folders must be specified here:
 
root_dir = ''; % folder with your scripts
out_dir = ''; % output data folder
mkdir(out_dir);

% Notes
% 1) For each subject it will acces the data in "_raw_data".
% 2) Encoding success accuracy will be collected for all participants in
%    the ENC and SEM groups --> mem_perf_sorted. This is a file with 1s and 0s
%    represented the chronological accuracy in encoding the word pairs from
%    Trials 1-24 during acquisition.
% 3) For all 33 trials (24 acquisition and 9 retention), the 
%    performance of all correct sequences will be recorded in keypresses per
%    second (KPS) --> ITI.
% 4) The number of correct sequences will be recorded using two methods. In the
%    first method, correct partial sequences at the end of the trial will
%    count be patially added on to the running total, such that completing
%    four elements of the sequence adds 0.8 to the total, three elements adds
%    0.6, and two elements adds 0.4. This will be stored in the "seq_score" variable. 
%    In the seocnd method, any correct partial sequence at the end of the trial 
%    counts as a full sequence ("seq_count").  
% 5) Online gains ("microonline") will be calculated for all 33 
%    acquisition trials as the difference in KPS between the last correct 
%    sequence (identified with "seq_count") and the first one.
% 6) Offline gains ("microoffline") will be calculated for all 32 offline 
%    periods as the difference in KPS between the first correct sequence 
%    on trial N+1 and the last correct sequence on trial N 
%    (identified with "seq_count"). Note that the 24th offline period is 
%    the 30 minute rest period before retention.
%
%    Output:
%    Number_of_correct_sequences.mat --> seq_score
%    Microonline_data.mat            --> microonline
%    Microofline_data.mat            --> microoffline
%    ENC_SEM_Memory_Performance.mat  --> mem_perf_sorted
%    ITI.mat                         --> ITI


%% Main loop (organizing participant data and calculating MOGS)
% Read in data
MOGS_data = readtable([root_dir, '/MOGS_data.txt'],'VariableNamingRule', 'preserve');

disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');
disp('         Collecting and organizing performance metrics         ');
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ');

% preallocated matrices
seq_score = zeros(15,33,3); 
mem_perf_sorted = zeros(24,15);
ITI = zeros(50,24,15,3);
first_cor_seq = zeros(33,15,3);
last_cor_seq = zeros(33,15,3);
microonline = zeros(33,15,3);
microoffline = zeros(32,15,3);

% Group loop
for grp=1:length(seq_score(1,1,:)) 
    % Gather subjects within group
    if grp == 1 % REST group
        subs = {'1001','1002','1003','1004','1005','1006','1007','1008', '1009','1010','1011','1012','1013','1014', '1015'}; 
    elseif grp == 2 % ENC group
        subs = {'3001','3002','3003','3004','3005','3006','3007','3008', '3009','3010','3011','3012','3013','3014', '3015'};
    else % SEM Group
        subs = {'4001','4002','4003','4004','4005','4006','4007','4008',  '4009','4010','4011','4012','4013','4014', '4015'}; 
    end
    disp([' -x-x-x-x-x-x-x- Group ' num2str(grp) ' -x-x-x-x-x-x-x- ']);
    % Subject loop
    for subject = 1:size(subs,2)

        % Indicate participant and read in data
        disp(['Processing subject ' num2str(subject) '...']);

        % Read in subject's data
        alldata = MOGS_data(MOGS_data.Subject == str2double(subs{subject}),:);

        % Retrieve the key variables related to task structure and
        % participant performance
        trl_data=alldata.("Block"); % Pull trial data
        resp_data=alldata.("TapSlide1.RESP"); % Pull response data
        iti_data=alldata.("TapSlide1.RTTime"); % Pull clock time data
        RT_data=alldata.("TapSlide1.RT"); % Pull RT times (RT for each slide)
        RetSlide_ACC=alldata.("RetSlide.ACC"); % Pull associative memory test accuracy data
        Retreival=alldata.("Retreival"); % Pull Retrieval trial items
        % Sort associative memory accuracy data by trial experience (1-24).

        % Gather associative memory data for ENC group
        if grp==2 || grp == 3

            % Retrieve the item number assigned to each word pair.
            mem_perf(:,1) = Retreival(~isnan(Retreival)); 

            % Retrieve the response accuracy for each word-pair experienced
            % during the associative memory recall test.
            mem_perf(:,2) = RetSlide_ACC(~isnan(RetSlide_ACC));

            % Retrieve the order that the word-pair occured during
            % acquisition from the beginning to end.
            [sorted_trial, sort_idx] = sort(mem_perf(:,1));

            % Sort the accuracy scores (0 or 1) by the order the word 
            % pairs were experienced during acquisition 
            mem_perf_sorted(:,subject, grp-1) = mem_perf(sort_idx, 2);
        end

        % This variable keeps track of the row in the participants data 
        % where the first tap occurs for every bloctrialk. It is initialized at
        % zero and increments to the first row in every trial
        tap_pre_trl=0;

        % Trial loop (all 33 trials including retention)
        for trl=1:33 

            % Determine the number of keypresses in trial
            tap_trl = size(find(trl_data==trl),1); 

            % Determine the row where the trial starts
            tap_pre_trl = tap_pre_trl+size(find(trl_data==trl-1),1);

            % Initialize keypress total at 1
            tapping=1; 

            % Initialize the count of correct sequences.
            seq_count=0;

            % This while loop will run through all trial taps.
            while (tapping <= tap_trl)

                % Increment the row
                tap = tap_pre_trl+tapping;

                % Identify if the row is not one of the last four taps in
                % the trial
                if(tapping<=tap_trl-4)

                    % If the last five button presses correctly complete
                    % the sequence (4-1-3-2-4)...
                    if (resp_data(tap,1)==4 && resp_data(tap+1,1)==1 && resp_data(tap+2,1)==3 && resp_data(tap+3,1)==2 && resp_data(tap+4,1)==4) 

                        % Increment the count of correct sequences
                        seq_count = seq_count + 1;

                        % Increment number of correct sequences (1 whole sequence)
                        seq_score(subject, trl,grp) = seq_score(subject, trl,grp)+1; 
                        
                        % % Calculate keypresses per second (time between
                        % first button press and last in seconds)
                        ITI(seq_count,trl,subject,grp)= 1/(((iti_data(tap+4,1)-iti_data(tap,1))/4)/1000); 

                        % Jump to the next row that starts a sequence
                        tapping=tapping+5; 

                    % If there was a mistake in the sequence...
                    else 
                        % Move to next row of data
                        tapping=tapping+1; 
                    end

                % If it's one of the final five taps in the trial...
                else 
                    % calculate how many taps are left in the trial
                    tapping_left=tap_trl-tapping+1;

                    % The "switch" here determines whether there are
                    % 4,3,2,or 1 taps left in the trial
                    switch tapping_left 

                        % 4 taps are left in the trial
                        case 4 

                            % If the next four repsonses are correct...
                            if (resp_data(tap,1)==4 && resp_data(tap+1,1)==1 && resp_data(tap+2,1)==3 && resp_data(tap+3,1)==2) 
                                
                                % Increment the count of sequences
                                seq_count = seq_count + 1;

                                % Add 0.8 to the total sequence number
                                seq_score(subject, trl,grp) = seq_score(subject, trl,grp)+0.8; 
    
                                % Calculate keypresses per second for
                                % those four keypresses
                                 ITI(seq_count,trl,subject,grp)=  ( (1/(((iti_data(tap+3,1)-iti_data(tap,1))/3)/1000)));
    
                                % Move to the end of the trial
                                tapping=tapping+4;                          
                            
                            % If any of the next four buttons is incorrect
                            else

                                % Increment the row by 1
                                tapping=tapping+1;
                            end

                        % 3 taps are left in the trial
                        case 3 

                             % If the next three responses are correct...
                            if (resp_data(tap,1)==4 && resp_data(tap+1,1)==1 && resp_data(tap+2,1)==3) 

                                % Increment the count of sequences
                                seq_count = seq_count + 1;

                                 % Increment number of correct sequences
                                % by 0.6
                                seq_score(subject, trl,grp) = seq_score(subject, trl,grp)+0.8; 
                                
                                % Calculate keypresses per second for
                                % those 3 keypresses
                                ITI(seq_count,trl,subject,grp)=  ( (1/(((iti_data(tap+2,1)-iti_data(tap,1))/2)/1000)));
    
                                % Calculate the KPS for those three
                                % keypresses
                                
                                % Move to the end of the trial
                                tapping=tapping+3;

                            % If any of the next three presses is incorrect
                            else

                                % Move to the next row
                                tapping=tapping+1;
                            end

                        % 2 taps are left in the trial
                        case 2 
                    
                            % ...and the next two keypresses are correct...
                            if (resp_data(tap,1)==4 && resp_data(tap+1,1)==1) % If the next two responses are correct

                                 % Increment the count of sequences
                                seq_count = seq_count + 1;
                                
                                % Increment the number of correct
                                % sequences by 0.4.
                                seq_score(subject, trl,grp) = seq_score(subject, trl,grp)+0.4; 
                                
                                % Calculate KPS for those two trials
                                ITI(seq_count,trl,subject,grp)=  (1/(((iti_data(tap+1,1)-iti_data(tap,1))/1)/1000));

                                % Go to the end of the trial
                                tapping=tapping+2;

                            % Only 1 keypress left
                            else

                                %Move to the end of the trial
                                tapping=tapping+1; 
                            end

                        case 1 % 1 tap is left in the trial

                            % Move to the final tap
                            tapping=tapping+1;

                    % End of number of taps that are left (switch)
                    end 
                
                % End of determining if there are enough taps left to make a full sequence (if)
                end 
    
             % End of trial data    
            end        

             % Check to see if there are any correct sequences in the
             % trial and if not, fill in the data with NaNs
            try 
                if nnz((ITI(:,trl,subject,grp))) == 0
    
                    % Fill in the value with NaNs
                    ITI(:,trl,subject,grp)= NaN; 
                    
                end
            catch
                ITI(:,trl,subject,grp)= NaN; 
            end

        % Done gathering performance metrics for the trial
        end

        % Check that all trials have a first and last trial that is
        % accurate
        for trl=1:33

            % Check that the first sequence is correct
            if ITI(1,trl,subject,grp) ~= 0

                % If so, get the KPS for that sequence
                first_cor_seq(trl,subject,grp)=ITI(1,trl,subject,grp);
            else

                % If not, throw an error signal
                disp(['error in first mogs. 0 in ITI trial = ', num2str(trl), ', sub=', num2str(subject)]);
            end

            % Check for last mogs
            if ITI(nnz((ITI(:,trl,subject,grp))),trl,subject,grp) ~= 0 
                
                % If there is a correct sequence
                last_cor_seq(trl,subject,grp)=ITI(nnz((ITI(:,trl,subject,grp))),trl,subject);
            else
                
                % If not, throw an error signal
                disp(['error in first mogs. 0 in ITI trial = ', num2str(trl), ', sub=', num2str(subject)]);
            end

        end

        % Calculate microonline gains
        for trl=1:33

            % Micro-online gains are the difference between first and last
            % sequence performance within a trial. Zeroes should be treated
            % as NaNs
            microonline(trl,subject,grp)=last_cor_seq(trl,subject,grp)-first_cor_seq(trl,subject,grp); 
        end

        % Calculate microoffline gains
        for trl=1:32

            % Micro-offline gains are the difference between first and last
            % sequence performance between consecutive trials
            microoffline(trl,subject,grp)=first_cor_seq(trl+1,subject,grp) - last_cor_seq(trl,subject,grp);
        end
    end % Subject
    clearvars -except grp microoffline microonline ITI mem_perf mem_perf_sorted home_dir out_dir root_dir seq_score seq_ccount MOGS_data
% Done calculating all performance metrics for all groups
end 

%% Saving and finishing up
cd(out_dir);
% Write out data

% Number of correct sequences by participant, trial, and group
save('Number_of_correct_sequences.mat',"seq_score","-v7.3") 

% Micro-online gains by trial, participant, and group
save('Microonline_data.mat',"microonline","-v7.3")

% Micro-offline gains by trial, participant, and group
save('Microoffline_data.mat',"microoffline","-v7.3")

% Memory performance by trial, participant, and group
save('ENC_SEM_Memory_Performance.mat',"mem_perf_sorted","-v7.3")

% Sequence performance time by sequence, trial, participant, and group
save("ITI.mat","ITI","-v7.3")

cd(root_dir);

disp(' - - - - - - - - - - - - - - - - - - - - - ');
disp('        Finished calculating MOGS!         ');
disp(' - - - - - - - - - - - - - - - - - - - - - ');