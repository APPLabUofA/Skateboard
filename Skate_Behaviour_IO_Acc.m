clear all
close all
ccc
%%
exp = 'Skateboard';
subs = {'100'	'101'	'102'	'103'	'104'	'106'	'107'	'108'	'109'	'110'...
    '111'	'112'	'113'	'114'	'115'	'116'	'117'	'119'	'120'...
    '122'	'123'	'125'	'126'	'127'	'129'	'130'...
    '131'	'132'	'133'	'134'	'135'	'136'	'137'};
is_goofy = [0,	0,	0,	0,	0,	0,	1,	0,	1,	1,	0,	0,	1,	1,...
    1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	0,	1,	1,    0,	1,	0];
dif_trig = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
            0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

%%
nsubs = length(subs);
conds = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
new_cond = {'P_In'; 'P_Out'; 'NP_In'; 'NP_Out'};

%preferred, clockwise - non-preffered, CCW
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\winter2019\';

% if ~exist([Pathname 'segments\'])
%     mkdir([Pathname 'segments\']);
% end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Marker Numbers
nStandard = 3;
nTarget = 5;
nFalseAlarm = 7;
nCorrectResponse = 9;

prop_correct = zeros(nsubs,length(new_cond));
prop_correctRej = zeros(nsubs,length(new_cond));
medianACC_correct = zeros(nsubs,length(new_cond));
medianRT_correct = zeros(nsubs,length(new_cond));
medianRT_falseAlarm = zeros(nsubs,length(new_cond));
%%
for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond} '.vhdr'];
        setname = Filename(1:end-5);
        
        EEG = pop_loadbv(Pathname, Filename);
        
        %% Find all event types and latencys
        
        %string trigger from brain recorder  (e.g. 'S  1')
        event_strings = {EEG.event.type}; %array of marker label strings
        %time since recording start in ms (integer)
        event_latency = [EEG.event.latency];
                
        %convert strings to integers so they are easier to work with
        event_markers = zeros(size(event_strings));
        event_markers(strcmp(event_strings,'S  3')) = nStandard;   %standard
        event_markers(strcmp(event_strings,'S  5')) = nTarget;   %target
        event_markers(strcmp(event_strings,'S  7')) = nFalseAlarm;   %false alarm
        event_markers(strcmp(event_strings,'S  9')) = nCorrectResponse;   %correct response
        
        if dif_trig(i_sub)
            nStandard = 1;
            nTarget = 2;
            nFalseAlarm = 3;
            nCorrectResponse = 4;
            
            event_markers(strcmp(event_strings,'S  1')) = nStandard;   %standard
            event_markers(strcmp(event_strings,'S  2')) = nTarget;   %target
            event_markers(strcmp(event_strings,'S  3')) = nFalseAlarm;   %false alarm
            event_markers(strcmp(event_strings,'S  4')) = nCorrectResponse;   %correct response
        end
        
        event_latency(event_markers == 0) = []; %remove any extra triggers
        event_markers(event_markers == 0) = []; % remove any extra triggers
        
        %% now step through the arrays and check for stuff
        
        %setup counters
        count_tones = 0;
        count_targets = 0;
        count_standards = 0;
        
        count_correct = 0;
        count_misses = 0;
        count_correctRej = 0;
        count_falseAlarm = 0;
        
        ACC_correct = [];
%         RT_correct = [];
%         RT_falseAlarm = [];
        
        %for every event
        for i_event = 1:length(event_markers)-1 %last one is a filler markers
            
            this_marker = event_markers(i_event);
            tone_time = event_latency(i_event);
            next_marker = event_markers(i_event+1);
            next_time = event_latency(i_event+1);
            potential_ACC = count_correct/50*100;
%             potential_RT = next_time-tone_time;
            
            %if it is a tone (|| means or)
            if this_marker == nTarget || this_marker == nStandard
                count_tones = count_tones + 1;
                fprintf('\n Tone Number: ') %\n is a new line
                fprintf(num2str(count_tones))
                fprintf(' --> ')
                
                fprintf('This marker: ')
                fprintf(num2str(this_marker))
                fprintf(' , ')
                
                fprintf('Next marker: ')
                fprintf(num2str(next_marker))
                fprintf(' , ')
                
                
            end
            
            if this_marker == nTarget 
                count_targets = count_targets + 1;
                
                
                %if correct response
                if next_marker == nCorrectResponse 
                    count_correct = count_correct + 1;
                    ACC_correct = [ACC_correct potential_ACC];
                    fprintf('Responded -- > Cummulative Accuracy = ')
                    fprintf(num2str(potential_ACC))
                    fprintf(' %')
                    
                    %if miss since next is another tone
                elseif next_marker == nStandard || next_marker == nTarget
                    count_misses = count_misses + 1;
                    fprintf('Did not respond')
                    
                    %anything else?
                else
                    fprintf('Not 9 or 3 or 5 (or 4, 1 or 2)')
                    
                end
                
            elseif this_marker == nStandard
                count_standards = count_standards + 1;
                
                %if correct rejection since next is another tone
                if next_marker == nStandard || next_marker == nTarget
                    count_correctRej = count_correctRej + 1;
                    fprintf('Correct Rejection')
                    
                    %if false alarm
                elseif next_marker == nFalseAlarm
                    %RT_falseAlarm = [RT_falseAlarm potential_RT];
                    count_falseAlarm = count_falseAlarm + 1;
                    fprintf('False Alarm')

                    
                    %anything else?
                else
                    fprintf('Not 7 or 3 or 5 (or 3 1 2')
                end
                
            end %if target,elseStandard
        end %every event
        
        %reassign conditions
        %IF GOOFY: P_CCW = preferred_in (new_cond_index = 1;)
        %IF GOOFY: P_CW = preferred_out (new_cond_index = 2;)
        %IF GOOFY: NP_CW = nonpreferred_in (new_cond_index = 3;)
        %IF GOOFY: NP_CCW = nonpreferred_out (new_cond_index = 4;)
        %IF regular: P_CW = preferred_in (new_cond_index = 1;)
        %IF regular: P_CCW = preferred_out (new_cond_index = 2;)
        %IF regular: NP_CCW = nonpreferred_in (new_cond_index = 3;)
        %IF regular: NP_CW = nonpreferred_out (new_cond_index = 4;)
        if is_goofy(i_sub)
            if strcmp(conds{i_cond},'P_CCW')
                new_cond_index = 1;
            elseif strcmp(conds{i_cond},'P_CW')
                new_cond_index = 2;
            elseif strcmp(conds{i_cond},'NP_CW')
                new_cond_index = 3;
            elseif strcmp(conds{i_cond},'NP_CCW')
                new_cond_index = 4;
            end
        elseif ~is_goofy(i_sub)
            if strcmp(conds{i_cond},'NP_CW')
                new_cond_index = 4;
            elseif strcmp(conds{i_cond},'NP_CCW')
                new_cond_index = 3;
            elseif strcmp(conds{i_cond},'P_CCW')
                new_cond_index = 2;
            elseif strcmp(conds{i_cond},'P_CW')
                new_cond_index = 1;
            end
        end
        
        prop_correct(i_sub,new_cond_index) = count_correct / count_targets;
        prop_correctRej(i_sub,new_cond_index) = count_correctRej / count_standards;
        medianACC_correct(i_sub,new_cond_index) = median(potential_ACC);
        %medianRT_falseAlarm(i_sub,new_cond_index) = median(RT_falseAlarm);
        
    end
end

%%
n_conditions = size(medianACC_correct,2);
%this the grand mean over subjects of the median RTs
grand_mean_ACC_Corr = mean(medianACC_correct);
%these are normal error bars (Standard Error)
grand_SE_ACC_Corr = std(medianACC_correct)/sqrt(nsubs);
%these are smaller within subject error bars
%made by subtracting away each subjects average
%from their other scores to remove between subject difference
sub_mean_ACC_Corr = mean(medianACC_correct,2); %average for each subject
%this subtracts each subjects average from their scores
%repmat repeats the matrix 4 times for each condition
mean_ACC_Corr_deviation = medianACC_correct - repmat(sub_mean_ACC_Corr,1,n_conditions);

%calculating standard errors for plots
SE_P = nanstd(mean(medianACC_correct(:,1:2),2))/sqrt(nsubs); 
SE_NP = nanstd(mean(medianACC_correct(:,3:4),2))/sqrt(nsubs);

%then take the standard error of those deviatoins from the mean
grand_withinSE_ACC_Corr = std(mean_ACC_Corr_deviation)/sqrt(nsubs);

%now do the same for proportion correct
grand_mean_prop_corr = nanmean(prop_correct);
grand_SE_prop_corr = nanstd(prop_correct)/sqrt(nsubs);
sub_mean_prop_corr = nanmean(prop_correct,2);
prop_corr_deviation = prop_correct - repmat(sub_mean_prop_corr,1,n_conditions);
grand_withinSE_prop_corr = nanstd(prop_corr_deviation)/sqrt(nsubs);

%% original 4 conditions plus proportion correct
% %plot it
conds_plot = {'P_FaceIn'; 'P_FaceOut';'NP_FaceIn'; 'NP_FaceOut'}; 
figure;
subplot (1,2,1)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(grand_mean_ACC_Corr,grand_withinSE_ACC_Corr);
ylim([50 110])
ylabel('Mean Accuracy')
title('Target Accuracy')
legend(conds)
subplot(1,2,2)
barweb(grand_mean_prop_corr,grand_withinSE_prop_corr);
ylim([.9 1])
ylabel('Proportion')
title('Proportion of Targets responded to')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STATS for all 4 conditions 
[p,tbl,stats] = anova1(medianACC_correct);

%%       

%side by side plots - PREFERENCE on X Axis
grand_meanACC_Pref = grand_mean_ACC_Corr(1:2);
grand_meanACC_NPref = grand_mean_ACC_Corr(3:4);
grand_wSE_ACC_Pref = grand_withinSE_ACC_Corr(1:2);
grand_wSE_ACC_NPref = grand_withinSE_ACC_Corr(3:4);

%side by side plots - FACING on X Axis
gran_meanACC_Face_IN = grand_mean_ACC_Corr(1:2:3);
gran_meanACC_Face_OUT = grand_mean_ACC_Corr(2:2:4);
grand_wSE_ACC_Face_IN = grand_withinSE_ACC_Corr(1:2:3);
grand_wSE_ACC_Face_OUT = grand_withinSE_ACC_Corr(2:2:4);

%side by side plots - PREFERENCE on X Axis
conds_plot = {'Facing In'; 'Facing Out'}; 
figure;
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
subplot(1,2,1)
barweb(grand_meanACC_Pref,grand_wSE_ACC_Pref);
ylim([80 100])
ylabel('Accuracy %')
xlabel ('Preferred')
title('Target Accuracy By Stance')
legend(conds_plot)
subplot(1,2,2)
conds_plot = {'Facing In'; 'Facing Out'}; 
barweb(grand_meanACC_NPref,grand_wSE_ACC_NPref);
ylim([80 100])
ylabel('Accuracy %')
xlabel('Non-preferred')
title('Target Accuracy By Stance')
legend(conds_plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATS for PREFERRED by stance
[h p ci stat] = ttest(medianACC_correct(:,1),medianACC_correct(:,2),.05,'both',1) 
%Stats for non-preferred by stance
[h p ci stat] = ttest(medianACC_correct(:,3),medianACC_correct(:,4),.05,'both',1) 


%side by side plots - FACING on X Axis
conds_plot = {'Preferred Stance'; 'Non-Preferred Stance'}; 
figure;
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
subplot(1,2,1)
barweb(gran_meanACC_Face_IN,grand_wSE_ACC_Face_IN);
ylim([80 100])
ylabel('Accuracy %')
xlabel ('Facing In')
title('Target Accuracy By Face Orientation')
legend(conds_plot)
subplot(1,2,2)
conds_plot = {'Preferred Stance'; 'Non-Preferred Stance'}; 
barweb(gran_meanACC_Face_OUT,grand_wSE_ACC_Face_OUT);
ylim([80 100])
ylabel('Accuracy %')
xlabel ('Facing Out')
title('Target Accuracy By Face Orientation')
legend(conds_plot)

%Stats for facing inside (pref vs non-pref)
[h p ci stat] = ttest (medianACC_correct(:,1),medianACC_correct(:,3),.05,'both',1) 

gran_meanACC_Face_OUT = grand_mean_ACC_Corr(2:2:4);
%Stats for facing inside (pref vs non-pref)
[h p ci stat] = ttest(medianACC_correct(:,2),medianACC_correct(:,4),.05,'both',1) 


%% global stance!!!!!!!!!!!!!!!!!!
%making a global comparison of preferred vs non-preferred 
global_meanACC_Pref = mean (grand_mean_ACC_Corr(1:2));
global_meanACC_NPref = mean(grand_mean_ACC_Corr(3:4));
global_mean_PNP = [global_meanACC_Pref; global_meanACC_NPref];

%standard errors
SE_P = nanstd(mean(medianACC_correct(:,1:2),2))/sqrt(nsubs); 
SE_NP = nanstd(mean(medianACC_correct(:,3:4),2))/sqrt(nsubs);

global_SE_PNP = [SE_P; SE_NP];

%ENSURE THAT WITHIN STANDARD ERROR WAS CALCULATED PROPERLY
figure;
conds_plot = {'Preferred Stance'; 'Non-Preferred Stance'}; 
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(global_mean_PNP,global_SE_PNP);
ylim([80 100])
ylabel('Accuracy %')
xlabel ('Stance')
title('Target Accuracy by Overall Preference')
legend(conds_plot)

%t-test
[h p ci stat] = ttest(mean(medianACC_correct(:,1:2),2),mean(medianACC_correct(:,3:4),2),.05,'both',1) 



%% 
%GLOBAL FACING ORIENTATION
global_meanACC_Face_IN = mean(grand_mean_ACC_Corr(1:2:3));
global_meanACC_Face_OUT = mean(grand_mean_ACC_Corr(2:2:4));
global_mean_IO = [global_meanACC_Face_IN; global_meanACC_Face_OUT];

SE_IN = nanstd(mean(medianACC_correct(:,1:2:3),2)/sqrt(nsubs)); 
SE_OUT = nanstd(mean(medianACC_correct(:,2:2:4),2)/sqrt(nsubs));
% global_wSE_ACC_Pref = mean (grand_withinSE_ACC_Corr(1:2));
% global_wSE_ACC_NPref = mean (grand_withinSE_ACC_Corr(3:4));
global_SE_IO = [SE_IN; SE_OUT];

%ENSURE THAT WITHIN STANDARD ERROR WAS CALCULATED PROPERLY
figure;
conds_plot = {'Facing In'; 'Facing Out'}; 
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(global_mean_IO,global_SE_IO);
ylim([80 100])
ylabel('Accuracy %')
xlabel ('Facing Orientation')
title('Target Accuracy by Overall Facing Orientation')
legend(conds_plot)

[h p ci stat] = ttest(mean(medianACC_correct(:,1:2:3),2),mean(medianACC_correct(:,2:2:4),2),.05,'both',1)

