clear all
close all
ccc
%%
exp = 'Skateboard';
subs = {'100' '101' '102' '103' '104' '106' '107' '108' '109' '110' '111'...
    '112' '113' '114' '115' '116' '117' '118' '119' '120' '122' '123' '124' '125'...
    '126' '127' '128' '129'};
is_goofy = [0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1,...
            0, 0, 0, 0, 0, 1, 1];
dif_trig = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
            0, 1, 1, 1, 1, 1, 1];
% 107 - Goofy
% 109 - Goofy
% 110 - Goofy
% 113 - Goofy
% 114 - Goofy
% 115 - Goofy
% 122 - Goofy
% 128 - Goofy
% 129 - Goofy
%%
nsubs = length(subs);
conds = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
new_cond = {'P_In'; 'P_Out'; 'NP_In'; 'NP_Out'};

%preferred, clockwise - non-preffered, CCW
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\winter2019\';

if ~exist([Pathname 'segments\'])
    mkdir([Pathname 'segments\']);
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Marker Numbers
nStandard = 3;
nTarget = 5;
nFalseAlarm = 7;
nCorrectResponse = 9;

prop_correct = zeros(nsubs,length(new_cond));
prop_correctRej = zeros(nsubs,length(new_cond));
medianACC_correct = zeros(nsubs,length(new_cond));
% medianRT_correct = zeros(nsubs,length(new_cond));
% medianRT_falseAlarm = zeros(nsubs,length(new_cond));
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
                    fprintf('Responded -- > Accuracy total = ')
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
        medianACC_correct(i_sub,new_cond_index) = median(ACC_correct);
        %medianRT_falseAlarm(i_sub,new_cond_index) = median(RT_falseAlarm);
        
    end
end

n_conditions = size(medianACC_correct,2)
%this the grand mean over subjects of the median RTs
grand_mean_ACC_Corr = mean(medianACC_correct)
%these are normal error bars (Standard Error)
grand_SE_ACC_Corr = std(medianACC_correct)/sqrt(nsubs)
%these are smaller within subject error bars
%made by subtracting away each subjects average
%from their other scores to remove between subject difference
sub_mean_ACC_Corr = mean(medianACC_correct,2); %average for each subject
%this subtracts each subjects average from their scores
%repmat repeats the matrix 4 times for each condition
mean_ACC_Corr_deviation = medianACC_correct - repmat(sub_mean_ACC_Corr,1,n_conditions);
%then take the standard error of those deviatoins from the mean
grand_withinSE_ACC_Corr = std(mean_ACC_Corr_deviation)/sqrt(nsubs)


%now do the same for proportion correct
grand_mean_prop_corr = mean(prop_correct);
grand_SE_prop_corr = std(prop_correct)/sqrt(nsubs);
sub_mean_prop_corr = mean(prop_correct,2);
prop_corr_deviation = prop_correct - repmat(sub_mean_prop_corr,1,n_conditions);
grand_withinSE_prop_corr = std(prop_corr_deviation)/sqrt(nsubs);

%plot it
conds_plot = {'Pref_FaceIn'; 'Pref_FaceOut';'NonPref_FaceIn'; 'NonPref_FaceOut'}; 
figure;
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
subplot(1,2,1)
barweb(grand_mean_RT_Corr,grand_withinSE_ACC_Corr);
ylim([450 525])
ylabel('Median RT (ms)')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)
subplot(1,2,2)
barweb(grand_mean_prop_corr,grand_withinSE_prop_corr);
ylim([.9 1])
ylabel('Proportion')
title('Proportion of Targets responded to')

%side by side plots - PREFERENCE on X Axis
gran_meanRT_FaceIn = grand_mean_RT_Corr(1:2);
grand_meanRT_FaceOut = grand_mean_RT_Corr(3:4);
grand_wSE_RT_FaceIn = grand_withinSE_ACC_Corr(1:2);
grand_wSE_RT_FaceOut = grand_withinSE_ACC_Corr(3:4);

conds_plot = {'FacingIn'; 'FacingOut'}; 
figure;
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
subplot(2,2,1)
barweb(gran_meanRT_FaceIn,grand_wSE_RT_FaceIn);
ylim([450 525])
ylabel('Median RT (ms)')
xlabel ('Preferred Stance')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)
subplot(2,2,2)
conds_plot = {'FacingIn'; 'FacingOut'}; 
barweb(grand_meanRT_FaceOut,grand_wSE_RT_FaceOut);
ylim([450 525])
ylabel('Median RT (ms)')
xlabel('Non-preferred Stance')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)
subplot(2,2,3)
barweb(grand_mean_prop_corr(1:2),grand_withinSE_prop_corr(1:2));
ylim([.9 1])
ylabel('Proportion')
xlabel ('Preferred Stance')
title('Proportion of Targets responded to')
subplot(2,2,4)
barweb(grand_mean_prop_corr(3:4),grand_withinSE_prop_corr(3:4));
ylim([.9 1])
ylabel('Proportion')
xlabel ('Non-preferred Stance')
title('Proportion of Targets responded to')

%side by side plots - FACING on X Axis
gran_meanRT_P_NP_IN = grand_mean_RT_Corr(1:2:3);
gran_meanRT_P_NP_OUT = grand_mean_RT_Corr(2:2:4);
grand_wSE_RT_P_NP_IN = grand_withinSE_ACC_Corr(1:2:3);
grand_wSE_RT_P_NP_OUT = grand_withinSE_ACC_Corr(2:2:4);

conds_plot = {'Preferred'; 'Non-Preferred'}; 
figure;
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
subplot(1,2,1)
barweb(gran_meanRT_P_NP_IN,grand_wSE_RT_P_NP_IN);
ylim([450 525])
ylabel('Median RT (ms)')
xlabel ('Facing In')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)
subplot(1,2,2)
conds_plot = {'Preferred'; 'Non-Preferred'}; 
barweb(gran_meanRT_P_NP_OUT,grand_wSE_RT_P_NP_OUT);
ylim([450 525])
ylabel('Median RT (ms)')
xlabel ('Facing Out')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)


