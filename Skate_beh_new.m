clear all
close all
ccc
%%
expacc = 'Skateboard';
% subsacc = {'100'	'101'	'102'	'103'	'104'	'106'	'107'	'108'	'109'	'110'...
%     '111'	'112'	'113'	'114'	'115'	'116'	'117'	'119'	'120'...
%     '122'	'123'	'125'	'126'	'127'	'129'	'130'...
%     '131'	'132'	'133'	'134'	'135'	'136'	'137'};
% is_goofyacc = [0,	0,	0,	0,	0,	0,	1,	0,	1,	1,	0,	0,	1,	1,...
%     1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	0,	1,	1,    0,	1,	0];
% dif_trigacc = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
%             0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
subsacc = {'100'	'101'	'102'	'103'	'104'	'106'	'108'	'109'	'110'...
    '111'	'112'	'113'	'115'	'116'	'117'	'119'	'120'...
    '122'	'123'	'125'	'126'	'127'	'129'	...
    '131'	'132'	'133'	'134'	'136'	'137'};
is_goofyacc = [0,	0,	0,	0,	0,	0,	0,	1,	1,	0,	0,	1,	...
    1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	1,    1,	0];        
dif_trigacc = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
            1, 1, 1, 1, 1, 1, 1, 1, 1,1];
%%
nsubsacc = length(subsacc);
condsacc = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
new_condacc = {'P_In'; 'P_Out'; 'NP_In'; 'NP_Out'};

%preferred, clockwise - non-preffered, CCW
ncondsacc = length(condsacc);
Pathnameacc = 'M:\Data\Skateboard\winter2019\';

% if ~exist([Pathname 'segments\'])
%     mkdir([Pathname 'segments\']);
% end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Marker Numbers
nStandardacc = 3;
nTargetacc = 5;
nFalseAlarmacc = 7;
nCorrectResponseacc = 9;

prop_correctacc = zeros(nsubsacc,length(new_condacc));
prop_correctRejacc = zeros(nsubsacc,length(new_condacc));
medianACC_correct = zeros(nsubsacc,length(new_condacc));
%medianRT_correct = zeros(nsubsacc,length(new_condacc));
%medianRT_falseAlarm = zeros(nsubsacc,length(new_condacc));
%%
for i_sub = 1:nsubsacc
    for i_cond = 1:ncondsacc
        
        Filenameacc = [subsacc{i_sub} '_' expacc '_' condsacc{i_cond} '.vhdr'];
        setnameacc = Filenameacc(1:end-5);
        
        EEG = pop_loadbv(Pathnameacc, Filenameacc);
        
        %% Find all event types and latencys
        
        %string trigger from brain recorder  (e.g. 'S  1')
        event_stringsacc = {EEG.event.type}; %array of marker label strings
        %time since recording start in ms (integer)
        event_latencyacc = [EEG.event.latency];
                
        %convert strings to integers so they are easier to work with
        event_markersacc = zeros(size(event_stringsacc));
        event_markersacc(strcmp(event_stringsacc,'S  3')) = nStandardacc;   %standard
        event_markersacc(strcmp(event_stringsacc,'S  5')) = nTargetacc;   %target
        event_markersacc(strcmp(event_stringsacc,'S  7')) = nFalseAlarmacc;   %false alarm
        event_markersacc(strcmp(event_stringsacc,'S  9')) = nCorrectResponseacc;   %correct response
        
        if dif_trigacc(i_sub)
            nStandardacc = 1;
            nTargetacc = 2;
            nFalseAlarmacc = 3;
            nCorrectResponseacc = 4;
            
            event_markersacc(strcmp(event_stringsacc,'S  1')) = nStandardacc;   %standard
            event_markersacc(strcmp(event_stringsacc,'S  2')) = nTargetacc;   %target
            event_markersacc(strcmp(event_stringsacc,'S  3')) = nFalseAlarmacc;   %false alarm
            event_markersacc(strcmp(event_stringsacc,'S  4')) = nCorrectResponseacc;   %correct response
        end
        
        event_latencyacc(event_markersacc == 0) = []; %remove any extra triggers
        event_markersacc(event_markersacc == 0) = []; % remove any extra triggers
        
        %% now step through the arrays and check for stuff
        
        %setup counters
        count_tonesacc = 0;
        count_targetsacc = 0;
        count_standardsacc = 0;
        
        count_correctacc = 0;
        count_missesacc = 0;
        count_correctRejacc = 0;
        count_falseAlarmacc = 0;
        
        ACC_correct = [];
%         RT_correct = [];
%         RT_falseAlarm = [];
        
        %for every event
        for i_event = 1:length(event_markersacc)-1 %last one is a filler markers
            
            this_markeracc = event_markersacc(i_event);
            tone_timeacc = event_latencyacc(i_event);
            next_markeracc = event_markersacc(i_event+1);
            next_timeacc = event_latencyacc(i_event+1);
            potential_ACC = count_correctacc/50*100;
%             potential_RT = next_time-tone_time;
            
            %if it is a tone (|| means or)
            if this_markeracc == nTargetacc || this_markeracc == nStandardacc
                count_tonesacc = count_tonesacc + 1;
                fprintf('\n Tone Number: ') %\n is a new line
                fprintf(num2str(count_tonesacc))
                fprintf(' --> ')
                
                fprintf('This marker: ')
                fprintf(num2str(this_markeracc))
                fprintf(' , ')
                
                fprintf('Next marker: ')
                fprintf(num2str(next_markeracc))
                fprintf(' , ')
                
                
            end
            
            if this_markeracc == nTargetacc 
                count_targetsacc = count_targetsacc + 1;
                
                
                %if correct response
                if next_markeracc == nCorrectResponseacc 
                    count_correctacc = count_correctacc + 1;
                    ACC_correct = [ACC_correct potential_ACC];
                    fprintf('Responded -- > Cummulative Accuracy = ')
                    fprintf(num2str(potential_ACC))
                    fprintf(' %')
                    
                    %if miss since next is another tone
                elseif next_markeracc == nStandardacc || next_markeracc == nTargetacc
                    count_missesacc = count_missesacc + 1;
                    fprintf('Did not respond')
                    
                    %anything else?
                else
                    fprintf('Not 9 or 3 or 5 (or 4, 1 or 2)')
                    
                end
                
            elseif this_markeracc == nStandardacc
                count_standardsacc = count_standardsacc + 1;
                
                %if correct rejection since next is another tone
                if next_markeracc == nStandardacc || next_markeracc == nTargetacc
                    count_correctRejacc = count_correctRejacc + 1;
                    fprintf('Correct Rejection')
                    
                    %if false alarm
                elseif next_markeracc == nFalseAlarmacc
                    %RT_falseAlarm = [RT_falseAlarm potential_RT];
                    count_falseAlarmacc = count_falseAlarmacc + 1;
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
        if is_goofyacc(i_sub)
            if strcmp(condsacc{i_cond},'P_CCW')
                new_cond_indexacc = 1;
            elseif strcmp(condsacc{i_cond},'P_CW')
                new_cond_indexacc = 2;
            elseif strcmp(condsacc{i_cond},'NP_CW')
                new_cond_indexacc = 3;
            elseif strcmp(condsacc{i_cond},'NP_CCW')
                new_cond_indexacc = 4;
            end
        elseif ~is_goofyacc(i_sub)
            if strcmp(condsacc{i_cond},'NP_CW')
                new_cond_indexacc = 4;
            elseif strcmp(condsacc{i_cond},'NP_CCW')
                new_cond_indexacc = 3;
            elseif strcmp(condsacc{i_cond},'P_CCW')
                new_cond_indexacc = 2;
            elseif strcmp(condsacc{i_cond},'P_CW')
                new_cond_indexacc = 1;
            end
        end
        
        prop_correctacc(i_sub,new_cond_indexacc) = count_correctacc / count_targetsacc;
        prop_correctRejacc(i_sub,new_cond_indexacc) = count_correctRejacc / count_standardsacc;
        medianACC_correct(i_sub,new_cond_indexacc) = median(potential_ACC);
        %medianRT_falseAlarm(i_sub,new_cond_index) = median(RT_falseAlarm);
        
    end
end

%%
n_conditionsacc = size(medianACC_correct,2);
%this the grand mean over subjects of the median RTs
grand_mean_ACC_Corr = mean(medianACC_correct);
%these are normal error bars (Standard Error)
grand_SE_ACC_Corr = std(medianACC_correct)/sqrt(nsubsacc);
%these are smaller within subject error bars
%made by subtracting away each subjects average
%from their other scores to remove between subject difference
sub_mean_ACC_Corr = mean(medianACC_correct,2); %average for each subject
%this subtracts each subjects average from their scores
%repmat repeats the matrix 4 times for each condition
mean_ACC_Corr_deviation = medianACC_correct - repmat(sub_mean_ACC_Corr,1,n_conditionsacc);

%calculating standard errors for plots
SE_Pacc = nanstd(mean(medianACC_correct(:,1:2),2))/sqrt(nsubsacc); 
SE_NPacc = nanstd(mean(medianACC_correct(:,3:4),2))/sqrt(nsubsacc);

%then take the standard error of those deviatoins from the mean
grand_withinSE_ACC_Corr = std(mean_ACC_Corr_deviation)/sqrt(nsubsacc);

%now do the same for proportion correct
grand_mean_prop_corracc = nanmean(prop_correctacc);
grand_SE_prop_corracc = nanstd(prop_correctacc)/sqrt(nsubsacc);
sub_mean_prop_corracc = nanmean(prop_correctacc,2);
prop_corr_deviationacc = prop_correctacc - repmat(sub_mean_prop_corracc,1,n_conditionsacc);
grand_withinSE_prop_corr = nanstd(prop_corr_deviationacc)/sqrt(nsubsacc);

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
legend(condsacc)
subplot(1,2,2)
barweb(grand_mean_prop_corracc,grand_withinSE_prop_corr);
ylim([.9 1])
ylabel('Proportion')
title('Proportion of Targets responded to')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STATS for all 4 conditions 
[p,tbl,stats] = anova1(medianACC_correct);

%repeated measures REAL STATS
s1 = 1:1:length(medianACC_correct(:,1));
s2 = 1:1:length(medianACC_correct(:,2));
s3 = 1:1:length(medianACC_correct(:,3));
s4 = 1:1:length(medianACC_correct(:,4));

cond1 = ones(length(s1),1)+0;
cond2 = ones(length(s2),1)+1;
cond3 = ones(length(s3),1)+2;
cond4 = ones(length(s4),1)+3;

Y = [medianACC_correct(:,1);medianACC_correct(:,2);medianACC_correct(:,3);medianACC_correct(:,4)];
S = [s1';s2';s3';s4'];
F1 = [cond1;cond2;cond3;cond4];
FACTNAMES = {'CONDITION'};
data_anova = [Y, F1, S];
stats = rmanova2(data_anova, 0.05,1,1);
RMAOV1(data_anova)







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
global_meanACC_PNP = [global_meanACC_Pref; global_meanACC_NPref];

%standard errors
SE_Pacc = nanstd(mean(medianACC_correct(:,1:2),2))/sqrt(nsubsacc); 
SE_NPacc = nanstd(mean(medianACC_correct(:,3:4),2))/sqrt(nsubsacc);

global_SEacc_PNP = [SE_Pacc; SE_NPacc];
%ENSURE THAT WITHIN STANDARD ERROR WAS CALCULATED PROPERLY
figure;
conds_plot = {'Preferred Stance'; 'Non-Preferred Stance'}; 
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(global_meanACC_PNP,global_SEacc_PNP);
ylim([80 100])
ylabel('Accuracy %')
xlabel ('Stance')
title('Accuracy')
legend(conds_plot)

%t-test
[h p ci stat] = ttest(mean(medianACC_correct(:,1:2),2),mean(medianACC_correct(:,3:4),2),.05,'both',1) 
difference_acc = mean(medianACC_correct(:,1:2),2)- mean(medianACC_correct(:,3:4),2)
[h p ci stat] = ttest(difference_acc,0,.05,'right',1) 
mean (difference_acc)



%%

%%
exp = 'Skateboard';
subs = {'100'	'101'	'102'	'103'	'104'	'106'	'108'	'109'	'110'...
    '111'	'112'	'113'	'115'	'116'	'117'	'119'	'120'...
    '122'	'123'	'125'	'126'	'127'	'129'	...
    '131'	'132'	'133'	'134'	'136'	'137'};
is_goofy = [0,	0,	0,	0,	0,	0,	0,	1,	1,	0,	0,	1,	...
    1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	1,    1,	0];        
dif_trig = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
            1, 1, 1, 1, 1, 1, 1, 1, 1,1];
       

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
new_conds = {'P_In'; 'P_Out'; 'NP_In'; 'NP_Out'};

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

prop_correct = zeros(nsubs,length(new_conds));
prop_correctRej = zeros(nsubs,length(new_conds));
medianRT_correct = zeros(nsubs,length(new_conds));
medianRT_falseAlarm = zeros(nsubs,length(new_conds));

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
        
        RT_correct = [];
        RT_falseAlarm = [];
        
        %for every event
        for i_event = 1:length(event_markers)-1 %last one is a filler markers
            
            this_marker = event_markers(i_event);
            tone_time = event_latency(i_event);
            next_marker = event_markers(i_event+1);
            next_time = event_latency(i_event+1);
            potential_RT = next_time-tone_time;
            
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
                    RT_correct = [RT_correct potential_RT];
                    fprintf('Responded -- > RT = ')
                    fprintf(num2str(potential_RT))
                    fprintf(' ms')
                    
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
                    RT_falseAlarm = [RT_falseAlarm potential_RT];
                    count_falseAlarm = count_falseAlarm + 1;
                    fprintf('False Alarm -- > RT = ')
                    fprintf(num2str(potential_RT))
                    fprintf(' ms')
                    
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
        medianRT_correct(i_sub,new_cond_index) = median(RT_correct);
        medianRT_falseAlarm(i_sub,new_cond_index) = median(RT_falseAlarm);
        
    end
end

%%

n_conditions = size(medianRT_correct,2);
%this the grand mean over subjects of the median RTs
grand_mean_RT_Corr = mean(medianRT_correct);
%these are normal error bars (Standard Error)
grand_SE_RT_Corr = std(medianRT_correct)/sqrt(nsubs);
%these are smaller within subject error bars
%made by subtracting away each subjects average
%from their other scores to remove between subject difference
sub_mean_RT_Corr = mean(medianRT_correct,2); %average for each subject
%this subtracts each subjects average from their scores
%repmat repeats the matrix 4 times for each condition
mean_RT_Corr_deviation = medianRT_correct - repmat(sub_mean_RT_Corr,1,n_conditions);
%then take the standard error of those deviatoins from the mean
grand_withinSE_RT_Corr = std(mean_RT_Corr_deviation)/sqrt(nsubs);


%now do the same for proportion correct
grand_mean_prop_corr = mean(prop_correct);
grand_SE_prop_corr = std(prop_correct)/sqrt(nsubs);
sub_mean_prop_corr = mean(prop_correct,2);
prop_corr_deviation = prop_correct - repmat(sub_mean_prop_corr,1,n_conditions);
grand_withinSE_prop_corr = std(prop_corr_deviation)/sqrt(nsubs);
%% original 4 conditions plus proportion correct
%plot it
conds_plot = {'Pref_FaceIn'; 'Pref_FaceOut';'NonPref_FaceIn'; 'NonPref_FaceOut'}; 
figure;
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
subplot(1,2,1)
barweb(grand_mean_RT_Corr,grand_withinSE_RT_Corr);
ylim([450 525])
ylabel('Median RT (ms)')
title('Target Reaction Time (w/i subject SE)')
legend(conds)
subplot(1,2,2)
barweb(grand_mean_prop_corr,grand_withinSE_prop_corr);
ylim([.9 1])
ylabel('Proportion')
title('Proportion of Targets responded to')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STATS for all 4 conditions 
[p,tbl,stats] = anova1(medianRT_correct);


%repeated measures REAL STATS
srt1 = 1:1:length(medianRT_correct(:,1));
srt2 = 1:1:length(medianRT_correct(:,2));
srt3 = 1:1:length(medianRT_correct(:,3));
srt4 = 1:1:length(medianRT_correct(:,4));

condrt1 = ones(length(srt1),1)+0;
condrt2 = ones(length(srt2),1)+1;
condrt3 = ones(length(srt3),1)+2;
condrt4 = ones(length(srt4),1)+3;

Yrt = [medianRT_correct(:,1);medianRT_correct(:,2);medianRT_correct(:,3);medianRT_correct(:,4)];
Srt = [srt1';srt2';srt3';srt4'];
F1rt = [cond1;cond2;cond3;cond4];
data_anovart = [Yrt, F1rt, Srt];
RMAOV1(data_anovart)
%%
%side by side plots - PREFERENCE on X Axis
grand_meanRT_Pref = grand_mean_RT_Corr (1:2);
grand_meanRT_NPref = grand_mean_RT_Corr(3:4);
grand_wSE_RT_Pref = grand_withinSE_RT_Corr(1:2);
grand_wSE_RT_NPref = grand_withinSE_RT_Corr(3:4);

conds_plot = {'FacingIn'; 'FacingOut'}; 
figure;
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
subplot(1,2,1)
barweb(grand_meanRT_Pref,grand_wSE_RT_Pref);
ylim([450 525])
ylabel('Median RT (ms)')
xlabel ('Preferred Stance')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)
subplot(1,2,2)
conds_plot = {'FacingIn'; 'FacingOut'}; 
barweb(grand_meanRT_NPref,grand_wSE_RT_NPref);
ylim([450 525])
ylabel('Median RT (ms)')
xlabel('Non-preferred Stance')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATS for PREFERRED by stance (facing in vs out)
[h p ci stat] = ttest(medianRT_correct(:,1),medianRT_correct(:,2),.05,'both',1) 
%Stats for non-preferred by stance
[h p ci stat] = ttest(medianRT_correct(:,3),medianRT_correct(:,4),.05,'both',1) 
%%
%side by side plots - FACING on X Axis
gran_meanRT_Face_IN = grand_mean_RT_Corr(1:2:3);
gran_meanRT_Face_OUT = grand_mean_RT_Corr(2:2:4);
grand_wSE_RT_Face_IN = grand_withinSE_RT_Corr(1:2:3);
grand_wSE_RT_Face_OUT = grand_withinSE_RT_Corr(2:2:4);

conds_plot = {'Preferred'; 'Non-Preferred'}; 
figure;
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
subplot(1,2,1)
barweb(gran_meanRT_Face_IN,grand_wSE_RT_Face_IN);
ylim([450 525])
ylabel('Median RT (ms)')
xlabel ('Facing In')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)
subplot(1,2,2)
conds_plot = {'Preferred'; 'Non-Preferred'}; 
barweb(gran_meanRT_Face_OUT,grand_wSE_RT_Face_OUT);
ylim([450 525])
ylabel('Median RT (ms)')
xlabel ('Facing Out')
title('Target Reaction Time (w/i subject SE)')
legend(conds_plot)

%%%%%%%%%%%%%%%%
%Stats for facing inside (pref vs non-pref)
[h p ci stat] = ttest (medianRT_correct(:,1),medianRT_correct(:,3),.05,'both',1) 


%Stats for facing inside (pref vs non-pref)
[h p ci stat] = ttest(medianRT_correct(:,2),medianRT_correct(:,4),.05,'both',1) 


%% 
% GLOBAL STANCE 
%% global stance!!!!!!!!!!!!!!!!!!
%making a global comparison of preferred vs non-preferred 
global_meanRT_Pref = mean (grand_mean_RT_Corr(1:2));
global_meanRT_NPref = mean(grand_mean_RT_Corr(3:4));
global_mean_PNP = [global_meanRT_Pref; global_meanRT_NPref];

SE_P = std(mean(medianRT_correct(:,1:2),2))/sqrt(nsubs); 
SE_NP = std(mean(medianRT_correct(:,3:4),2))/sqrt(nsubs);
grand_withinSE_RT_Corr = std(mean_RT_Corr_deviation)/sqrt(nsubs);

global_SE_PNP = [7.8540; 10.0362];


%ENSURE THAT WITHIN STANDARD ERROR WAS CALCULATED PROPERLY
figure;
conds_plot = {'Preferred Stance'; 'Non-Preferred Stance'}; 
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(global_mean_PNP,global_SE_PNP);
ylim([450 535])
ylabel('Median RT (ms)')
xlabel ('Stance')
title('Target RT by Overall Preference')
legend(conds_plot)

%t-test
RT_pref_dif = mean(medianRT_correct(:,[1,2]),2)- mean(medianRT_correct(:,[3,4]),2)
[h p ci stat] = ttest(mean(medianRT_correct(:,[1,2]),2),mean(medianRT_correct(:,[3,4]),2),.05,'both',1) 
[h p ci stat] = ttest(RT_pref_dif,0,.05,'left',1) 



%%

%different pannels
%original conds acc
conds_plot_acc = {'Preferred Clockwise '; 'Preferred Counterclockwise';'Non-Preferred Clockwise'; 'Non-Preferred Counterclockwise'}; 
figure;
subplot (2,2,1)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(grand_mean_ACC_Corr,grand_withinSE_ACC_Corr);
ylim([85 100])
ylabel('Mean Accuracy')
xlabel ('Conditions')
%title('Target Accuracy')
legend(conds_plot_acc)

%global stance acc
subplot (2,2,2)
conds_plot_PNP = {'Preferred Stance'; 'Non-Preferred Stance'}; 
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(global_meanACC_PNP,global_SEacc_PNP);
ylim([85 100])
ylabel('Accuracy %')
xlabel ('Stance')
%title('Target Accuracy by Overall Preference')
legend(conds_plot_PNP)

%original conds RT
% original 4 conditions plus proportion correct
%plot it
conds_plotRT = {'Preferred Clockwise '; 'Preferred Counterclockwise';'Non-Preferred Clockwise'; 'Non-Preferred Counterclockwise'}; 
subplot(2,2,3)
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(grand_mean_RT_Corr,grand_withinSE_RT_Corr);
ylim([480 540])
ylabel('Median RT (ms)')
xlabel ('Conditions')
%title('Target Reaction Time (w/i subject SE)')
legend(conds_plotRT)

%RT by preference
subplot(2,2,4)
conds_plotRT_PNP = {'Preferred Stance'; 'Non-Preferred Stance'}; 
set(gcf,'color','w');
set(gcf, 'Position',  [100, 500, 1000, 400])
barweb(global_mean_PNP,global_SE_PNP);
ylim([480 540])
ylabel('Median RT (ms)')
xlabel ('Stance')
%title('Target RT by Overall Preference')
legend(conds_plotRT_PNP)



