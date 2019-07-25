clear all
close all
ccc

exp = 'Skateboard';
subs = {'100' '102' '103' '104' '107' '108' '109' ...
    '111' '113' '116' '117' '118' '119'};% ...
    ...'125' '126'};
% subs = {'100' '101' '102' '103' '104' '106' '107' '108' '109' ...
%     '110' '111' '112' '113' '114' '115' '116' '117' '118' '119' ...
%     '120' '122' '123' '124' '125' '126' '127' '128' '129'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%I HAVE NOT INTEGRATED A VERSION OF THE CODE THAT PREPROCESSES BOTH SET
%OF TRIGGERS AT ONCE. SINCE WE RAN INTO THE TRIGGER ISSUE LATER IN
%THE STUDY I AM NOT WORRYING ABOUT PREPROCESSING ALL THE DATA THAT'S BEEN
%SUCESSFULLY PREPROCESSED. HENCE, SUBJECTS 124 TO 129 WERE PREPROCESSED
%SEPARATELY AND CHANGING THE TRIGGER VALUES MANUALLY. DANIEL.
% These settings are meant to create a new set of data stored by preferred
% vs non-preferred stances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsubs = length(subs);
conds = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
nconds = length(conds);
new_conds = {'preferred'; 'non-preferred'};
Pathname = 'M:\Data\Skateboard\winter2019\';

if  exist([Pathname 'segments\'])
    mkdir([Pathname 'segments_P_NP\']);
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond} '.vhdr'];
        setname = Filename(1:end-5);
        
        EEG = pop_loadbv(Pathname, Filename);
        
        % get electrode locations
        EEG=pop_chanedit(EEG, 'load',{'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced' 'filetype' 'autodetect'});
        
        % arithmetically rereference to linked mastoid
        for x=1:EEG.nbchan-2
            EEG.data(x,:) = (EEG.data(x,:)-((EEG.data(EEG.nbchan-2,:))*.5));
        end
        
        %Filter the data with low pass of 30
        EEG = pop_eegfilt( EEG, .1, 0, [], 0);  %high pass filter
        EEG = pop_eegfilt( EEG, 0, 30, [], 0);  %low pass filter
        
        
        all_events = length(EEG.event)
        for i_event = 2:all_events
            EEG.event(i_event).type = num2str(str2num((EEG.event(i_event).type(2:end))));
        end
        
        %epoch
        
        EEG = pop_epoch( EEG, {  '3'  '5'  }, [-.2  1], 'newname',  sprintf('%s epochs' , setname), 'epochinfo', 'yes'); %Changed from [-.2 1] to [-1 2]. DR
        EEG = pop_rmbase( EEG, [-200    0]);
        
        %         eeglab redraw
        
        %    Artifact rejection, trials with range >500 uV
        EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],-500,500,EEG.xmin,EEG.xmax,0,1);
        
        %   EMCP occular correction
        temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
        selection_cards = {'3','5' }; %different bin names, each condition should be separate
        EEG = gratton_emcp(EEG,selection_cards,{'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
        EEG.emcp.table %this prints out the regression coefficients
        EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data
        %    Artifact rejection, trials with range >250 uV
        EEG = pop_rmbase( EEG, [-200 0]); %baseline again since this changed it
        EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)-2],-200,200,EEG.xmin,EEG.xmax,0,1);
        
        tempEEG = EEG;
        
                 
        if strcmp(conds{i_cond},'P_CCW') || strcmp(conds{i_cond},'P_CW')
            ThisIsTheName = [subs{i_sub} '_' exp '_' new_conds{1} '_Corrected_Standard.set']
        elseif strcmp(conds{i_cond},'NP_CCW') || strcmp(conds{i_cond},'NP_CW')
            ThisIsTheName = [subs{i_sub} '_' exp '_' new_conds{2} '_Corrected_Standard.set']
        end
        
        
        EEG = pop_selectevent( tempEEG, 'type',3 ,'renametype','Standard','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',ThisIsTheName);
        EEG = pop_saveset( EEG, 'filename',ThisIsTheName,'filepath',[Pathname 'segments_P_NP\']);
        
        
        
        if strcmp(conds{i_cond},'P_CCW') || strcmp(conds{i_cond},'P_CW')
            ThisIsTheName1 = [subs{i_sub} '_' exp '_' new_conds{1} '_Corrected_Target.set']
        elseif strcmp(conds{i_cond},'NP_CCW') || strcmp(conds{i_cond},'NP_CW')
            ThisIsTheName1 = [subs{i_sub} '_' exp '_' new_conds{2} '_Corrected_Target.set']
        end
        
                    
        EEG = pop_selectevent( tempEEG, 'type',5,'renametype','Target','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',ThisIsTheName1);
        EEG = pop_saveset( EEG, 'filename',ThisIsTheName1,'filepath',[Pathname 'segments_P_NP\']);
        
        
    end
end


