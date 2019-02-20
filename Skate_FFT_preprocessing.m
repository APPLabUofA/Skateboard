% Strategy will be, load in preprocessed epochs -1000 to 2000 ms, 
% for each epoch go through each channel, run an FFT on the -1000 to 0 data 
%(prestimulus), get the spectra, 
% Then subtract the spectra for left and right side of the head
% Then average those over epochs
% 
% If you start making a script to load in the preprocessed .set files,...
...and indicate where these steps should go, I can help with writing a draft of a script to run this analysis
    
%%
clear all
close all
ccc

exp = 'Skateboard';
%subs = {'100' '101' '102' '103' '104' '106' '107' '108' '109' '110' '111'...
...'112' '113' '114' '115' '016'};
subs = {'117'}; %to test on just one sub

nsubs = length(subs);
conds = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};%preferred, clockwise - non-preffered, CCW
%conds = {'G1'; 'G2';'R1';'R2'};
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\winter2019\';

if exist([Pathname 'segments\'])
    mkdir([Pathname 'segmentsFFT\']);
end
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond} '.vhdr'];
        setname = Filename(1:end-5)
        
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
               
        EEG = pop_epoch( EEG, {  '3'  '5'  }, [-1  2], 'newname',  sprintf('%s epochs' , setname), 'epochinfo', 'yes'); %Changed from [-.2 1] to [-1 2]. DR
        EEG = pop_rmbase( EEG, [-200    0]);
        
        %         eeglab redraw
             
        
        %    Artifact rejection, trials with range >500 uV
        EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)],-1000,1000,EEG.xmin,EEG.xmax,0,1);
        
        %   EMCP occular correction
        temp_ocular = EEG.data(end-1:end,:,:); %to save the EYE data for after
        selection_cards = {'3','5' }; %different bin names, each condition should be separate
        EEG = gratton_emcp(EEG,selection_cards,{'VEOG'},{'HEOG'}); %this assumes the eye channels are called this
        EEG.emcp.table %this prints out the regression coefficients
        EEG.data(end-1:end,:,:) = temp_ocular; %replace the eye data
        %    Artifact rejection, trials with range >250 uV
        EEG = pop_rmbase( EEG, [-200 0]); %baseline again since this changed it
        EEG = pop_eegthresh(EEG,1,[1:size(EEG.data,1)-2],-500,500,EEG.xmin,EEG.xmax,0,1);
        
        tempEEG =   EEG;
        
        
        %now select the corrected trials
        EEG = pop_selectevent( tempEEG, 'type',5,'renametype','Target','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Target']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Target.set'],'filepath',[Pathname 'segmentsFFT\']);
        
        
        EEG = pop_selectevent( tempEEG, 'type',3 ,'renametype','Standard','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Standard']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Standard.set'],'filepath',[Pathname 'segmentsFFT\']);
        
        
    end 
end 
