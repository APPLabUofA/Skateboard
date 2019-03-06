%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    REGULAR IN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all
close all
ccc

exp = 'Skateboard';

% loading regulars first
subs = {'100' '101' '102' '103' '104' '106' '108' '111'...
  '112' '116' '117' '118' '119' '120'};
%subs = {'100'}; %to test on just one sub

 nsubs = length(subs);
%conds = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds =  {'P_CW';'NP_CCW'}; %facing inside for regulars
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\winter2019\';

if exist([Pathname 'segments\'])
    mkdir([Pathname 'segmentsIO\']);
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
        
        if i_cond ==1
        tempEEG1 =   EEG;
        elseif i_cond==2
        tempEEG2 = EEG;
        
        tempEEG = pop_mergeset(tempEEG1,tempEEG2,1);
        
        %now select the corrected trials
        EEG = pop_selectevent( tempEEG, 'type',5,'renametype','Target','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Target']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_facing_in_Corrected_Target.set'],'filepath',[Pathname 'segmentsIO\']);
        
        
        EEG = pop_selectevent( tempEEG, 'type',3 ,'renametype','Standard','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Standard']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_facing_in_Corrected_Standard.set'],'filepath',[Pathname 'segmentsIO\']);
        
        end
        
    end 
end 
commandhistory
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    REGULAR OUT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all
close all
ccc

exp = 'Skateboard';

% loading regulars first
subs = {'100' '101' '102' '103' '104' '106' '108' '111'...
  '112' '116' '117' '118' '119' '120'};
%subs = {'100'}; %to test on just one sub

 nsubs = length(subs);
%conds = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds =  {'P_CCW';'NP_CW'}; %facing OUTSIDE for regulars
conds2 = {'facing_out1'; 'facing_out2'};

nconds = length(conds);
Pathname = 'M:\Data\Skateboard\winter2019\';

if exist([Pathname 'segments\'])
    mkdir([Pathname 'segmentsIO\']);
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
        
        if i_cond ==1
        tempEEG1 =   EEG;
        elseif i_cond==2
        tempEEG2 = EEG;
        
        tempEEG = pop_mergeset(tempEEG1,tempEEG2,1);
        
        
        %now select the corrected trials
        EEG = pop_selectevent( tempEEG, 'type',5,'renametype','Target','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Target']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_facing_out_Corrected_Target.set'],'filepath',[Pathname 'segmentsIO\']);
        
        
        EEG = pop_selectevent( tempEEG, 'type',3 ,'renametype','Standard','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Standard']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_facing_out_Corrected_Standard.set'],'filepath',[Pathname 'segmentsIO\']);
        end
        
    end 
end 
commandhistory
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    GOOFY IN
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all
close all
ccc

exp = 'Skateboard';

% GOOFY SUBS
subs = {'107' '109' '110' '113' '114' '115' '122'};
%subs = {'100'}; %to test on just one sub

 nsubs = length(subs);
%conds = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds =  {'P_CCW';'NP_CW'}; %facing in - GOOFY
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\winter2019\';

if exist([Pathname 'segments\'])
    mkdir([Pathname 'segmentsIO\']);
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
        
        if i_cond ==1
        tempEEG1 =   EEG;
        elseif i_cond==2
        tempEEG2 = EEG;
        
        tempEEG = pop_mergeset(tempEEG1,tempEEG2,1);  
        
        %now select the corrected trials
        EEG = pop_selectevent( tempEEG, 'type',5,'renametype','Target','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Target']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_facing_in_Corrected_Target.set'],'filepath',[Pathname 'segmentsIO\']);
        
        
        EEG = pop_selectevent( tempEEG, 'type',3 ,'renametype','Standard','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Standard']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_facing_in_Corrected_Standard.set'],'filepath',[Pathname 'segmentsIO\']);
        
        end 
    end 
end 
commandhistory

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    GOOFY OUT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all
close all
ccc

exp = 'Skateboard';

% GOOFY SUBS
subs = {'107' '109' '110' '113' '114' '115' '122'};
%subs = {'100'}; %to test on just one sub

 nsubs = length(subs);
%conds = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds =  {'P_CW';'NP_CCW'}; %facing OUT - GOOFY
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\winter2019\';

if exist([Pathname 'segments\'])
    mkdir([Pathname 'segmentsIO\']);
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
        
        if i_cond ==1
        tempEEG1 =   EEG;
        elseif i_cond==2
        tempEEG2 = EEG;
        
        tempEEG = pop_mergeset(tempEEG1,tempEEG2,1);
        
        
        %now select the corrected trials
        EEG = pop_selectevent( tempEEG, 'type',5,'renametype','Target','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Target']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_facing_out_Corrected_Target.set'],'filepath',[Pathname 'segmentsIO\']);
        
        
        EEG = pop_selectevent( tempEEG, 'type',3 ,'renametype','Standard','deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_editset(EEG, 'setname',[subs{i_sub} '_' exp '_' conds{i_cond} '_Corrected_Standard']);
        EEG = pop_saveset( EEG, 'filename',[subs{i_sub} '_' exp '_facing_out_Corrected_Standard.set'],'filepath',[Pathname 'segmentsIO\']);
        end
        
    end 
end 
commandhistory