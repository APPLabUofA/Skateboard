

clear all
close all
% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
%% Input Information

 exp = 'Skateboard';

% *************************************************************************
subs = {'100'	'101'	'102'	'103'	'104'	'106'	'108'	'109'	'110'...
    '111'	'112'	'113'	'115'	'116'	'117'	'119'	'120'...
    '122'	'123'	'125'	'126'	'127'	'129'	...
    '131'	'132'	'133'	'134'	'136'	'137'};
%subs = {'100' '101' '102'}; %to test on just one sub 
nsubs = length(subs); % number of subjects being run
% *************************************************************************
conds =  {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
%conds = {'P_CW';'NP_CW'};
nconds = length(conds);
% *************************************************************************

Pathname = 'M:\Data\Skateboard\Winter2019/';

% Location of electrode information
electrode_loc = 'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced';

% *************************************************************************
% A few electrodes
 electrode = [13 15];
 elec_names = {'Fz';'Pz'};

% Multiple electrodes
% electrode = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
% elec_names = {'Oz';'P7';'T7';'P3';'C3';'F3';'Pz';'Cz';'Fz';'P4';'C4';'F4';'P8';'T8';'FP2'};

% *************************************************************************
% Pick the type of trials
%   1 = standards
%   2 = targets
%   3 = standards and targets
perms = 3;
trialevent = {'Standards';'Targets'};
%trialevent = {'Standards'};
%trialevent = {'Targets'};

% *************************************************************************
% If using eeglab to plot ersp
elab_plot = 'Off'; % No

% *************************************************************************
% Set baseline
baseln = [-1000 -500];
%baseln = [NaN];

% *************************************************************************
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%time frequency parameters
tf_epochslim = [-1  1.9];
tf_cycles = 3;
dum1 = [];
dum2 = [];
timesout = 200;
analfreq = [4 20];
%
%-------------------------------------------------------------------------
%% TF Analysis
% -------------------------------------------------------------------------
% clear ersp freqs itc powbase times

i_count = 0;

for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        i_count = i_count + 1; % counter to select data from ALLEEG
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
            
        if perms == 1 % Load standards data
            EEG = pop_loadset('filename',[Filename '_fft_Standards.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_fft_JK\');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            
            for i_chan = 1:length(electrode)
                EEG = eeg_checkset(EEG);
                % Caption if creating plots with pop_newtimef
%                 i_cap = strcat(exp, '_', subs(i_sub), '_', elec_names(i_chan),...
%                     '_', conds(i_cond));
                % If plotting with the pop_newtimef function
                if strcmp(elab_plot, 'Yes')
                    figure
                end
                % FFT with output ersp & itc in the order of each subj, cond, 
                % trial type (perm), and channel.
                [ersp(i_sub,i_cond,1,i_chan,:,:),itc(i_sub,i_cond,1,i_chan,:,:),powbase,times,freqs] =...
                    pop_newtimef(EEG, 1, i_chan, tf_epochslim*1000, tf_cycles, 'topovec', i_chan,...
                    'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'baseline', baseln,...
                    'freqs', analfreq, 'plotersp', elab_plot, 'plotitc', elab_plot,...
                    'padratio', 8, 'timesout', timesout);
            end
            
        elseif perms == 2 % Load targets data
            EEG = pop_loadset('filename',[Filename '_fft_Targets.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_fft_JK\');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            
            for i_chan = 1:length(electrode)
                EEG = eeg_checkset(EEG);
                % Caption if creating plots with pop_newtimef
%                 i_cap = strcat(exp, '_', subs(i_sub), '_', elec_names(i_chan),...
%                     '_', conds(i_cond));
                % If plotting with the pop_newtimef function
                if strcmp(elab_plot, 'Yes')
                    figure
                end
                % FFT with output ersp & itc in the order of each subj, cond, 
                % trial type (perm), and channel.
                [ersp(i_sub,i_cond,1,i_chan,:,:),itc(i_sub,i_cond,1,i_chan,:,:),powbase,times,freqs] =...
                    pop_newtimef(EEG, 1, i_chan, tf_epochslim*1000, tf_cycles, 'topovec', i_chan,...
                    'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'baseline', baseln,...
                    'freqs', analfreq, 'plotphase', 'Off', 'plotersp', elab_plot,...
                    'plotitc', elab_plot, 'padratio', 8, 'timesout', timesout);
            end
            
        elseif perms == 3
            
            % Load standards data
             EEG = pop_loadset('filename',[Filename '_fft_Standards.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_fft_JK\');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            for i_chan = 1:length(electrode)
                EEG = eeg_checkset(EEG);
                % Caption if creating plots with pop_newtimef
%                 i_cap = strcat(exp, '_', subs(i_sub), '_', elec_names(i_chan),... commented this out for now DR
%                     '_', conds(i_cond));
                % If plotting with the pop_newtimef function
                if strcmp(elab_plot, 'Yes')
                    figure
                end
                % FFT with output ersp & itc in the order of each subj, cond, 
                % trial type (perm), and channel.
                [ersp(i_sub,i_cond,1,i_chan,:,:),itc(i_sub,i_cond,1,i_chan,:,:),powbase,times,freqs] =...
                    pop_newtimef(EEG, 1, i_chan, tf_epochslim*1000, tf_cycles, 'topovec', i_chan,...
                    'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'baseline', baseln,...
                    'freqs', analfreq, 'plotphase', 'Off', 'plotersp', elab_plot,...
                    'plotitc', elab_plot, 'padratio', 8, 'timesout', timesout);
            end
            
            % Load targets data
            EEG = pop_loadset('filename',[Filename '_fft_Targets.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_fft_JK\');
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            for i_chan = 1:length(electrode)
                EEG = eeg_checkset(EEG);
                % Caption if creating plots with pop_newtimef
%                 i_cap = strcat(exp, '_', subs(i_sub), '_', elec_names(i_chan),...
%                     '_', conds(i_cond));
                % If plotting with the pop_newtimef function
                if strcmp(elab_plot, 'Yes')
                    figure
                end
                % FFT with output ersp & itc in the order of each subj, cond, 
                % trial type (perm), and channel.
                [ersp(i_sub,i_cond,2,i_chan,:,:),itc(i_sub,i_cond,2,i_chan,:,:),powbase,times,freqs] =...
                    pop_newtimef(EEG, 1, i_chan, tf_epochslim*1000, tf_cycles, 'topovec', i_chan,...
                    'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'baseline', baseln,...
                    'freqs', analfreq, 'plotphase', 'Off', 'plotersp', elab_plot,...
                    'plotitc', elab_plot, 'padratio', 8, 'timesout', timesout);
            end
            
        end     
          
    end
 
end

eeglab redraw

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -------------------------------------------------------------------------
%% Plotting functions

%  Skate_ERSP_Plot(ALLEEG, 1, EEG, elec_names, electrode, electrode_loc,...
%      ersp, exp, Filename, freqs, itc, Pathname, perms, powbase, subs, times,...
%      trialevent)
%  
% Skate_ERSP_Plot(freqs,ersp)

% -------------------------------------------------------------------------
% /////////////////////////////////////////////////////////////////////////
% -----------------------------------------------------------------------

%%
%QUICK PLOTS

CLim = [-1.5 1.5];

for i_event = 1:length(trialevent)

    for ch_ersp = 1:length(electrode)

        % ERSP values by electrode
        % (participants x conditions x events x electrodes x frequencies x timepoints)
        pref_ersp_chan  =  squeeze(mean(mean(ersp(:,1:2,i_event,ch_ersp,:,:),1),2));
        npref_ersp_chan = squeeze(mean(mean(ersp(:,3:4,i_event,ch_ersp,:,:),1),2));

        figure; 
        colormap('jet')

        % Subplot 1: out
        subplot(3,1,1); 
        imagesc(times,freqs,npref_ersp_chan,CLim);
        title({['ERSP: ' char(trialevent{i_event}) ', NonPref']; char(elec_names(ch_ersp))}); 
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
        ylabel('Freq (Hz)');
        colorbar
        % Subplot 1: in
        subplot(3,1,2); 
        imagesc(times,freqs,pref_ersp_chan,CLim);
        title({['ERSP: ' char(trialevent{i_event}) ', Pref']; char(elec_names(ch_ersp))}); 
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
        ylabel('Freq (Hz)');
        colorbar
        % Subplot 1: out-in
        subplot(3,1,3); 
        imagesc(times,freqs,npref_ersp_chan-pref_ersp_chan,CLim); 
        title({['ERSP: ' char(trialevent{i_event}) ', np-p']; char(elec_names(ch_ersp))});
        set(gca,'Ydir','Normal')
        line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5) 
        ylabel('Freq (Hz)'); xlabel('Time (ms)');
        colorbar

        %Overall subplot title
    %     supertitle(['ERSP: ' char(elec_names(ch_ersp))],'FontSize',12)

        clear in_ersp_chan out_ersp_chan 
    end

    clear ch_ersp 

end

clear i_event

%% ***************************************************************
%PLOTS USING DIFFERENCE ERSPS AND AVERAGED BY PREFERENCE 
%erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
%time_window = find(EEG.times>-200,1)-1:find(EEG.times>1000,1)-2;

%averaging ERSP windows and subtracting targets-standards
% (participants x conditions x events x electrodes x frequencies x timepoints)

% ersp_dif = squeeze(ersp(:,:,2,:,:,:)-ersp(:,:,1,:,:,:));

pref_ersp_stand  =  squeeze(mean(mean(ersp(:,[1,2],1,2,:,:),1),2));
pref_ersp_targ  =  squeeze(mean(mean(ersp(:,[1,2],2,2,:,:),1),2));
pref_ersp_dif = squeeze(pref_ersp_targ-pref_ersp_stand);

npref_ersp_stand = squeeze(mean(mean(ersp(:,[3,4],1,2,:,:),1),2));
npref_ersp_targ = squeeze(mean(mean(ersp(:,[3,4],2,2,:,:),1),2));
npref_ersp_dif = squeeze(npref_ersp_targ-npref_ersp_stand);


CLim = [-1.5 1.5];
figure;
colormap('jet')

% Subplot 1: pref
subplot(3,1,1);
imagesc(times,freqs,pref_ersp_stand,CLim);
title('ERSP:Standards, Preferred @Pz');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[12 12],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: in
subplot(3,1,2);
imagesc(times,freqs,pref_ersp_targ,CLim);
title('ERSP:Targets, Preferred');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[12 12],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: out-in
subplot(3,1,3);
imagesc(times,freqs,pref_ersp_dif,CLim);
title('ERSP:Targets-Standards, Preferred');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[12 12],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar

%NON-PREF
CLim = [-1.5 1.5];
figure;
colormap('jet')

% Subplot 1: pref
subplot(3,1,1);
imagesc(times,freqs,npref_ersp_stand,CLim);
title('ERSP:Standards, Non Preferred @ Pz');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[12 12],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: in
subplot(3,1,2);
imagesc(times,freqs,npref_ersp_targ,CLim);
title('ERSP:Targets, Non Preferred');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[12 12],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: out-in
subplot(3,1,3);
imagesc(times,freqs,npref_ersp_dif,CLim);
title('ERSP:Targets-Standards, Non Preferred');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[12 12],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar

%preferred unpreferred difference wave 
pnp_grand_diff_pz = squeeze(pref_ersp_dif-npref_ersp_dif);
figure
subplot(3,1,1);
imagesc(times,freqs,pnp_grand_diff_pz,CLim);
title('ERSP:Preferred-Unpreferred, Target-Standards');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[12 12],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar

%% *******************************************************************************
%  THETA PLOTS 

pref_ersp_stand_fz  =  squeeze(mean(mean(ersp(:,[1,2],1,1,:,:),1),2));
pref_ersp_targ_fz  =  squeeze(mean(mean(ersp(:,[1,2],2,1,:,:),1),2));
pref_ersp_dif_fz = squeeze(pref_ersp_targ_fz-pref_ersp_stand_fz);

npref_ersp_stand_fz = squeeze(mean(mean(ersp(:,[3,4],1,1,:,:),1),2));
npref_ersp_targ_fz = squeeze(mean(mean(ersp(:,[3,4],2,1,:,:),1),2));
npref_ersp_dif_fz = squeeze(npref_ersp_targ_fz-npref_ersp_stand_fz);


CLim = [-1.5 1.5];
figure;
colormap('jet')

% Subplot 1: pref
subplot(3,1,1);
imagesc(times,freqs,pref_ersp_stand_fz,CLim);
title('ERSP:Standards, Preferred @Fz');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[4 4],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: in
subplot(3,1,2);
imagesc(times,freqs,pref_ersp_targ_fz,CLim);
title('ERSP:Targets, Preferred');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[4 4],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: out-in
subplot(3,1,3);
imagesc(times,freqs,pref_ersp_dif_fz,CLim);
title('ERSP:Targets-Standards, Preferred');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[4 4],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar

%NON-PREF
CLim = [-1.5 1.5];
figure;
colormap('jet')

% Subplot 1: pref
subplot(3,1,1);
imagesc(times,freqs,npref_ersp_stand_fz,CLim);
title('ERSP:Standards, Non Preferred @ Fz');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[4 4],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: in
subplot(3,1,2);
imagesc(times,freqs,npref_ersp_targ_fz,CLim);
title('ERSP:Targets, Non Preferred');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[4 4],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)');
colorbar
% Subplot 1: out-in
subplot(3,1,3);
imagesc(times,freqs,npref_ersp_dif_fz,CLim);
title('ERSP:Targets-Standards, Non Preferred');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[4 4],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar

pnp_grand_diff_fz = squeeze(pref_ersp_dif_fz-npref_ersp_dif_fz);
figure
subplot(3,1,1);
imagesc(times,freqs,pnp_grand_diff_fz,CLim);
title('ERSP:Preferred-Unpreferred, Target-Standards');
set(gca,'Ydir','Normal')
line([0 0],[min(freqs) max(freqs)],'Color','m','LineStyle','--','LineWidth',1.5)
line([min(times) max(times)],[4 4],'Color','k','LineStyle','--','LineWidth',.1)
line([min(times) max(times)],[8 8],'Color','k','LineStyle','--','LineWidth',.1)
ylabel('Freq (Hz)'); xlabel('Time (ms)');
colorbar

