ccc
exp = 'Skateboard';
subs = {'100'	'101'	'102'	'103'	'106'	'108'	'109'	'110'	'111'	'112'	'115'...
        '116'	'117'	'118'	'119'	'120'	'121'	'122'	'123'	'124'	'125'	'126'...
        '127'	'128'	'129'	'131'	'132'	'133'	'134'	'135'	'136'	'137'};
is_goofy = [0,	0,	0,	0,	0,	0,	1,	1,	0,	0,	1,	0,	0,	0,	0,	0,	0,	1,	0,...
            0,	0,	0,	0,	1,	1,	0,	0,	1,	1,	0,	1,	0];
        
nsubs = length(subs); 
conds =  {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds_lab = {'Preferred Clockwise'; 'Preferred Counterclockwise'; 'Non-preferred Clockwise'; 'Non-preferred Counterclockwise'};
pref_lab = {'Preferred'; 'Non-Preferred'};
facing_lab = {'Facing Inside'; 'Facing Outside'};
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\Winter2019/';
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_Corrected_Standard.set'],'filepath','M:\Data\Skateboard\Winter2019\segmentsFFT\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


    end
end
eeglab redraw

electrode = [1:16]; 

wavenumber = 12; %wavelet cycles
freqs = [.1:.5:30]; %wavelet frequencies - made freq bin larger
power_out = [];
i_count = 0;
n_electrode = EEG.nbchan;
for i_sub = 1:2%nsubs 
    fprintf(['Subject - ' num2str(i_sub) '. \n']);
    for i_cond = 1:nconds
        fprintf(['Condition - ' conds_lab{i_cond} '. \n']);
        i_count = i_count+1; %which data set in ALLEEG to use
        disp(ALLEEG(i_count).setname);
        power = [];
        n_trials = ALLEEG(i_count).trials;
        for i_trial = 1:n_trials
            for i_electrode = 1:length([electrode]) %pairs of contralateral electrodes
                tempdata = ALLEEG(i_count).data(electrode(i_electrode), :,i_trial) ;
                [temp_power temp_times temp_phase] = BOSC_tf(tempdata,freqs,EEG.srate,wavenumber);
                power(:,:,i_electrode,i_trial) = log(temp_power); %take log ?
            end
        end
        power_out(:,:,:,i_sub,i_cond) = nanmean(power,4); %save the power data, averaging across trials
    end
        if is_goofy(i_sub) == 1
        power_out_face_in(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[2,3]),5));
        power_out_face_out(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,4]),5));
    elseif is_goofy(i_sub) == 0
        power_out_face_in(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,4]),5));
        power_out_face_out(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[2,3]),5));
    end
    power_out_p(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,2]),5)); %averages preferred blocks 1:2
    power_out_np(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[3,4]),5)); %save the power data
end

%%
%%%now make topographies at each frequency bin
text_size = 8;
axis_label_size = 8;
title_size = 8;
line_width = 1;

F = 0:.5:30;
F = freqs;
delta_bins = find(F >= 0 & F <= 3);
theta_bins = find(F >= 4 & F <= 7);
alpha_bins = find(F >= 8 & F <= 12);
beta_bins = find(F >= 13 & F <= 30);

elec_locs = 'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced';%% Collapse over time to get power spectra, average over subjects

%%%P%%%
topo1 = squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,1),1),2),4));
topo1(17:18,1) = NaN;

%%%NP%%%
topo2 = squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,4),1),2),4));
topo2(17:18,1) = NaN;

plot_title = 'Alpha';


%%%Difference%%%
topo3 = topo1-topo2;

min_lim = min([topo1;topo2])
max_lim = max([topo1;topo2])

diff_min_lim = min([topo1-topo2])
diff_max_lim = max([topo1-topo2])


figure;
topoplot(topo1,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
set(gca,'FontSize',text_size,'FontWeight', 'bold');
title(['No Headset : ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0],'FontSize',axis_label_size);
set(t,'LineWidth',line_width);

figure;
topoplot(topo2,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
set(gca,'FontSize',text_size,'FontWeight', 'bold');
title(['Headset : ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0],'FontSize',axis_label_size);
set(t,'LineWidth',line_width);

figure;
topoplot(topo3,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[-0.5 0.5])
set(gca,'FontSize',text_size,'FontWeight', 'bold');
title(['No Headset - Headset: ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0],'FontSize',axis_label_size);
set(t,'LineWidth',line_width);

% 
% %%
% time_window = find(EEG.times>350,1)-1:find(EEG.times>550,1)-2;
% power_spectra = squeeze(mean(power_out,2)); %collapse over time can pick time range here in 2nd dim.
% power_spectra = squeeze(mean(power_out(:,time_window,:,:,:),2)); %collapse over time can pick time range here in 2nd dim.
% 
% %frequency x electrode x subject x condition
% power_spectra_p = squeeze(mean(power_out_p,2));
% power_spectra_np = squeeze(mean(power_out_np,2));
% power_spectra_face_in = squeeze(mean(power_out_face_in,2)); 
% power_spectra_face_out = squeeze(mean(power_out_face_out,2)); 
% 
% 
%                        
%                         
% power_spectra_mean = squeeze(mean(power_spectra,2));
% power_spectra_se = squeeze(std(power_spectra,[],2)/sqrt(nsubs));
% power_spectra_mean_p = squeeze(mean(power_spectra_p,2));
% power_spectra_se_p = squeeze(std(power_spectra_p,[],2)/sqrt(nsubs));
% power_spectra_mean_np = squeeze(mean(power_spectra_np,2));
% power_spectra_se_np = squeeze(std(power_spectra_np,[],2)/sqrt(nsubs));
% power_spectra_mean_facing_in = squeeze(mean(power_spectra_face_in,2));
% power_spectra_se_facing_in = squeeze(std(power_spectra_face_in,[],2)/sqrt(nsubs));
% power_spectra_mean_facing_out = squeeze(mean(power_spectra_face_out,2));
% power_spectra_se_facing_out = squeeze(std(power_spectra_face_out,[],2)/sqrt(nsubs));
% 
% %% Plot spectra accross conditions
% % close all
% figure;
% boundedline(freqs,power_spectra_mean(:,1),power_spectra_se(:,1),'k',...
%     freqs,power_spectra_mean(:,2),power_spectra_se(:,2),'m',...
%     freqs,power_spectra_mean(:,3),power_spectra_se(:,3),'g',...
%     freqs,power_spectra_mean(:,4),power_spectra_se(:,4),'b');
% axis tight
% xlim([0 30])
% % ylim ([0 3000000])
% xlabel('Frequency (Hz)');
% ylabel('Power (uV^2)');
% title('Power PZ');
% legend(conds_lab,'Location','NorthEast');
% 
% %plot spectra by preference 
% figure;
% boundedline (freqs, power_spectra_mean_p, power_spectra_se_p,'r',...
%              freqs, power_spectra_mean_np, power_spectra_se_np,'b');
% axis tight
% xlim([0 30])
% % ylim ([0 3000000])
% xlabel('Frequency (Hz)');
% ylabel('Power (uV^2)');
% title('Power by Preference');
% legend(pref_lab,'Location','NorthEast'); 
% 
% %plot spectra by facing orientation 
% figure;
% boundedline (freqs, power_spectra_mean_facing_in, power_spectra_se_facing_in,'k',...
%              freqs, power_spectra_mean_facing_out, power_spectra_se_facing_out,'g');
% axis tight
% xlim([0 30])
% % ylim ([0 3000000])
% xlabel('Frequency (Hz)');
% ylabel('Power (uV^2)');
% title('Power by Orientation');
% legend(facing_lab,'Location','NorthEast');
