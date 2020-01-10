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

electrode = 15; 

wavenumber = 12; %wavelet cycles
freqs = [1:1:30]; %wavelet frequencies
power_out = [];
i_count = 0;
n_electrode = EEG.nbchan;
for i_sub = 1:nsubs 
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
                power(:,:,i_electrode,i_trial) = temp_power; %take log ?
            end
        end
        power_out(:,:,:,i_sub,i_cond) = nanmean(power,4); %save the power data, averaging across trials
    end
        if is_goofy(i_sub) == 1
        power_out_face_in(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[2,3]),5))
        power_out_face_out(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,4]),5))
    elseif is_goofy(i_sub) == 0
        power_out_face_in(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,4]),5))
        power_out_face_out(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[2,3]),5))
    end
    power_out_p(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,2]),5)); %averages preferred blocks 1:2
    power_out_np(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[3,4]),5)); %save the power data
end

%% Collapse over time to get power spectra, average over subjects


power_spectra = squeeze(mean(power_out,2)); %collapse over time can pick time range here in 2nd dim.
%frequency x electrode x subject x condition
power_spectra_p = squeeze(mean(power_out_p,2));
power_spectra_np = squeeze(mean(power_out_np,2));
power_spectra_face_in = squeeze(mean(power_out_face_in,2)); 
power_spectra_face_out = squeeze(mean(power_out_face_out,2)); 


                       
                        
power_spectra_mean = squeeze(mean(power_spectra,2));
power_spectra_se = squeeze(std(power_spectra,[],2)/sqrt(nsubs));
power_spectra_mean_p = squeeze(mean(power_spectra_p,2));
power_spectra_se_p = squeeze(std(power_spectra_p,[],2)/sqrt(nsubs));
power_spectra_mean_np = squeeze(mean(power_spectra_np,2));
power_spectra_se_np = squeeze(std(power_spectra_np,[],2)/sqrt(nsubs));
power_spectra_mean_facing_in = squeeze(mean(power_spectra_face_in,2));
power_spectra_se_facing_in = squeeze(std(power_spectra_face_in,[],2)/sqrt(nsubs));
power_spectra_mean_facing_out = squeeze(mean(power_spectra_face_out,2));
power_spectra_se_facing_out = squeeze(std(power_spectra_face_out,[],2)/sqrt(nsubs));

%% Plot spectra accross conditions
close all
figure;
boundedline(freqs,power_spectra_mean(:,1),power_spectra_se(:,1),'k',...
    freqs,power_spectra_mean(:,2),power_spectra_se(:,2),'m',...
    freqs,power_spectra_mean(:,3),power_spectra_se(:,3),'g',...
    freqs,power_spectra_mean(:,4),power_spectra_se(:,4),'b');
axis tight
xlim([0 30])
ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title('Power PZ');
legend(conds_lab,'Location','NorthEast');

%plot spectra by preference 
figure;
boundedline (freqs, power_spectra_mean_p, power_spectra_se_p,'r',...
             freqs, power_spectra_mean_np, power_spectra_se_np,'b');
axis tight
xlim([0 30])
ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title('Power by Preference');
legend(pref_lab,'Location','NorthEast'); 

%plot spectra by facing orientation 
figure;
boundedline (freqs, power_spectra_mean_facing_in, power_spectra_se_facing_in,'k',...
             freqs, power_spectra_mean_facing_out, power_spectra_se_facing_out,'g');
axis tight
xlim([0 30])
ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title('Power by Orientation');
legend(facing_lab,'Location','NorthEast');
