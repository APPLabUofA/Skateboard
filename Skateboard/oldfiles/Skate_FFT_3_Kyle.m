%CONDITIONS...
...preferred, clockwise - non-preffered, CCW
%%    
ccc
%
exp = 'Skateboard';
subs = {'100' '101' '102'};
%subs = {'101'}; %to test on just one sub 

nsubs = length(subs); 
conds =  {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds_lab = {'Preferred Clockwise'; 'Preferred Counterclockwise'; 'Non-preferred Clockwise'; 'Non-preferred Counterclockwise'};
%{'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\Winter2019';
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_Corrected_Standard.set'],'filepath','M:\Data\Skateboard\Winter2019\segmentsFFT\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


    end
end
eeglab redraw

%% loop through the sets and do a wavelet analysis
% Kyle TODO next - NEXT - loop this through electrodes
% Kyle TODO next - compare left vs right hemisphere


wavenumber = 3; %wavelet cycles
electrode = 1;% pz
freqs = [.1:.5:30]; %wavelet frequencies
power_out = [];
i_count = 0;

for i_sub = 1:nsubs 
    fprintf(['Subject - ' num2str(i_sub) '. \n']);
    for i_cond = 1:nconds
        fprintf(['Condition - ' conds_lab{i_cond} '. \n']);

        i_count = i_count+1; %which data set in ALLEEG to use
        n_trials = ALLEEG(i_count).trials;
        for i_trial = 1:n_trials
            tempdata = ALLEEG(i_count).data(electrode, :,i_trial) ;
            power = [];
            phase = [];
            [temp_power temp_times temp_phase] = BOSC_tf(tempdata,freqs,EEG.srate,wavenumber);
            power(:,:,i_trial) = log(temp_power); %take log 
            
        end
        power_out(:,:,i_sub,i_cond) = nanmean(power,3); %save the power data
 
        
    end
 
end


%% plot spectrograms for each subject (rows) and condition (columns)
 figure;
 i_count = 0;
for i_sub = 1:nsubs
    
    for i_cond = 1:nconds
        i_count = i_count+1;
        
        subplot(nsubs,nconds,i_count)
        imagesc(EEG.times,freqs,squeeze(power_out(:,:,i_sub,i_cond)));
        set(gca,'Ydir','normal');
        title(['Sub ' num2str(i_sub) '-' conds_lab{i_cond}])
    end
end
xlabel('Time (s)');
ylabel('Frequency (Hz)');



%% Collapse over time to get power spectra, average over subjects

power_spectra = squeeze(mean(power_out(:,1:find(EEG.times>0,1)-1,:,:),2)); %collapse over time before 0 ms

mean_power_out = squeeze(mean(power_spectra,2));
stderr_power_out = squeeze(std(power_out,[],2))./sqrt(nsubs);


%% Plot spectra accross conditions
figure; 

boundedline(freqs,mean_power_out(:,1),stderr_power_out(:,1), 'b', ...
    freqs,mean_power_out(:,2),stderr_power_out(:,2), 'g', ...
    freqs,mean_power_out(:,3),stderr_power_out(:,3), 'r', ...
    freqs,mean_power_out(:,4),stderr_power_out(:,4), 'k'); 
axis tight
xlim([0 30])

xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title('Entire trial power Standard Trials Electrode Pz');
legend(conds_lab,'Location','NorthEast'); 
