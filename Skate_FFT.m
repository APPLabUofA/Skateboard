%CONDITIONS...
%...preferred, clockwise - non-preffered, CCW
%%    
ccc
%
exp = 'Skateboard';
subs = {'104' '106' '107' '108' '109' '110' '111'...
    '112' '113' '114' '115' '116' '117' };
%subs = {'100'}; %to test on just one sub 

nsubs = length(subs); 
conds =  {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds_lab = {'Preferred Clockwise'; 'Preferred Counterclockwise'; 'Non-preferred Clockwise'; 'Non-preferred Counterclockwise'};
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

%1:2 -->   F3:F4
%3:4 -->   T7:T8
%5:6 -->   C3:C4
%7:8 -->   P7:P8
%9:10  --> P3:P4
%11:12 --> O1:O2

left_electrode = 5
right_electrode = 6

wavenumber = 12; %wavelet cycles
freqs = [1:1:30]; %wavelet frequencies
power_out = [];
i_count = 0;
n_electrode = EEG.nbchan;
for i_sub = 1:nsubs 
    fprintf(['Subject - ' num2str(i_sub) '. \n']);
    for i_cond = 1:nconds
        fprintf(['Condition - ' conds_lab{i_cond} '. \n']);
        power = [];
        i_count = i_count+1; %which data set in ALLEEG to use
        n_trials = ALLEEG(i_count).trials;
        for i_trial = 1:n_trials
            for i_electrode = [left_electrode right_electrode] %pairs of contralateral electrodes
                tempdata = ALLEEG(i_count).data(i_electrode, :,i_trial) ;
                [temp_power temp_times temp_phase] = BOSC_tf(tempdata,freqs,EEG.srate,wavenumber);
                power(:,:,i_electrode,i_trial) = temp_power; %take log ?
            end
            
        end
        power_out(:,:,:,i_sub,i_cond) = nanmean(power,4); %save the power data
 
        
    end
 
end

%% plot spectrograms for each subject (rows) and condition (columns)
 figure;
 i_count = 0;
 CLim = [-1e6 1e6];
for i_sub = 1:nsubs
    maxval = max(max(max(max(abs(power_out(:,:,:,i_sub,:))))))/10; %computer a colorbar plot range for each subject based 1/10th the max value
    CLim = [-1*maxval maxval];

    for i_cond = 1:nconds
        i_count = i_count+1;
        
        subplot(nconds,nsubs,i_count)
        imagesc(EEG.times,freqs,...
                squeeze(power_out(:,:,left_electrode,i_sub,i_cond))-...
                squeeze(power_out(:,:,right_electrode,i_sub,i_cond)),CLim); %
        colormap(redblue)
        set(gca,'Ydir','normal');
        title(['Sub ' num2str(i_sub) '-' conds_lab{i_cond} '- Left-Right'])
        colorbar
    end
end
xlabel('Time (s)');
ylabel('Frequency (Hz)');



%% Collapse over time to get power spectra, average over subjects


power_spectra = squeeze(mean(power_out,2)); %collapse over time can pick time range here in 2nd dim.
%frequency x electrode x subject x condition

power_spectra_diff = squeeze( power_spectra(:,left_electrode,:,:)-...
                            power_spectra(:,right_electrode,:,:));

power_spectra_diff_mean = squeeze(mean(power_spectra_diff,2));
power_spectra_diff_se = squeeze(std(power_spectra_diff,[],2)/sqrt(nsubs));

%% Plot spectra accross conditions
figure;
boundedline(freqs,power_spectra_diff_mean(:,1),power_spectra_diff_se(:,1),'r',...
            freqs,power_spectra_diff_mean(:,2),power_spectra_diff_se(:,2),'b',...
            freqs,power_spectra_diff_mean(:,3),power_spectra_diff_se(:,3),'g',...
            freqs,power_spectra_diff_mean(:,4),power_spectra_diff_se(:,4)),'r';
axis tight
xlim([0 30])

xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title('Power (Left - Right)');
legend(conds_lab,'Location','NorthEast');





