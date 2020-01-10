%CONDITIONS...
%...preferred, clockwise - non-preffered, CCW
%%    
ccc
%
exp = 'Skateboard';
% subs = {'100'	'101'	'102'	'103'	'104'	'106'	'107'	'108'	'109'	'110'...
%     '111'	'112'	'113'	'114'	'115'	'116'	'117'	'119'	'120'...
%     '122'	'123'	'125'	'126'	'127'	'129'	'130'...
%     '131'	'132'	'133'	'134'	'135'	'136'	'137'}; %removed 121 and 128 due to preprocessing issues
%subs = {'100';'101'}; %to test on just one sub 
 %after getting rid of some loud subs
subs = {'100'	'101'	'102'	'103'	'104'	'106'	'107' '108'	'109'	'110'...
          	'114'	'115'	'116'	'117'	'119'	'120' ...
          '122'	'123'		'126'	'127' 	'129'	'130'...
    '131'		'133'	'134'	'135'	'137'};
nsubs = length(subs); 
conds =  {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds_lab = {'Preferred Clockwise'; 'Preferred Counterclockwise'; 'Non-preferred Clockwise'; 'Non-preferred Counterclockwise'};
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\Winter2019/';
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_Corrected_Standard.set'],'filepath','M:\Data\Skateboard\Winter2019\segmentsFFTtest\');
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

left_electrode = 5;
right_electrode = 6;
electrodes = [left_electrode,right_electrode];

wavenumber = 15; %wavelet cycles
freqs = [1:0.5:30]; %wavelet frequencies
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
            for i_electrode = 1:length([left_electrode right_electrode]) %pairs of contralateral electrodes
                tempdata = ALLEEG(i_count).data(electrodes(i_electrode), :,i_trial) ;
                [temp_power temp_times temp_phase] = BOSC_tf(tempdata,freqs,EEG.srate,wavenumber);
                power(:,:,i_electrode,i_trial) = log(temp_power); %take log ?
            end
        end
        power_out(:,:,:,i_sub,i_cond) = nanmean(power,4); %save the power data, averaging across trials
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
                squeeze(power_out(:,:,1,i_sub,i_cond))-...
                squeeze(power_out(:,:,2,i_sub,i_cond)),CLim); %
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

power_spectra_diff = squeeze( power_spectra(:,1,:,:)-...
                            power_spectra(:,2,:,:));

power_spectra_diff_mean = squeeze(mean(power_spectra_diff,2));
power_spectra_diff_se = squeeze(std(power_spectra_diff,[],2)/sqrt(nsubs));

%% Plot spectra accross conditions
close all
figure;
boundedline(freqs,power_spectra_diff_mean(:,1),power_spectra_diff_se(:,1),'k',...
            freqs,power_spectra_diff_mean(:,2),power_spectra_diff_se(:,2),'m',...
            freqs,power_spectra_diff_mean(:,3),power_spectra_diff_se(:,3),'g',...
            freqs,power_spectra_diff_mean(:,4),power_spectra_diff_se(:,4),'b');
axis tight
xlim([0 30])

xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title('Power (Left - Right)');
legend(conds_lab,'Location','NorthEast');
%%
for i_sub = 1:nsubs
     figure; 
boundedline (freqs,mean(power_spectra_diff(:,i_sub,3),3), 0, 'b');
axis tight
xlim([0 30])
ylim([-1 1])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title(subs {i_sub});
end 

%%
close all
for i_sub = 28
figure; 
boundedline (freqs,power_spectra_diff(:,i_sub,1), 0, 'k');
axis tight
xlim([0 30])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title(subs {i_sub});

     figure; 
boundedline (freqs,power_spectra_diff(:,i_sub,2), 0, 'm');
axis tight
xlim([0 30])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title(subs {i_sub});

     figure; 
boundedline (freqs,power_spectra_diff(:,i_sub,3), 0, 'g');
axis tight
xlim([0 30])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title(subs {i_sub});

     figure; 
boundedline (freqs,power_spectra_diff(:,i_sub,4), 0, 'b');
axis tight
xlim([0 30])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title(subs {i_sub});
end




