ccc
exp = 'Skateboard';
subs = {'100'	'101'	'102'	'103'	'104'	'106'	'108'	'109'	'110'...
    '111'	'112'	'113'	'115'	'116'	'117'	'119'	'120'...
    '122'	'123'	'125'	'126'	'127'	'129'	...
    '131'	'132'	'133'	'134'	'136'	'137'};
is_goofy = [0,	0,	0,	0,	0,	0,	0,	1,	1,	0,	0,	1,	...
    1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	1,    1,	0];

nsubs = length(subs);
conds =  {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds_lab = {'Preferred Clockwise'; 'Preferred Counterclockwise'; 'Non-preferred Clockwise'; 'Non-preferred Counterclockwise'};
pref_lab = {'Preferred'; 'Non-Preferred'};
facing_lab = {'Facing Inside'; 'Facing Outside'};
clock_lab = {'Clockwise'; 'Counterclock-Wise'};
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\Winter2019/';
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_fft_Standards.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_fft_JK\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        
    end
end
eeglab redraw

% electrode = 15;
 electrode = [1:16];

wavenumber = 12; %wavelet cycles
freqs = [.1:.5:30]; %wavelet frequencies - made freq bin larger
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
                power(:,:,i_electrode,i_trial) = log(temp_power); %take log ?
            end
        end
        power_out(:,:,:,i_sub,i_cond) = nanmean(power,4); %save the power data, averaging across trials
    end
%     if is_goofy(i_sub) == 1
%         power_out_face_in(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[2,3]),5));
%         power_out_face_out(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,4]),5));
%     elseif is_goofy(i_sub) == 0
%         power_out_face_in(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,4]),5));
%         power_out_face_out(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[2,3]),5));
%     end
    power_out_p(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,2]),5)); %averaged preferred blocks 1:2
    power_out_np(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[3,4]),5)); %averaged block 3:4
    power_out_cw(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[1,3]),5)); %clockwise conditions
    power_out_ccw(:,:,:,i_sub) = squeeze(nanmean(power_out(:,:,:,i_sub,[2,4]),5));%ccw merged 
end
%%
%BASELINE

subs_baseline = {'106'	 '108'	'109' '110'...
    '111'  '113' '115'	'116'	'117'	'118'	'119' '120'...
    '121' '122' '123'	'124'	'125'	'126'	'127' '128' '129' ...
    '131' '132' '133'	'134' '135' '136' '137'};

nsubs_baseline = length(subs_baseline);
conds_baseline = {'alpha_C';'alpha_O'};%preferred, clockwise - non-preffered, CCW
nconds_baseline = length(conds_baseline);

Pathname = 'M:\Data\Skateboard\winter2019\';

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for i_sub = 1:nsubs_baseline
    for i_cond = 1:nconds_baseline
        
        Filename = [subs_baseline{i_sub} '_' exp '_' conds_baseline{i_cond} '.vhdr'];
        setname = Filename(1:end-5)
        
        EEG = pop_loadbv(Pathname, Filename);
        
        % get electrode locations
        EEG=pop_chanedit(EEG, 'load',{'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced' 'filetype' 'autodetect'});
        
        % arithmetically rereference to linked mastoid
        for x=1:EEG.nbchan-2
            EEG.data(x,:) = (EEG.data(x,:)-((EEG.data(EEG.nbchan-2,:))*.5));
        end
        
        %%run each electrode through BOSC%%%
        srate = EEG.srate;
        for i_chan = 1:16
            % [bosc_spectra(:,:,i_sub,i_chan,i_cond)] = log(BOSC_tf(EEG.data(i_chan,1:180*srate),F,srate,wavenum));
            [bosc_spectra_baseline(:,i_sub,i_chan,i_cond)] = mean(log(BOSC_tf(EEG.data(i_chan,1:180*srate),freqs,srate,wavenumber)),2);
            
            
        end
    end
end

%% Collapse over time to get power spectra, average over subjects


power_spectra = squeeze(mean(power_out,2)); %collapse over time can pick time range here in 2nd dim.
%frequency x electrode x subject x condition

% time_window = find(EEG.times>350,1)-1:find(EEG.times>550,1)-2;
%here I can plot the spectra of the frequency I want using a time window
%for pre-post stimulus presentation
% power_spectra = squeeze(mean(power_out(:,time_window,:,:,:),2)); %collapse over time can pick time range here in 2nd dim.

power_spectra_p = squeeze(mean(power_out_p,2));
power_spectra_np = squeeze(mean(power_out_np,2));
power_spectra_face_in = squeeze(mean(power_out_face_in,2));
power_spectra_face_out = squeeze(mean(power_out_face_out,2));
power_spectra_cw = squeeze(mean(power_out_cw,2));
power_spectra_ccw = squeeze(mean(power_out_ccw,2));


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
power_spectra_mean_cw = squeeze(mean(power_spectra_cw,2));
power_spectra_se_cw = squeeze(std(power_spectra_cw,[],2)/sqrt(nsubs));
power_spectra_mean_ccw = squeeze(mean(power_spectra_ccw,2));
power_spectra_se_ccw = squeeze(std(power_spectra_ccw,[],2)/sqrt(nsubs));

%%
%prestimulus spectra 
time_window = find(EEG.times>-1000,1)-1:find(EEG.times>0,1)-2;
power_spectra = squeeze(mean(power_out(:,time_window,:,:,:),2)); %collapse over time can pick time range here in 2nd dim.
power_spectra_p = squeeze(mean(power_out_p(:,time_window,:,:),2)); %collapse over time can pick time range here in 2nd dim.
power_spectra_np = squeeze(mean(power_out_np(:,time_window,:,:),2)); %collapse over time can pick time range here in 2nd dim.
power_spectra_cw = squeeze(mean(power_out_cw(:,time_window,:,:),2)); %collapse over time can pick time range here in 2nd dim.
power_spectra_ccw = squeeze(mean(power_out_ccw(:,time_window,:,:),2)); %collapse over time can pick time range here in 2nd dim.


power_spectra_mean = squeeze(mean(power_spectra,2));
power_spectra_se = squeeze(std(power_spectra,[],2)/sqrt(nsubs));
power_spectra_mean_p = squeeze(mean(power_spectra_p,2));
power_spectra_se_p = squeeze(std(power_spectra_p,[],2)/sqrt(nsubs));
power_spectra_mean_np = squeeze(mean(power_spectra_np,2));
power_spectra_se_np = squeeze(std(power_spectra_np,[],2)/sqrt(nsubs));

power_spectra_mean_cw = squeeze(mean(power_spectra_cw,2));
power_spectra_se_cw = squeeze(std(power_spectra_cw,[],2)/sqrt(nsubs));
power_spectra_mean_ccw = squeeze(mean(power_spectra_ccw,2));
power_spectra_se_ccw = squeeze(std(power_spectra_ccw,[],2)/sqrt(nsubs));

%% Plot spectra accross conditions
close all
figure;
boundedline(freqs,power_spectra_mean(:,1),power_spectra_se(:,1),'k',...
    freqs,power_spectra_mean(:,2),power_spectra_se(:,2),'m',...
    freqs,power_spectra_mean(:,3),power_spectra_se(:,3),'g',...
    freqs,power_spectra_mean(:,4),power_spectra_se(:,4),'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title('Power FZ');
legend(conds_lab,'Location','NorthEast');

%plot spectra by preference
figure;
subplot (1,2,1);
boundedline (freqs, power_spectra_mean_p, power_spectra_se_p,'r',...
    freqs, power_spectra_mean_np, power_spectra_se_np,'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title('Power by Preference');
legend(pref_lab,'Location','NorthEast');
%%
%F3/F4 Preferred 
    
%plot spectra by preference

power_spectra_p_f3 = mean(power_spectra_p(:,1,:),2);
power_spectra_mean_pF3 = squeeze(mean(power_spectra_p_f3,3));
power_spectra_se_pF3 = squeeze(std(power_spectra_p_f3,[],3)/sqrt(nsubs));
power_spectra_p_f4 = mean(power_spectra_p(:,2,:),2);
power_spectra_mean_pF4 = squeeze(mean(power_spectra_p_f4,3));
power_spectra_se_pF4 = squeeze(std(power_spectra_p_f4,[],3)/sqrt(nsubs));
power_spectra_p_P7 = mean(power_spectra_p(:,7,:),2);
power_spectra_mean_pP7 = squeeze(mean(power_spectra_p_P7,3));
power_spectra_se_pP7 = squeeze(std(power_spectra_p_P7,[],3)/sqrt(nsubs));
power_spectra_p_P8 = mean(power_spectra_p(:,8,:),2);
power_spectra_mean_pP8 = squeeze(mean(power_spectra_p_P8,3));
power_spectra_se_pP8 = squeeze(std(power_spectra_p_P8,[],3)/sqrt(nsubs));

power_spectra_np_f3 = mean(power_spectra_np(:,1,:),2);
power_spectra_mean_npF3 = squeeze(mean(power_spectra_np_f3,3));
power_spectra_se_npF3 = squeeze(std(power_spectra_np_f3,[],3)/sqrt(nsubs));
power_spectra_np_f4 = mean(power_spectra_np(:,2,:),2);
power_spectra_mean_npF4 = squeeze(mean(power_spectra_np_f4,3));
power_spectra_se_npF4 = squeeze(std(power_spectra_np_f4,[],3)/sqrt(nsubs));
power_spectra_np_P7 = mean(power_spectra_np(:,7,:),2);
power_spectra_mean_npP7 = squeeze(mean(power_spectra_np_P7,3));
power_spectra_se_npP7 = squeeze(std(power_spectra_np_P7,[],3)/sqrt(nsubs));
power_spectra_np_P8 = mean(power_spectra_np(:,8,:),2);
power_spectra_mean_npP8 = squeeze(mean(power_spectra_np_P8,3));
power_spectra_se_npP8 = squeeze(std(power_spectra_np_P8,[],3)/sqrt(nsubs));
%F3 and F4 preferred
figure;
subplot(2,3,1)
boundedline (freqs, power_spectra_mean_pF3, power_spectra_se_pF3,'b',...
    freqs, power_spectra_mean_pF4, power_spectra_se_pF4,'g');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'F3', 'F4'},'Location','northeast')

%F3 and F4 NON preferred

subplot(2,3,2)
boundedline (freqs, power_spectra_mean_npF3, power_spectra_se_npF3,'r',...
    freqs, power_spectra_mean_npF4, power_spectra_se_npF4,'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'F3', 'F4'},'Location','northeast')

figure
%P7 and P8 PREF
subplot(2,3,4)
boundedline (freqs, power_spectra_mean_pP7,power_spectra_se_pP7,'b',...
    freqs, power_spectra_mean_pP8, power_spectra_se_pP8,'g');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'P7', 'P8'},'Location','northeast')

%P7 and P8 NPREF
subplot(2,3,5)
boundedline (freqs, power_spectra_mean_npP7,power_spectra_se_npP7,'r',...
    freqs, power_spectra_mean_npP8, power_spectra_se_npP8,'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'P7', 'P8'},'Location','northeast')


%  non parametric test
fttest = find(freqs>8,1)-1:find(freqs>12,1)-2;

%F3 and F4 preferred
%[p h stat] = ranksum(power_spectra_mean_pF3 ,power_spectra_mean_pF4,'alpha', 0.05, 'tail','both')
[p h stat] = ranksum(squeeze(mean(mean(power_spectra_p_f3(alpha_bins,:,:),1),2)) ,squeeze(mean(mean(power_spectra_p_f4(alpha_bins,:,:),1),2)),'alpha', 0.05, 'tail','both')
%[p h stat] = ranksum(squeeze(nanmean(nanmean(power_out_p(alpha_bins,:,1,:),1),2)),squeeze(nanmean(nanmean(power_out_p(alpha_bins,:,2,:),1),2)),'alpha', 0.05, 'tail','both') 

%F3 and F4 non-preferred
[p h stat] = ranksum(squeeze(mean(mean(power_spectra_np_f3(alpha_bins,:,:),1),2)) ,squeeze(mean(mean(power_spectra_np_f4(alpha_bins,:,:),1),2)),'alpha', 0.05, 'tail','both')
%P7 and P8 PREF
[p h stat] = ranksum(squeeze(mean(mean(power_spectra_p_P7(alpha_bins,:,:),1),2)) ,squeeze(mean(mean(power_spectra_p_P8(alpha_bins,:,:),1),2)),'alpha', 0.05, 'tail','both')
%P7 and P8 non-PREF
[p h stat] = ranksum(squeeze(mean(mean(power_spectra_np_P7(alpha_bins,:,:),1),2)) ,squeeze(mean(mean(power_spectra_np_P8(alpha_bins,:,:),1),2)),'alpha', 0.05, 'tail','both')

%GRAND PREFERENCE 
[p h stat] = ranksum(squeeze(mean(power_spectra_p_left(alpha_bins,:,:),1)) ,squeeze(mean(power_spectra_p_right(alpha_bins,:,:),1)),'alpha', 0.05, 'tail','both')
% p_lr_dif = mean(squeeze(mean(power_spectra_p_left(alpha_bins,:,:),1)) - squeeze(mean(power_spectra_p_right(alpha_bins,:,:),1)));
[p h stat] = ranksum(squeeze(mean(power_spectra_np_left(alpha_bins,:,:),1)) ,squeeze(mean(power_spectra_np_right(alpha_bins,:,:),1)),'alpha', 0.05, 'tail','both')
% np_lr_dif = mean(squeeze(mean(power_spectra_np_left(alpha_bins,:,:),1)) - squeeze(mean(power_spectra_np_right(alpha_bins,:,:),1)));

p_lr_dif = squeeze(mean(power_spectra_p_left(alpha_bins,:,:),1)) - squeeze(mean(power_spectra_p_right(alpha_bins,:,:),1));
np_lr_dif = squeeze(mean(power_spectra_np_left(alpha_bins,:,:),1)) - squeeze(mean(power_spectra_np_right(alpha_bins,:,:),1));
p_fz_dif = squeeze (mean(power_spectra_p_fz(alpha_bins,:,:),1)) - squeeze (mean(power_spectra_np_fz(alpha_bins,:,:),1))
p_pz_dif = squeeze (mean(power_spectra_p_pz(alpha_bins,:,:),1)) - squeeze (mean(power_spectra_np_pz(alpha_bins,:,:),1))

%TTEST 
[h p ci stat] = ttest(p_lr_dif,0,.05,'both',1) 
[h p ci stat] = ttest(np_lr_dif,0,.05,'both',1) 
[h p ci stat] = ttest(p_fz_dif,0,.05,'both',1) 
[h p ci stat] = ttest(p_pz_dif,0,.05,'both',1) 
mean (p_fz_dif)



%%
%Grand preferred left vs right 
power_spectra_p_left = mean(power_spectra_p(:,[1,3,5,7,9,11],:),2);
power_spectra_p_fz = mean(power_spectra_p(:,[13],:),2);
power_spectra_p_pz = mean(power_spectra_p(:,[15],:),2);

power_spectra_p_left_mean = squeeze(mean(power_spectra_p_left,3));
power_spectra_se_p_left = squeeze(std(power_spectra_p_left,[],3)/sqrt(nsubs));
power_spectra_p_right = mean(power_spectra_p(:,[2,4,6,8,10,12],:),2);
power_spectra_p_right_mean = squeeze(mean(power_spectra_p_right,3));
power_spectra_se_p_right = squeeze(std(power_spectra_p_right,[],3)/sqrt(nsubs));

power_spectra_np_left = mean(power_spectra_np(:,[1,3,5,7,9,11],:),2);
power_spectra_np_fz = mean(power_spectra_np(:,[13],:),2);
power_spectra_np_pz = mean(power_spectra_np(:,[15],:),2);
power_spectra_np_left_mean = squeeze(mean(power_spectra_np_left,3));
power_spectra_se_np_left = squeeze(std(power_spectra_np_left,[],3)/sqrt(nsubs));
power_spectra_np_right = mean(power_spectra_np(:,[2,4,6,8,10,12],:),2);
power_spectra_np_right_mean = squeeze(mean(power_spectra_np_right,3));
power_spectra_se_np_right = squeeze(std(power_spectra_np_right,[],3)/sqrt(nsubs));

%preferred left vs right hemisphere
figure;
subplot(2,3,1)
boundedline (freqs, power_spectra_p_left_mean, power_spectra_se_p_left,'b',...
    freqs, power_spectra_p_right_mean, power_spectra_se_p_right,'g');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'Left', 'Right'},'Location','northeast', 'Autoupdate', 'off')
fill([8;8;12;12],[11;15;15;11],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':')


%non-preferred left vs right hemisphere
subplot(2,3,2)
boundedline (freqs, power_spectra_np_left_mean, power_spectra_se_np_left,'r',...
    freqs, power_spectra_np_right_mean, power_spectra_se_np_right,'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'Left', 'Right'},'Location','northeast', 'Autoupdate', 'off')
fill([8;8;12;12],[11;15;15;11],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':')


%%

%Ipsilateral comparisons Left hemisphere vs right hemisphere CW
power_spectra_cw_left = mean(power_spectra_cw(:,[1,3,5,7,9,11],:),2); 
power_spectra_cw_left_mean = mean(power_spectra_cw_left,3);
power_spectra_cw_left_se = squeeze(std(power_spectra_cw_left,[],3)/sqrt(nsubs));
power_spectra_cw_right = mean(power_spectra_cw(:,[2,4,6,8,10,12],:),2); 
power_spectra_cw_right_mean = mean(power_spectra_cw_right,3);
power_spectra_cw_right_se = squeeze(std(power_spectra_cw_right,[],3)/sqrt(nsubs));

%Ipsilateral comparisons Left hemisphere vs right hemisphere CCW
power_spectra_ccw_left = mean(power_spectra_ccw(:,[1,3,5,7,9,11],:),2); 
power_spectra_ccw_left_mean = mean(power_spectra_ccw_left,3);
power_spectra_ccw_left_se = squeeze(std(power_spectra_ccw_left,[],3)/sqrt(nsubs));
power_spectra_ccw_right = mean(power_spectra_ccw(:,[2,4,6,8,10,12],:),2); 
power_spectra_ccw_right_mean = mean(power_spectra_ccw_right,3);
power_spectra_ccw_right_se = squeeze(std(power_spectra_ccw_right,[],3)/sqrt(nsubs));

%Ipsi clockwise 
figure;
subplot(2,3,1)
boundedline (freqs, power_spectra_cw_left_mean,power_spectra_cw_left_se,'g',...
    freqs, power_spectra_cw_right_mean, power_spectra_cw_right_se,'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'Left', 'Right'},'Location','northeast')

%Ipsi counterclockwise 
subplot(2,3,2)
boundedline (freqs, power_spectra_ccw_left_mean,power_spectra_ccw_left_se,'r',...
    freqs, power_spectra_ccw_right_mean, power_spectra_ccw_right_se,'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'Left', 'Right'},'Location','northeast')

%Parieto-occipital 
power_spectra_cw_left_po = mean(mean(power_spectra_cw(:,[7,9,11],:),2),3); 
power_spectra_cw_left_se_po = squeeze(std(mean(power_spectra_cw(:,[7,9,11],:),2),[],3)/sqrt(nsubs));
power_spectra_cw_right_po = mean(mean(power_spectra_cw(:,[8,10,12],:),2),3); 
power_spectra_cw_right_se_po = squeeze(std(mean(power_spectra_cw(:,[8,10,12],:),2),[],3)/sqrt(nsubs));
power_spectra_ccw_left_po = mean(mean(power_spectra_ccw(:,[7,9,11],:),2),3); 
power_spectra_ccw_left_se_po = squeeze(std(mean(power_spectra_ccw(:,[7,9,11],:),2),[],3)/sqrt(nsubs));
power_spectra_ccw_right_po = mean(mean(power_spectra_ccw(:,[8,10,12],:),2),3); 
power_spectra_ccw_right_se_po = squeeze(std(mean(power_spectra_ccw(:,[8,10,12],:),2),[],3)/sqrt(nsubs));
%Ipsi clockwise 
figure;
subplot(2,3,1)
boundedline (freqs, power_spectra_cw_left_po,power_spectra_cw_left_se_po,'g',...
    freqs, power_spectra_cw_right_po, power_spectra_cw_right_se_po,'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'Left', 'Right'},'Location','northeast', 'Autoupdate', 'off')
fill([8;8;12;12],[11;16;16;11],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':')

%Ipsi counterclockwise 
subplot(2,3,2)
boundedline (freqs, power_spectra_ccw_left_po,power_spectra_ccw_left_se_po,'r',...
    freqs, power_spectra_ccw_right_po, power_spectra_ccw_right_se_po,'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
%title('Power by Preference');
legend({'Left', 'Right'},'Location','northeast', 'Autoupdate', 'off')
fill([8;8;12;12],[11;16;16;11],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':')

%%

    %clock and counterclockwise spectra

figure;
boundedline (freqs, power_spectra_mean_cw, power_spectra_se_cw,'r',...
    freqs, power_spectra_mean_ccw, power_spectra_se_ccw,'b');
axis tight
xlim([0 30])
% ylim ([0 3000000])
xlabel('Frequency (Hz)');
ylabel('Power (uV^2)');
title('Power by Clock Orientation');
legend(clock_lab,'Location','NorthEast');


%%
%Baseline spectra
electrode = 15;
figure;
[hl,hr] = boundedline(...
    freqs,squeeze(mean(bosc_spectra_baseline(:,:,electrode,1),2)),std(bosc_spectra_baseline(:,:,electrode,1),[],2)/sqrt(length(subs_baseline)),'r',...
    freqs,squeeze(mean(bosc_spectra_baseline(:,:,electrode,2),2)),std(bosc_spectra_baseline(:,:,electrode,2),[],2)/sqrt(length(subs_baseline)),'b');
set(hl,'linewidth',3);
set(gca,'FontSize',14,'FontWeight', 'bold','linewidth',3,'box','off','color','none','Layer','Top');
legend({'Closed';'Open'});
legend('boxoff');
ylabel('Power (uV^2)','FontSize', 14,'FontWeight', 'bold');
xlabel('Frequency (Hz)','FontSize', 14,'FontWeight', 'bold');
title(['Spectra: Grand-Average: ' EEG.chanlocs(electrode).labels],'FontSize', 16,'FontWeight', 'bold');
%%
%Topographies

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

elec_locs = 'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced';

%Original conds
%%%cond1%%%
topo_p_cw = squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,1),1),2),4));
topo_p_cw(17:18,1) = NaN;
%%%cond2%%%
topo_p_ccw = squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,2),1),2),4));
topo_p_ccw(17:18,1) = NaN;
%%%cond3%%%
topo_np_cw = squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,3),1),2),4));
topo_np_cw(17:18,1) = NaN;
%%%cond4%%%
topo_np_ccw = squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,4),1),2),4));
topo_np_ccw(17:18,1) = NaN;


%preferred & non-preferred

topopref = squeeze(nanmean(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,[1,2]),1),2),4),5));
topopref(17:18,1) = NaN;
topo_npref = squeeze(nanmean(nanmean(nanmean(nanmean(power_out(alpha_bins,:,:,:,[3,4]),1),2),4),5));
topo_npref(17:18,1) = NaN;


%Clockwise and counterclockwise
topocw = squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,:,:),1),2),4));
topocw(17:18,1) = NaN;
topoccw = squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,:,:),1),2),4));
topoccw(17:18,1) = NaN;

%%
time_window = find(EEG.times>-1000,1)-1:find(EEG.times>0,1)-2;
%Clockwise and counterclockwise prestim 
topocw_ps = squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,time_window,:,:),1),2),4));
topocw_ps(17:18,1) = NaN;
topoccw_ps = squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,time_window,:,:),1),2),4));
topoccw_ps(17:18,1) = NaN;
topopref_ps = squeeze(nanmean(nanmean(nanmean(nanmean(power_out(alpha_bins,time_window,:,:,[1,2]),1),2),4),5));
topopref_ps(17:18,1) = NaN;
topo_npref_ps = squeeze(nanmean(nanmean(nanmean(nanmean(power_out(alpha_bins,time_window,:,:,[3,4]),1),2),4),5));
topo_npref_ps(17:18,1) = NaN;

plot_title = 'Power Topography';
%%%Difference%%%
%topo3 = topo1-topo2;
min_lim = min([topo_p_cw;topo_p_ccw]);
max_lim = max([topo_p_cw;topo_p_ccw]);
% diff_min_lim = min([topo1-topo2]);
% diff_max_lim = max([topo1-topo2]);
%%
%original conds
figure;
subplot(1,4,1);
topoplot(topo_p_cw,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
    set(gca,'Color',[1 1 1]);
title('Frequency Range 8-12 Hz');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0]);
%set(t,'LineWidth',line_width);

subplot(1,4,2);
topoplot(topo_p_ccw,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
    set(gca,'Color',[1 1 1]);
%title(['Preferred CCW : ' plot_title]);
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
%set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0]);
%set(t,'LineWidth',line_width);

subplot(1,4,3);
topoplot(topo_np_cw,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
    set(gca,'Color',[1 1 1]);
%title(['Non-Preferred CW: ' plot_title]);
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0]);
%set(t,'LineWidth',line_width);

subplot(1,4,4);
topoplot(topo_np_ccw,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
    set(gca,'Color',[1 1 1]);
%title(['Non-Preferred CWW: ' plot_title]);
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0]);
%set(t,'LineWidth',line_width);



%%
%preference topographies

figure;
subplot(2,1,1);
topoplot(topopref,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
%set(gca,'FontSize',text_size,'FontWeight', 'bold');
%title(['Preferred Stance: ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0]);
%set(t,'LineWidth',line_width);
%non-preferred
subplot(2,1,2)
topoplot(topo_npref,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
%set(gca,'FontSize',text_size,'FontWeight', 'bold');
%title(['Non-Preferred Stance: ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0]);
%set(t,'LineWidth',line_width);

% PRE STIMULUS 
figure;
subplot(2,1,1);
topoplot(topopref_ps,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
set(gca,'FontSize',text_size,'FontWeight', 'bold');
title(['Preferred Stance: ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0],'FontSize',axis_label_size);
set(t,'LineWidth',line_width);
%non-preferred
subplot(2,1,2)
topoplot(topo_npref_ps,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
set(gca,'FontSize',text_size,'FontWeight', 'bold');
title(['Non-Preferred Stance: ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0],'FontSize',axis_label_size);
set(t,'LineWidth',line_width);

%%
%clockwise
figure;
subplot(1,2,1)
topoplot(topocw,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
set(gca,'FontSize',text_size,'FontWeight', 'bold');
title(['Clockwise: ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0],'FontSize',axis_label_size);
set(t,'LineWidth',line_width);
%counterclockwise
subplot(1,2,2)
topoplot(topoccw,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
set(gca,'FontSize',text_size,'FontWeight', 'bold');
title(['Counterclockwise: ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0],'FontSize',axis_label_size);
set(t,'LineWidth',line_width);

% prestimulus 

figure;
subplot(1,2,1)
topoplot(topocw_ps,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
set(gca,'FontSize',text_size,'FontWeight', 'bold');
title(['Clockwise: ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0],'FontSize',axis_label_size);
set(t,'LineWidth',line_width);
%counterclockwise
subplot(1,2,2)
topoplot(topoccw_ps,elec_locs, 'whitebk','on','plotrad',.6,'maplimits',[min_lim max_lim])
set(gca,'FontSize',text_size,'FontWeight', 'bold');
title(['Counterclockwise: ' plot_title],'FontSize', title_size,'FontWeight', 'bold');
t = colorbar('peer',gca);
ylp = get(get(t,'ylabel'), 'Position');
ext=get(get(t,'ylabel'),'Extent');
set(get(t,'ylabel'),'String', 'Power (uV^2)','Rotation',270,'Position',ylp+[3 0 0],'FontSize',axis_label_size);
set(t,'LineWidth',line_width);

%% Contralateral Pairs Mean Power BAR
delta_bins = find(F >= 0 & F <= 3);
theta_bins = find(F >= 4 & F <= 7);
alpha_bins = find(F >= 8 & F <= 12);
beta_bins = find(F >= 13 & F <= 30);

for i_elec = 1:2:12
    figure;
    for i_cond = 1:4
        subplot(2,2,i_cond);
        barweb([...
            squeeze(nanmean(nanmean(nanmean(power_out(delta_bins,:,i_elec,:,i_cond),1),2),4)),...
            squeeze(nanmean(nanmean(nanmean(power_out(delta_bins,:,i_elec+1,:,i_cond),1),2),4));...
            squeeze(nanmean(nanmean(nanmean(power_out(theta_bins,:,i_elec,:,i_cond),1),2),4)),...
            squeeze(nanmean(nanmean(nanmean(power_out(theta_bins,:,i_elec+1,:,i_cond),1),2),4));...
            squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,i_elec,:,i_cond),1),2),4)),...
            squeeze(nanmean(nanmean(nanmean(power_out(alpha_bins,:,i_elec+1,:,i_cond),1),2),4));...
            squeeze(nanmean(nanmean(nanmean(power_out(beta_bins,:,i_elec,:,i_cond),1),2),4)),...
            squeeze(nanmean(nanmean(nanmean(power_out(beta_bins,:,i_elec+1,:,i_cond),1),2),4))],...
            [std(squeeze(nanmean(nanmean(power_out(delta_bins,:,i_elec,:,i_cond),1),2)))/sqrt(nsubs),...
            std(squeeze(nanmean(nanmean(power_out(delta_bins,:,i_elec+1,:,i_cond),1),2)))/sqrt(nsubs);
            std(squeeze(nanmean(nanmean(power_out(theta_bins,:,i_elec,:,i_cond),1),2)))/sqrt(nsubs),...
            std(squeeze(nanmean(nanmean(power_out(theta_bins,:,i_elec+1,:,i_cond),1),2)))/sqrt(nsubs);...
            std(squeeze(nanmean(nanmean(power_out(alpha_bins,:,i_elec+1,:,i_cond),1),2)))/sqrt(nsubs),...
            std(squeeze(nanmean(nanmean(power_out(alpha_bins,:,i_elec+1,:,i_cond),1),2)))/sqrt(nsubs);...
            std(squeeze(nanmean(nanmean(power_out(beta_bins,:,i_elec,:,i_cond),1),2)))/sqrt(nsubs),...
            std(squeeze(nanmean(nanmean(power_out(beta_bins,:,i_elec+1,:,i_cond),1),2)))/sqrt(nsubs)]);
            %,[5,5;5,5;5,5;5,5]);
        set(gca,'FontSize',12,'FontWeight', 'bold','linewidth',3,'box','off','color','none','Layer','Top',...
            'XTickLabel',{'Delta';'Theta';'Alpha';'Beta'});
        
        title([conds_lab{i_cond}]);
        xlabel('Frequency Bins','FontSize', 12,'FontWeight', 'bold');
        ylabel('Mean Power','FontSize', 12,'FontWeight', 'bold');
        ylim([12 16]);
        legend({EEG.chanlocs(i_elec).labels;EEG.chanlocs(i_elec+1).labels});
        legend('boxoff');
    end
end

power_spectra_p_left = mean(power_spectra_p(:,[1,3,5,7,9,11],:),2);
power_spectra_p_left_mean = squeeze(mean(power_spectra_p_left,3));
power_spectra_se_p_left = squeeze(std(power_spectra_p_left,[],3)/sqrt(nsubs));
power_spectra_p_right = mean(power_spectra_p(:,[2,4,6,8,10,12],:),2);
power_spectra_p_right_mean = squeeze(mean(power_spectra_p_right,3));
power_spectra_se_p_right = squeeze(std(power_spectra_p_right,[],3)/sqrt(nsubs));

power_spectra_np_left = mean(power_spectra_np(:,[1,3,5,7,9,11],:),2);
power_spectra_np_left_mean = squeeze(mean(power_spectra_np_left,3));
power_spectra_se_np_left = squeeze(std(power_spectra_np_left,[],3)/sqrt(nsubs));
power_spectra_np_right = mean(power_spectra_np(:,[2,4,6,8,10,12],:),2);
power_spectra_np_right_mean = squeeze(mean(power_spectra_np_right,3));
power_spectra_se_np_right = squeeze(std(power_spectra_np_right,[],3)/sqrt(nsubs));


%By preference contralateral comparisons LEFT VS RIGHT
p_mean_left =  mean(power_spectra_p_left_mean);
p_mean_right = mean (power_spectra_p_right_mean);
p_sd_left = mean(power_spectra_se_p_left);
p_sd_right = mean(power_spectra_se_p_right);
np_mean_left =  mean(power_spectra_np_left_mean);
np_mean_right = mean (power_spectra_np_right_mean);
np_sd_left = mean(power_spectra_se_np_left);
np_sd_right = mean(power_spectra_se_np_right);

figure;
barweb([p_mean_left,p_mean_right; np_mean_left,np_mean_right],...
        [p_sd_left, p_sd_right;np_sd_left, np_sd_right ])

set(gca,'box','off','color','none','Layer','Top',...
    'XTickLabel',{'Preferred'; 'Non-Preferred'});
%title('Contralateral Alpha Power');
%xlabel('Alpha Power','FontSize', 12,'FontWeight', 'bold');
ylim([12 14]);
ylabel('Mean Power');
legend({'Left', 'Right'});
legend('boxoff');

p_lr_dif = squeeze(mean(power_spectra_p_left(alpha_bins,:,:),1)) - squeeze(mean(power_spectra_p_right(alpha_bins,:,:),1));
np_lr_dif = squeeze(mean(power_spectra_np_left(alpha_bins,:,:),1)) - squeeze(mean(power_spectra_np_right(alpha_bins,:,:),1));
mean (squeeze(mean(power_spectra_p_left(alpha_bins,:,:),1)))
mean (squeeze(mean(power_spectra_p_right(alpha_bins,:,:),1)))
%TTEST 
[h p ci stat] = ttest(p_lr_dif,0,.05,'both',1) 
[h p ci stat] = ttest(np_lr_dif,0,.05,'both',1) 

%Overall lateralization by preference 
[h p ci stat] = ttest(p_lr_dif,np_lr_dif,.05,'both',1) 


%ipsilateral comparisons clockwise vs counterclock-wise 
figure;
    barweb([...
        squeeze(nanmean(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[1,3,5,7,9,11],:),1),2),3),4)),...
        squeeze(nanmean(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[2,4,6,8,10,12],:),1),2),3),4));...
        squeeze(nanmean(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[1,3,5,7,9,11],:),1),2),3),4)),...
        squeeze(nanmean(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[2,4,6,8,10,12],:),1),2),3),4))],...
        [...
        std(squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[1,3,5,7,9,11],:),1),2),3)))/sqrt(nsubs),...
        std(squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[2,4,6,8,10,12],:),1),2),3)))/sqrt(nsubs);...
        std(squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[1,3,5,7,9,11],:),1),2),3)))/sqrt(nsubs),...
        std(squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[2,4,6,8,10,12],:),1),2),3)))/sqrt(nsubs)]);
    %,[5,5;5,5;5,5;5,5]);
    set(gca,'box','off','color','none','Layer','Top',...
        'XTickLabel',{'CW'; 'CCW'});
    
    %title('Ipsilateral Power by Clock Orientation');
    %xlabel('Alpha Power');
    ylabel('Mean Power');
    ylim([12 14]);
    legend({'Left Hemisphere';'Right Hemisphere'});
    legend('boxoff');
    
%Using parietooccipital regions only 
figure;
    barweb([...
        squeeze(nanmean(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[7,9,11],:),1),2),3),4)),...
        squeeze(nanmean(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[8,10,12],:),1),2),3),4));...
        squeeze(nanmean(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[7,9,11],:),1),2),3),4)),...
        squeeze(nanmean(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[8,10,12],:),1),2),3),4))],...
        [...
        std(squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[7,9,11],:),1),2),3)))/sqrt(nsubs),...
        std(squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[8,10,12],:),1),2),3)))/sqrt(nsubs);...
        std(squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[7,9,11],:),1),2),3)))/sqrt(nsubs),...
        std(squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[8,10,12],:),1),2),3)))/sqrt(nsubs)]);
    %,[5,5;5,5;5,5;5,5]);
    set(gca,'box','off','color','none','Layer','Top',...
        'XTickLabel',{'CW'; 'CCW'});
    
    %title('Ipsilateral Power by Clock Orientation');
    %xlabel('Alpha Power');
    ylabel('Mean Power');
    ylim([12.5 14]);
    legend({'Left Hemisphere';'Right Hemisphere'});
    legend('boxoff');
    %%
%ipsilateral CW left vs right hemisphere
[h p ci stat] = ttest(squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[1,3,5,7,9,11],:),1),2),3)),squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[2,4,6,8,10,12],:),1),2),3)),.05,'both',1) 
%ipsilateral CCW left vs right hemisphere
[h p ci stat] = ttest(squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[1,3,5,7,9,11],:),1),2),3)),squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[2,4,6,8,10,12],:),1),2),3)),.05,'both',1) 



%parieto-occipital comparisons 
[h p ci stat] = ttest(squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[7,9,11],:),1),2),3)),squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[8,10,12],:),1),2),3)),.05,'both',1) 
%ipsilateral CCW left vs right hemisphere
[h p ci stat] = ttest(squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[7,9,11],:),1),2),3)),squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[8,10,12],:),1),2),3)),.05,'both',1) 


cw_lr_dif = squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[7,9,11],:),1),2),3)) - squeeze(nanmean(nanmean(nanmean(power_out_cw(alpha_bins,:,[8,10,12],:),1),2),3));
ccw_lr_dif = squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[7,9,11],:),1),2),3)) - squeeze(nanmean(nanmean(nanmean(power_out_ccw(alpha_bins,:,[8,10,12],:),1),2),3));

%TTEST 
[h p ci stat] = ttest(cw_lr_dif,0,.05,'both',1) 
[h p ci stat] = ttest(ccw_lr_dif,0,.05,'both',1)

%clock orientation and lateralization 
[h p ci stat] = ttest(cw_lr_dif,ccw_lr_dif,.05,'both',1) 
