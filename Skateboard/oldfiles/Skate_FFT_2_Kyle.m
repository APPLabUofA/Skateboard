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
Pathname = 'M:\Data\Skateboard\';
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_Corrected_Standard.set'],'filepath','M:\Data\Skateboard\segmentsFFT\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );


    end
end
eeglab redraw

%%

electrode = 7;% pz


i_count = 0;
for i_sub = 1:nsubs 
   
    for i_cond = 1:nconds
        
        i_count = i_count+1; %which data set in ALLEEG to use
        n_trials = ALLEEG(i_count).trials;
        for i_trial = 1:n_trials 
               test = ALLEEG(i_count).data(:, :,:) ;
               power = [];
               phase = [];
                 for i_pick = 1:pick_trials
                     tempdat = test(electrode,:,i_pick);
                    [power(:,i_pick) phase(:,i_pick) freqs] = kyle_fft(tempdat,EEG.srate,30);
%                     power(2:end,i_pick) = power(2:end,i_pick) - mean(power(2:end,i_pick,1)); %subtract mean spectra
                end
                power_out(:,i_sub,i_cond,i_perm) = mean(power(2:end,:),2);
                
        end
 
    end
end
eeglab redraw

mean_power_out = squeeze(mean(power_out,2));
stderr_power_out = squeeze(std(power_out,[],2))./sqrt(nsubs);

figure; 
boundedline(freqs(2:end),mean_power_out(:,1),stderr_power_out(:,1), 'b', freqs(2:end),mean_power_out(:,2),stderr_power_out(:,2), 'r' ); axis tight


%% Check for significance

%Alpha frequencies: 8-12Hz
subfreqs1 = (mean(power_out(16:24,:,1,:),1))'
subfreqs2 = mean(power_out(16:24,:,2,:),1)'
%two-tailed
[h p ci test] = ttest(subfreqs1,subfreqs2)
mdiff = mean(subfreqs1)-mean(subfreqs2)

%%
%high frequencies: >15Hz
subfreqs1 = mean(power_out(30:end,:,1,:),1)
subfreqs2 = mean(power_out(30:end,:,2,:),1)

[h p ci test] = ttest(subfreqs1,subfreqs2)
mdiff = mean(subfreqs1)-mean(subfreqs2)


%Beta frequencies: >23-30Hz (normally 35max)
subfreqs1 = mean(power_out(46:end,:,1,:),1)
subfreqs2 = mean(power_out(46:end,:,2,:),1)

[h p ci test] = ttest(subfreqs1,subfreqs2)
mdiff = mean(subfreqs1)-mean(subfreqs2)



%delta frequencies: >1-4Hz (normally 35max)
subfreqs1 = mean(power_out(2:8,:,1,:),1)
subfreqs2 = mean(power_out(2:8,:,2,:),1)

[h p ci test] = ttest(subfreqs1,subfreqs2)
mdiff = mean(subfreqs1)-mean(subfreqs2)


%theta frequencies: >4-8Hz (normally 35max)
subfreqs1 = mean(power_out(8:16,:,1,:),1)
subfreqs2 = mean(power_out(8:16,:,2,:),1)

[h p ci test] = ttest(subfreqs1,subfreqs2)
mdiff = mean(subfreqs1)-mean(subfreqs2)

