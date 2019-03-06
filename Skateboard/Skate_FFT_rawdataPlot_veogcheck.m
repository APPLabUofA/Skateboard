clear all
close all
ccc

exp = 'Skateboard';
subs = {'100';'101'; '102'};

nsubs = length(subs);
conds = {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};%preferred, clockwise - non-preffered, CCW

nconds = length(conds);
Pathname = 'M:\Data\Skateboard\';


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
i_count = 0;
figure;
for i_sub = 1:nsubs
    for i_cond = 1:nconds
        i_count = i_count +1;
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond} '.vhdr'];

        
        EEG = pop_loadbv(Pathname, Filename);
          
        
        subplot(nsubs,nconds,i_count);
        
        
        for i_chan = 15:18
            hold on;
            plot(EEG.times,EEG.data(i_chan,:)-mean(EEG.data(i_chan,:)));
        end
    end
end
legend({EEG.chanlocs(15:18).labels})