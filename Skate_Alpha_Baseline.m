% Strategy will be, load in preprocessed epochs -1000 to 2000 ms,
% for each epoch go through each channel, run an FFT on the -1000 to 0 data
%(prestimulus), get the spectra,
% Then subtract the spectra for left and right side of the head
% Then average those over epochs
%

%%
clear all
close all
ccc

exp = 'Skateboard';
subs = {'106'	 '108'	'109' '110'...
         '111'  '113' '115'	'116'	'117'	'118'	'119' '120'...
         '121' '122' '123'	'124'	'125'	'126'	'127' '128' '129' ...
         '131' '132' '133'	'134' '135' '136' '137'};

% subs = {'106'	 '108'	'109' '110' '111' '113' '115'	'116'	'117'	'118'  }; %to test on just one sub

nsubs = length(subs);
conds = {'alpha_C';'alpha_O'};%preferred, clockwise - non-preffered, CCW
nconds = length(conds);

F = 0.1:.5:30;
wavenum = 12;

Pathname = 'M:\Data\Skateboard\winter2019\';
filename = 'BOSC_skate_log.mat';

if ~exist([Pathname 'segmentsFFT\'])
    mkdir([Pathname 'segmentsFFT\']);
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
        
        %%run each electrode through BOSC%%%
        srate = EEG.srate;
        for i_chan = 15
           % [bosc_spectra(:,:,i_sub,i_chan,i_cond)] = log(BOSC_tf(EEG.data(i_chan,1:180*srate),F,srate,wavenum));
             [bosc_spectra(:,i_sub,i_chan,i_cond)] = mean(log(BOSC_tf(EEG.data(i_chan,1:180*srate),F,srate,wavenum)),2);
%             kyle said to get rid of time altogether and average to keep
%             file smaller

        end
    end
end


%%
%%%plot grand-average bosc%%%
electrode = 15;
figure;
[hl,hr] = boundedline(...
    F,squeeze(mean(bosc_spectra(:,:,electrode,1),2)),std(bosc_spectra(:,:,electrode,1),[],2)/sqrt(length(subs)),'r',...
    F,squeeze(mean(bosc_spectra(:,:,electrode,2),2)),std(bosc_spectra(:,:,electrode,2),[],2)/sqrt(length(subs)),'b');
set(hl,'linewidth',3);
set(gca,'FontSize',14,'FontWeight', 'bold','linewidth',3,'box','off','color','none','Layer','Top');
legend({'Closed';'Open'});
legend('boxoff');
ylabel('Power (uV^2)','FontSize', 14,'FontWeight', 'bold');
xlabel('Frequency (Hz)','FontSize', 14,'FontWeight', 'bold');
title(['Spectra: Grand-Average: ' EEG.chanlocs(electrode).labels],'FontSize', 16,'FontWeight', 'bold');

% electrode = 15;
% figure;
% [hl,hr] = boundedline(...
%     F,squeeze(mean(bosc_spectra(:,:,electrode,1),2),3)),std(mean(bosc_spectra(:,:,electrode,1),2),[],3)/sqrt(length(subs)),'r',...
%     F,squeeze(mean(mean(bosc_spectra(:,:,electrode,2),2),3)),std(mean(bosc_spectra(:,:,electrode,2),2),[],3)/sqrt(length(subs)),'b');
% set(hl,'linewidth',3);
% set(gca,'FontSize',14,'FontWeight', 'bold','linewidth',3,'box','off','color','none','Layer','Top');
% legend({'Closed';'Open'});
% legend('boxoff');
% ylabel('Power (uV^2)','FontSize', 14,'FontWeight', 'bold');
% xlabel('Frequency (Hz)','FontSize', 14,'FontWeight', 'bold');
% title(['Spectra: Grand-Average: ' EEG.chanlocs(electrode).labels],'FontSize', 16,'FontWeight', 'bold');
