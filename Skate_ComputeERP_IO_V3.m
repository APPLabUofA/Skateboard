%%
ccc
%%

exp = 'Skateboard';
subs = {'100' '101' '102' '103' '104' '106' '107' '108' '109' ...
    '110' '111' '112' '113' '114' '115' '116' '117' '118' '119' ...
    '120' '122' '123' '124' '125' '126' '127' '128' '129'};
is_goofy = [0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1,...
    0, 0, 0, 0, 0, 1, 1];

nsubs = length(subs);
conds =  {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};

nconds = length(conds);
Pathname = 'M:\Data\Skateboard\Winter2019\'; %M:\Data\Skateboard\Winter2019
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_Corrected_Target.set'],'filepath','M:\Data\Skateboard\Winter2019\segments\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = pop_loadset('filename',[Filename '_Corrected_Standard.set'],'filepath','M:\Data\Skateboard\Winter2019\segments\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        
    end
end
eeglab redraw

%%
%subject erps
electrode = 15;%this is PZ in this electrode map
erp_out = [];
for i_sub = 1:nsubs
    fprintf(['Subject ' num2str(i_sub)])
    for i_cond = 1:nconds
        %average over trials (3rd dimension)
        erp_out(:,1,:,i_cond,i_sub) = mean(ALLEEG(1+ 2*((i_sub-1)*nconds+(i_cond-1))).data,3)'; %Targets
        erp_out(:,2,:,i_cond,i_sub) = mean(ALLEEG(2+ 2*((i_sub-1)*nconds+(i_cond-1))).data,3)'; %standards
        ALLEEG(1+ 2*((i_sub-1)*nconds+(i_cond-1))).setname
        ALLEEG(2+ 2*((i_sub-1)*nconds+(i_cond-1))).setname
        
    end
end

%%
%generating new variables to plot them
%erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));


% if is_goofy(i_sub)% this is to separate variables by goofy vs not goofy. 
%     erp_P_in_targ{1} = squeeze(erp_out(:,1,:,2,logical(is_goofy))); %{'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
%     erp_P_out_targ{1} = squeeze(erp_out(:,1,:,1,logical(is_goofy)));
%     erp_NP_in_targ{1} = squeeze(erp_out(:,1,:,3,logical(is_goofy)));
%     erp_NP_out_targ{1} = squeeze(erp_out(:,1,:,4,logical(is_goofy)));
%     
%     erp_P_in_stand{1} = squeeze(erp_out(:,2,:,2,logical(is_goofy))); %{'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
%     erp_P_out_stand{1} = squeeze(erp_out(:,2,:,1,logical(is_goofy)));
%     erp_NP_in_stand{1} = squeeze(erp_out(:,2,:,3,logical(is_goofy)));
%     erp_NP_out_stand{1} = squeeze(erp_out(:,2,:,4,logical(is_goofy)));
%     
% elseif ~is_goofy(i_sub)
%     erp_P_in_targ{2} = squeeze(erp_out(:,1,:,1,~logical(is_goofy))); %{'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
%     erp_P_out_targ{2} = squeeze(erp_out(:,1,:,2,~logical(is_goofy)));
%     erp_NP_in_targ{2} = squeeze(erp_out(:,1,:,4,~logical(is_goofy)));
%     erp_NP_out_targ{2} = squeeze(erp_out(:,1,:,3,~logical(is_goofy)));
%     
%     erp_P_in_stand{2} = squeeze(erp_out(:,2,:,1,~logical(is_goofy))); %{'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
%     erp_P_out_stand{2} = squeeze(erp_out(:,2,:,2,~logical(is_goofy)));
%     erp_NP_in_stand{2} = squeeze(erp_out(:,2,:,4,~logical(is_goofy)));
%     erp_NP_out_stand{2} = squeeze(erp_out(:,2,:,3,~logical(is_goofy)));
 
for i_sub = 1:nsubs
    if is_goofy(i_sub)
        erp_P_in_targ(:,:,i_sub) = squeeze(erp_out(:,1,:,2,i_sub)); %{'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
        erp_P_out_targ(:,:,i_sub) = squeeze(erp_out(:,1,:,1,i_sub));
        erp_NP_in_targ(:,:,i_sub) = squeeze(erp_out(:,1,:,3,i_sub));
        erp_NP_out_targ(:,:,i_sub) = squeeze(erp_out(:,1,:,4,i_sub));

        erp_P_in_stand(:,:,i_sub) = squeeze(erp_out(:,2,:,2,i_sub)); %{'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
        erp_P_out_stand(:,:,i_sub) = squeeze(erp_out(:,2,:,1,i_sub));
        erp_NP_in_stand(:,:,i_sub) = squeeze(erp_out(:,2,:,3,i_sub));
        erp_NP_out_stand(:,:,i_sub) = squeeze(erp_out(:,2,:,4,i_sub));

    elseif ~is_goofy(i_sub)

        erp_P_in_targ(:,:,i_sub) = squeeze(erp_out(:,1,:,1,i_sub)); %{'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
        erp_P_out_targ(:,:,i_sub) = squeeze(erp_out(:,1,:,2,i_sub));
        erp_NP_in_targ(:,:,i_sub) = squeeze(erp_out(:,1,:,4,i_sub));
        erp_NP_out_targ(:,:,i_sub) = squeeze(erp_out(:,1,:,3,i_sub));

        erp_P_in_stand(:,:,i_sub) = squeeze(erp_out(:,2,:,1,i_sub)); %{'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
        erp_P_out_stand(:,:,i_sub) = squeeze(erp_out(:,2,:,2,i_sub));
        erp_NP_in_stand(:,:,i_sub) = squeeze(erp_out(:,2,:,4,i_sub));
        erp_NP_out_stand(:,:,i_sub) = squeeze(erp_out(:,2,:,3,i_sub));

    end
end
%%
%% PLOTS

%preferred facing IN
figure('Color',[1 1 1]);
subplot (2,2,1);
boundedline(EEG.times,squeeze(mean(erp_P_in_targ(:,electrode,:),3)),squeeze(std(erp_P_in_targ(:,electrode,:),[],3))./sqrt(nsubs),'r',...
    EEG.times,squeeze(mean(erp_P_in_stand(:,electrode,:),3)),squeeze(std(erp_P_in_stand(:,electrode,:),[],3))./sqrt(nsubs),'k');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
legend('Targets','Standards','Location','northwest');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Facing In - Preferred');
xlabel('Time (ms)');
ylabel('Voltage (uV)');

%non-preferred facing IN
subplot (2,2,2);
boundedline(EEG.times,squeeze(mean(erp_NP_in_targ(:,electrode,:),3)),squeeze(std(erp_NP_in_targ(:,electrode,:),[],3))./sqrt(nsubs),'b',...
    EEG.times,squeeze(mean(erp_NP_in_stand(:,electrode,:),3)),squeeze(std(erp_NP_in_stand(:,electrode,:),[],3))./sqrt(nsubs),'k');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
legend('Targets','Standards','Location','northwest');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Facing In - Non-Preferred');
xlabel('Time (ms)');
ylabel('Voltage (uV)');

%preferred facing OUT
% figure('Color',[1 1 1]);
subplot (2,2,3);
boundedline(EEG.times,squeeze(mean(erp_P_out_targ(:,electrode,:),3)),squeeze(std(erp_P_out_targ(:,electrode,:),[],3))./sqrt(nsubs),'r',...
    EEG.times,squeeze(mean(erp_P_out_stand(:,electrode,:),3)),squeeze(std(erp_P_out_stand(:,electrode,:),[],3))./sqrt(nsubs),'k');    set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
legend('Targets', 'Standards','Location','northwest');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Facing Out - Preferred');
xlabel('Time (ms)');
ylabel('Voltage (uV)');

%non-preferred out
subplot (2,2,4);
boundedline(EEG.times,squeeze(mean(erp_NP_out_targ(:,electrode,:),3)),squeeze(std(erp_NP_out_targ(:,electrode,:),[],3))./sqrt(nsubs),'b',...
    EEG.times,squeeze(mean(erp_NP_out_stand(:,electrode,:),3)),squeeze(std(erp_NP_out_stand(:,electrode,:),[],3))./sqrt(nsubs),'k');    set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
legend('Targets', 'Standards','Location','northwest');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Facing Out - Non-Preferred');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%All SUBJECTS side by side plots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%preferred facing IN
figure('Color',[1 1 1]);
for i_sub = 1:nsubs
    hold on
    subplot(6,5,i_sub)
    boundedline(EEG.times,squeeze(mean(erp_P_in_targ(:,electrode,i_sub),3)),squeeze(std(erp_P_in_targ(:,electrode,i_sub),[],3))./sqrt(nsubs),'r',...
        EEG.times,squeeze(mean(erp_P_in_stand(:,electrode,i_sub),3)),squeeze(std(erp_P_in_stand(:,electrode,i_sub),[],3))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    %     legend('Targets','Standards','Location','northwest');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(['pref-in',subs(i_sub)]);
    %xlabel('Time (ms)');
    ylabel('Voltage (uV)');
end

%non-preferred facing IN
figure('Color',[1 1 1]);
for i_sub = 1:nsubs
    hold on
    subplot(6,5,i_sub)
    boundedline(EEG.times,squeeze(mean(erp_NP_in_targ(:,electrode,i_sub),3)),squeeze(std(erp_NP_in_targ(:,electrode,i_sub),[],3))./sqrt(nsubs),'b',...
        EEG.times,squeeze(mean(erp_NP_in_stand(:,electrode,i_sub),3)),squeeze(std(erp_NP_in_stand(:,electrode,i_sub),[],3))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    %legend('Targets','Standards','Location','northwest');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(['nonpref-in',subs(i_sub)]);
    %xlabel('Facing In NonPref - Time (ms)');
    ylabel('Voltage (uV)');
end

%preferred facing OUT
figure('Color',[1 1 1]);
for i_sub = 1:nsubs
    hold on
    subplot(6,5,i_sub)
    boundedline(EEG.times,squeeze(mean(erp_P_out_targ(:,electrode,i_sub),3)),squeeze(std(erp_P_out_targ(:,electrode,i_sub),[],3))./sqrt(nsubs),'r',...
        EEG.times,squeeze(mean(erp_P_out_stand(:,electrode,i_sub),3)),squeeze(std(erp_P_out_stand(:,electrode,i_sub),[],3))./sqrt(nsubs),'k');    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    %legend('Targets', 'Standards','Location','northwest');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(['pref-out',subs(i_sub)]);
    %xlabel('Facing Out Pref - Time (ms)');
    ylabel('Voltage (uV)');
end

%non-preferred out
figure('Color',[1 1 1]);
for i_sub = 1:nsubs
    hold on
    subplot (6,5,i_sub);
    boundedline(EEG.times,squeeze(mean(erp_NP_out_targ(:,electrode,i_sub),3)),squeeze(std(erp_NP_out_targ(:,electrode,i_sub),[],3))./sqrt(nsubs),'b',...
        EEG.times,squeeze(mean(erp_NP_out_stand(:,electrode,i_sub),3)),squeeze(std(erp_NP_out_stand(:,electrode,i_sub),[],3))./sqrt(nsubs),'k');    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    %legend('Targets', 'Standards','Location','northwest');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(['nonpref-out',subs(i_sub)]);
    %xlabel('Facing Out NonPref - Time (ms)');
    ylabel('Voltage (uV)');
end