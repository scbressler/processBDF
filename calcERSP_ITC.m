%calcERSP_ITC.m
clear all; close all;

% Add appropriate EEGLab functions to path
addpath ~/Documents/MATLAB/eeglab13_4_4b/functions/timefreqfunc/;
addpath ~/Documents/MATLAB/eeglab13_4_4b/functions/sigprocfunc/;
addpath ~/Documents/MATLAB/eeglab13_4_4b/functions/adminfunc/;
addpath ~/Documents/MATLAB/eeglab13_4_4b/functions/guifunc/;
addpath ~/Documents/MATLAB/eeglab13_4_4b/functions/popfunc/;
%%
Fs = 4096;
freqs = [6 60]; % wavelet analysis frequency cutoffs
CH = 1:32; % channels to analyze

saveTFData = 0; % Save time-frequency data?
printTFPlots = 0; % Save plots as .eps files?

Ch = [11:20]; % use electrodes PO7 and P08 only (for now...)
wCh = [14 18];
basePath = 'processedDATA/ERP/1-100Hz/';
modality = 'Visual';
N1 = 0.130
YLIM = [0 0.5];

load chanlocs32.mat;
chNames = {chanlocs(CH).labels};

c = dir(strcat(basePath,'*_DATA.mat'));
fnames = {c.name}';
% t = t/1000; % convert time vector from ms to seconds

%% Run it per subject
poolobj = parpool;
for k = 1:length(fnames)

disp(strcat('Loading file: ',fnames{k},'...'));
load(strcat(basePath,fnames{k}));

% Organize all epochs into attend conditions (lead, lag, passive)
erpLagAll = ERP.erp(:,CH,SCORE.attLagR(:));
erpLedAll = ERP.erp(:,CH,SCORE.attLeadL(:));
erpPasAll = ERP.erp(:,CH,SCORE.attPassive(:));

% Artifact rejection
threshold = 100; % remove epochs with any channel (CH) exceeding ±100 µV.

gTr(:,1,k) = squeeze(max(max(abs(erpLedAll(:,Ch,:))),[],2))<threshold;
gTr(:,2,k) = squeeze(max(max(abs(erpLagAll(:,Ch,:))),[],2))<threshold;
gTr(:,3,k) = squeeze(max(max(abs(erpPasAll(:,Ch,:))),[],2))<threshold;

erpLedAll = erpLedAll(:,:,gTr(:,1,k));
erpPasAll = erpPasAll(:,:,gTr(:,2,k));
erpLagAll = erpLagAll(:,:,gTr(:,3,k));

erpLedAll = permute(erpLedAll,[2 1 3]);
erpPasAll = permute(erpPasAll,[2 1 3]);
erpLagAll = permute(erpLagAll,[2 1 3]);

nLed(k) = size(erpLedAll,3);
nPas(k) = size(erpPasAll,3);
nLag(k) = size(erpLagAll,3);

Na(k,:) = [nLed(k),nPas(k),nLag(k)];
N(k) = min(Na(k,:)); % temporarily set to lowest number N for all subjects

% Randomly select epochs for fair comparison
rkLag = randsample(nLag(k),N(k));
rkLed = randsample(nLed(k),N(k));
rkPas = randsample(nPas(k),N(k));

tRange = [ERP.t(1) ERP.t(end)]*1000;
% Run ERSP and ITC analysis on time range, tRange (Line 70)
tic;
parfor ch = 1:length(CH)
    [lagERSP(:,:,ch,k),lagITC(:,:,ch,k),mbase(:,:,ch,k),t(ch,:),f(ch,:)] = ...
     newtimef( erpLagAll(ch,:,rkLag), 7169, tRange, 2048, [3 0.5],...
     'padratio',1,'baseline',[0],'freqs',freqs,...
     'plotmean','off','plotersp','off','plotitc','off','verbose','on');
 
     [ledERSP(:,:,ch,k),ledITC(:,:,ch,k),mbase(:,:,ch,k)] = ...
     newtimef( erpLedAll(ch,:,rkLed), 7169, tRange, 2048, [3 0.5],...
     'padratio',1,'baseline',[0],'freqs',freqs,...
     'plotmean','off','plotersp','off','plotitc','off','verbose','on');
    
     [pasERSP(:,:,ch,k),pasITC(:,:,ch,k),mbase(:,:,ch,k)] = ...
     newtimef( erpPasAll(ch,:,rkPas), 7169, tRange, 2048, [3 0.5],...
     'padratio',1,'baseline',[0],'freqs',freqs,...
     'plotmean','off','plotersp','off','plotitc','off','verbose','on');
end
tRun = toc;
end
delete(poolobj);
t = t(1,:); f = f(1,:); % necessary step if using parfor

%% Generate plots
sLag = 0.75:0.75:3;
sLed = 0.6:0.6:3;
sPas = 0:3;

t13 = [sLag(1:end-1);sLag(1:end-1);sLag(2:end)];
t14 = [sLed(1:end-1);sLed(1:end-1);sLed(2:end)];
t1n = [sPas(1:end-1);sPas(1:end-1);sPas(2:end)];

gN = find(N>=24);

for k = gN%1:length(fnames)
    figure; colormap jet;
    set(gcf,'Position',[66 1 1855 1121]);
    sp(1) = subplot(2,3,1);
    p1 = imagesc(t/1000,f,ledERSP(:,:,wCh(1),k));
    colorbar;
    set(gca,'YDir','normal');
    xlabel('Time(sec)'); ylabel('Freq(Hz)');
    title(strcat(strrep(fnames{k}(1:5),'_',' '),' ERSP:Attend Leading'));
    line(0.1+t13(1:2,:),ylim,'Color','b','LineStyle','--');
    line(0.1+t14(1:2,:),ylim,'Color','r','LineStyle','--');
    line(0.1+t1n(1:2,:),ylim,'Color',0.6*ones(1,3),'LineStyle','--');
    
    sp(4) = subplot(2,3,4);
    p4 = imagesc(t/1000,f,abs(ledITC(:,:,wCh(1),k)));
    colorbar;
    set(gca,'YDir','normal','CLim',[0 0.5]);
    xlabel('Time(sec)'); ylabel('Freq(Hz)');
    title(strcat(strrep(fnames{k}(1:5),'_',' '),' ITC:Attend Leading'));
    line(0.1+t13(1:2,:),ylim,'Color','b','LineStyle','--');
    line(0.1+t14(1:2,:),ylim,'Color','r','LineStyle','--');
    line(0.1+t1n(1:2,:),ylim,'Color',0.6*ones(1,3),'LineStyle','--');
    
    sp(2) = subplot(2,3,2);
    p2 = imagesc(t/1000,f,lagERSP(:,:,wCh(1),k));
    colorbar;
    set(gca,'YDir','normal');
    xlabel('Time(sec)');
    title(strcat(strrep(fnames{k}(1:5),'_',' '),' ERSP:Attend Lagging'));
    line(0.1+t13(1:2,:),ylim,'Color','b','LineStyle','--');
    line(0.1+t14(1:2,:),ylim,'Color','r','LineStyle','--');
    line(0.1+t1n(1:2,:),ylim,'Color',0.6*ones(1,3),'LineStyle','--');
    
    sp(5) = subplot(2,3,5);
    p4 = imagesc(t/1000,f,abs(lagITC(:,:,wCh(1),k)));
    colorbar;
    set(gca,'YDir','normal','CLim',[0 0.5]);
    xlabel('Time(sec)'); 
    title(strcat(strrep(fnames{k}(1:5),'_',' '),'ITC:Attend Lagging'));
    line(0.1+t13(1:2,:),ylim,'Color','b','LineStyle','--');
    line(0.1+t14(1:2,:),ylim,'Color','r','LineStyle','--');
    line(0.1+t1n(1:2,:),ylim,'Color',0.6*ones(1,3),'LineStyle','--');
    
    sp(3) = subplot(2,3,3);
    p2 = imagesc(t/1000,f,pasERSP(:,:,wCh(1),k));
    colorbar;
    set(gca,'YDir','normal');
    xlabel('Time(sec)'); 
    title(strcat(strrep(fnames{k}(1:5),'_',' '),' ERSP:Attend Passive'));
    line(0.1+t13(1:2,:),ylim,'Color','b','LineStyle','--');
    line(0.1+t14(1:2,:),ylim,'Color','r','LineStyle','--');
    line(0.1+t1n(1:2,:),ylim,'Color','w','LineStyle','--');
    
    sp(6) = subplot(2,3,6);
    p4 = imagesc(t/1000,f,abs(pasITC(:,:,wCh(1),k)));
    colorbar;
    set(gca,'YDir','normal','CLim',[0 0.5]);
    xlabel('Time(sec)'); 
    title(strcat(strrep(fnames{k}(1:5),'_',' '),'ITC:Attend Passive'));
    line(0.1+t13(1:2,:),ylim,'Color','b','LineStyle','--');
    line(0.1+t14(1:2,:),ylim,'Color','r','LineStyle','--');
    line(0.1+t1n(1:2,:),ylim,'Color',0.6*ones(1,3),'LineStyle','--');
    
    if printTFPlots
    print(['-f',num2str(k)],'-r600','-depsc',['Figures/TF_plots/',fnames{gN(k)}(1:5),'_AllEpochs_TF']);
    end

end

%% Plot low frequency averaged ITCs
f0 = 1:4;
figure;
nBS = 5000; % number of bootstrap draws
thr = 0.95;
tt = [t/1000,fliplr(t/1000)];
n = 1; %length(gN);

for k = 1:2
    PAS = squeeze(mean(abs(pasITC(f0,:,wCh(k),gN))));
    LED = squeeze(mean(abs(ledITC(f0,:,wCh(k),gN))));
    LAG = squeeze(mean(abs(lagITC(f0,:,wCh(k),gN))));
    
    fprintf('Running NONPARAMETRIC permutation stats on %s\n',chNames{wCh(k)});
    h(1,:) = permutationTest(LED',PAS',nBS,thr);
    h(2,:) = permutationTest(LAG',PAS',nBS,thr);
    h = double(h); h(h==0) = NaN;
    
    subplot(2,1,k)
    pa0 = patch(tt',cat(1,mean(PAS,2)+std(PAS,0,2)/sqrt(n),flipud(mean(PAS,2)-std(PAS,0,2)/sqrt(n))),'k');
    hold on;
    set(pa0,'LineStyle','none','FaceColor',0.8*ones(1,3),'FaceAlpha',0.4); hold on;
    pa1 = patch(tt',cat(1,mean(LAG,2)+std(LAG,0,2)/sqrt(n),flipud(mean(LAG,2)-std(LAG,0,2)/sqrt(n))),'b');
    set(pa1,'LineStyle','none','FaceAlpha',0.2);
    pa2 = patch(tt',cat(1,mean(LED,2)+std(LED,0,2)/sqrt(n),flipud(mean(LED,2)-std(LED,0,2)/sqrt(n))),'r');
    set(pa2,'LineStyle','none','FaceAlpha',0.2);
    
    p02 = plot(t/1000,mean(mean(abs(pasITC(f0,:,wCh(k),gN)),4)),'k');
    hold on;
    p01 = plot(t/1000,mean(mean(abs(ledITC(f0,:,wCh(k),gN)),4)),'r');
    p03 = plot(t/1000,mean(mean(abs(lagITC(f0,:,wCh(k),gN)),4)),'b');
    
    plot(t/1000,0.055*h(1,:),'r','LineWidth',3);
    plot(t/1000,0.045*h(2,:),'b','LineWidth',3);
    
    xlabel('Time(sec)'); ylabel('ITC'); ylim(YLIM);
    plotN100Lines(N1);
    title(sprintf('%s Across-Subject Averaged ITC (n=%1.0f), Channel %s',...
          modality,length(fnames),chNames{wCh(k)}));
    legend([p01 p03 p02],'Attend Lead (left)','Attend Lag (right)','Passive',...
                         'Location','NorthWest');
end
    
%% Plot individual subject low frequency ITCs 
for m = 1:length(fnames)
    figure;
for k = 1:2
    subplot(2,1,k)
    p02 = plot(t/1000,mean(abs(pasITC(f0,:,wCh(k),m))),'Color',0.75*ones(1,3),'LineWidth',2);
    hold on;
    p01 = plot(t/1000,mean(abs(ledITC(f0,:,wCh(k),m))),'r');
    p03 = plot(t/1000,mean(abs(lagITC(f0,:,wCh(k),m))),'b');
    xlabel('Time(sec)'); ylabel('ITC');
    plotN100Lines(0.130);
    title(sprintf('%s Subject %s ITC, Channel %s',...
          modality,fnames{m}(1:5),chNames{wCh(k)}));
    legend([p01 p03 p02],'Attend Lead','Attend Lag','Passive',...
                         'Location','NorthWest');
end
end

%% Plot across-subject averaged ITCs
% mLedITC = mean(abs(ledITC),4);
% mPasITC = mean(abs(pasITC),4);
% mLagITC = mean(abs(lagITC),4);
% 
% figure;
% for k = 1:2
%     subplot(1,2,k);
%     imagesc(t/1000,f,ledERSP(:,:,wCh,gN(k)));
%     colorbar;
%     set(gca,'YDir','normal');
% end

%% Topoplots of low frequency ITCs
mLedITC = squeeze(mean(abs(ledITC(f0,:,:,:))));
mPasITC = squeeze(mean(abs(pasITC(f0,:,:,:))));
mLagITC = squeeze(mean(abs(lagITC(f0,:,:,:))));
clim = [0 0.4];

figure;
for k = 1:3
tk = find(t/1000>=(sLag(k)+N1),1,'first');
spX(k) = subplot(3,4,k);
topoplot(mean(mLagITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('Attend Lag: t=%1.4f',t(tk)/1000))
end
for k = 1:4
tk = find(t/1000>=(sLed(k)+N1),1,'first');
spX(4+k) = subplot(3,4,4+k);
topoplot(mean(mLedITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('Attend Lead: t=%1.4f',t(tk)/1000))
end
for k = 1:3
tk = find(t/1000>=(sPas(k)+N1),1,'first');
spX(8+k) = subplot(3,4,8+k);
topoplot(mean(mPasITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('Passive: t=%1.4f',t(tk)/1000))
end

cb = colorbar;
set(spX(11),'Position',[0.5422 0.1100 0.1566 0.2157]);

% More topoplots
% 6-10 Hz ITC repsonses to onset of LAGGING STIMULUS
figure;
for k = 1:3
tk = find(t/1000>=(sLag(k)+N1),1,'first');
sp1(k) = subplot(3,3,k);
topoplot(mean(mLagITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('LAG>AttLag: t=%1.4f',t(tk)/1000));

sp1(3+k) = subplot(3,3,3+k);
topoplot(mean(mLedITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('LAG>AttLead: t=%1.4f',t(tk)/1000));

sp1(6+k) = subplot(3,3,6+k);
topoplot(mean(mPasITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('LAG>Pass: t=%1.4f',t(tk)/1000));
end

tpPos = get(sp1(6+k),'Position');
colorbar;
set(sp1(6+k),'Position',tpPos);

% 6-10 Hz ITC repsonses to onset of LEADING STIMULUS
figure;
for k = 1:4
tk = find(t/1000>=(sLed(k)+N1),1,'first');
sp2(k) = subplot(3,4,k);
topoplot(mean(mLedITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('LEAD>AttLead: t=%1.4f',t(tk)/1000));

sp2(4+k) = subplot(3,4,4+k);
topoplot(mean(mLagITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('LEAD>AttLag: t=%1.4f',t(tk)/1000));

sp2(8+k) = subplot(3,4,8+k);
topoplot(mean(mPasITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('LEAD>Pass: t=%1.4f',t(tk)/1000));
end
tpPos = get(sp2(8+k),'Position');
colorbar;
set(sp2(8+k),'Position',tpPos);

% 6-10 Hz ITC repsonses to onset of PASSIVE/IGNORED/CENTER STIMULUS
figure;

for k = 1:3
tk = find(t/1000>=(sPas(k)+N1),1,'first');
sp3(k) = subplot(3,3,k);
topoplot(mean(mLagITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('PASS>AttLag: t=%1.4f',t(tk)/1000));

sp3(3+k) = subplot(3,3,3+k);
topoplot(mean(mLedITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('PASS>AttLead: t=%1.4f',t(tk)/1000));
tk = find(t/1000>=(sPas(k)+N1),1,'first');

sp3(6+k) = subplot(3,3,6+k);
topoplot(mean(mPasITC(tk,:,:),3),chanlocs,'maplimits',clim);
title(sprintf('PASS>Pass: t=%1.4f',t(tk)/1000))
end

tpPos = get(sp3(6+k),'Position');
colorbar;
set(sp3(6+k),'Position',tpPos);

%% Save data
if saveTFData
save(sprintf('processedDATA/TFData_%s.mat',modality),...
     'lagERSP','ledERSP','pasERSP','lagITC','ledITC','pasITC',...
     't','f','f0','Na','N','CH','N1','fnames','freqs','gN','gTr');
end
