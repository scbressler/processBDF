function ERP = processBDF(bdfFilename,trg,tEpo,cf)
% Make sure EEGLab is not in the path
% edits:
%   2015-04-01: adjusted epoching scheme SCB
%   2016-04-25: generalized version to apply broadly to many EEG
%               applications (for WRNMMC)

if nargin==3
    cf = [1 30];    
elseif nargin==2
    tEpo = input('Enter peristimulus epoch times: ');
    cf = [1 30]; 
elseif nargin==1
    trg = input('Which trigger(s) do you want to epoch? ');
    tEpo = input('Enter peristimulus epoch times (in seconds): ');
    cf = [1 30];     
elseif nargin==0
    [fname,pname] = uigetfile('*.bdf');
    bdfFilename = [pname,fname];
    trg = input('Which trigger(s) do you want to epoch? ');
    tEpo = input('Enter peristimulus epoch times: ');
    cf = [1 30]; 
end

% Flags for additional processing: [0=no, 1=yes(default)]
removeEyeblinks = 1;    % remove eye blink artifacts using SSP
baselineCorrect = 1;    % baseline correction of prestimulus average

% Filter parameters
cf1 = cf(1);    % BPF cutoff frequency in Hz
cf2 = cf(2);   % BPF cutoff frequency in Hz

% Get root file name:
[subjDir,subjID,~] = fileparts(bdfFilename);

% Save data to:
svPath = sprintf('processedDATA/%1.0f-%1.0fHz/',cf1,cf2);

%% FROM .bdf FILES
% Open and read the bdf file
fprintf(['Loading file: ',bdfFilename,'...\n'])
h = openbdf(bdfFilename);
d = readbdf(h,1:h.Head.NRec);

if (h.Head.NS<=5)
    fprintf('4-channel setup...\n')
    channels = [1:2];
elseif (h.Head.NS>5 & h.Head.NS<=41)
    fprintf('32-channel setup...\n');
    channels = [1:32];
elseif (h.Head.NS>41 & h.Head.NS<=73)
    fprintf('64-channel setup...\n');
    channels = [1:64];
end

% Set EXTERNAL ELECTRODES indices: mastoid & EOG
mastRef = length(channels) + [1 2];
EOGch = length(channels) + [3:5];

fs = h.Head.SampleRate(1);
rawBDF = d.Record([channels,mastRef,EOGch],:)';

%% Process
% Get and organize event triggers
triggers = mod(d.Event.triggers,256);
diodes = d.Event.buttons(:,[4 3 2 1]);
diodes = cleanupDiodes(max(0,cat(1,zeros(1,4),diff(diodes))));

OnOff = [0;diff(triggers)];

% Find trigger and diode event samples
di1 = find(diodes(:,1));   
di2 = find(diodes(:,2));
di3 = find(diodes(:,3));
di4 = find(diodes(:,4));

% Find start of epoch trigger(s)
%   Note: you can enter multiple "start of trial" triggers
trCue = []; % initialize variable
for k = 1:length(trg)
    trCue = cat(1,trCue,find(OnOff==trg(k)));
end

trCue = sort(trCue);
trialType = OnOff(trCue);
    
% Bandpass filter continuous raw dataset
nOrder = 2*fs; % number of FIR taps
fprintf('Filtering data (%1.0f to %1.0f Hz)...\n',cf1,cf2);
[data,Hd] = BPF(rawBDF,fs,nOrder,cf1,cf2);

% Re-reference to mastoid
data = data-repmat(mean(data(:,mastRef),2),1,size(data,2));

% Determine if there are multiple recording blocks:
%   You can get discontinuity artifacts that will affect the eyeblink
%   removal processing
%
%   Assumes the "record" triggers = 255
trRec = find(OnOff==255);
for k = 2:length(trRec)
    data(trRec(k)-fix(2*fs):trRec(k)+fix(2*fs),:) = 0;
end

% Zero-out transient response of FIR filter
data(1:fix(2*fs),:) = 0;

% Remove eyeblinks using SSP
%  ref: Uusatilo & Ilmoniemi
if removeEyeblinks
    fprintf('Removing eyeblinks...\n');
    EOG = abs(data(:,1)-data(:,EOGch(1)));
    EOG(1:2*fs) = 0; % zero-out transient response of FIR filter
    [f,x] = ecdf(EOG);
    cdfThresh = x(find(f>0.95,1));
    thresh = 25;
    
    if(cdfThresh>thresh)
        EOGthresh = cdfThresh;
    else
        EOGthresh = thresh;
    end
    
    EOGsamples = EOG>EOGthresh;
    blinkData = data(EOGsamples,:);
    blinkCov = cov(blinkData);
    [E,S] = eig(blinkCov);
    maxEigVec = find(diag(S)==max(diag(S)));
    C = E(:,maxEigVec);
    P = eye(size(data,2)) - C*C';

    data = (P*data')';
end

% Epoch to desired trigger events
tPre = tEpo(1); tPost = tEpo(2);
fprintf('Epoching %1.3f sec to +%1.3f sec...\n',tEpo);
preroll = fix(tPre*fs);
postroll = fix(tPost*fs);
    
for k = 1:length(trCue)
    erp(:,:,k) = data(trCue(k)+preroll:trCue(k)+postroll,:);
    epDiodes(:,:,k) = diodes(trCue(k)+preroll:trCue(k)+postroll,:);
    epTrigs(:,k) = triggers(trCue(k)+preroll:trCue(k)+postroll,:);
end

t = (0:size(erp,1)-1)/fs;
t = t+(tPre);

% Correct for baseline offset
if baselineCorrect
    fprintf('Removing baseline re: first %1.0f ms...\n',-tPre*1000);
    tBL = find(t<0.000);
    erp = erp-repmat(mean(erp(tBL,:,:)),[size(erp,1),1,1]);
end

%% Compile ERP data structure
ERP.subjID = subjID;
ERP.erp = erp; % EEG data [n x channels x trials]
ERP.trialType = full(trialType); % list of trial type by event trigger
ERP.triggers = epTrigs; % list of all within-epoch trigger events
ERP.diodes = epDiodes; % list of all within-epoch diode events
ERP.t = t; % epoch time vector (seconds)
ERP.preStimTime = tPre; % pre-stimulus time range for epoch
ERP.postStimTime = tPost; % post-stimulus time range for epoch
ERP.fs = fs; % data sampling frequency (Hz)
ERP.eyeblinkRej = removeEyeblinks; % [0=no, 1=yes]
ERP.eyeblinkRejThreshold = EOGthresh; % threshold level for eyeblink detection (µV)
ERP.baselineCorrection = baselineCorrect; % [0=no, 1=yes]
ERP.filterObj = Hd; % filter object
ERP.filterCutOffs = [cf1,cf2]; % filter passband cutoff frequencies (Hz)
ERP.dataDimensions = '[n x channels x trials]';

end



