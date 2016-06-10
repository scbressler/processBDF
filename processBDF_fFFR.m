function FFR = processBDF_fFFR(bdfFilepath)
% Make sure EEGLab is not in the path

if nargin<1
    [fname,pname] = uigetfile('rawDATA/*.bdf');
    bdfFilename = [pname,fname];
end

removeEyeBlink = 1;
saveData = 0; % to save to processedData/ directory

% Get root file name:
[subjDir,subjFileID,~] = fileparts(bdfFilename);
subjID = subjFileID(1:5);

%% FROM .bdf FILES
% Open and read the bdf file\
fprintf('Opening BDF file...\n');
h = openbdf(bdfFilename);
d = readbdf(h,1:h.Head.NRec);

% [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16,... % L
%  17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32,... % R
%  33 34 35 36 37 38 39 40]                            % Externals
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
mastoids = channels(end)+[1 2];
llEOG = channels(end) + [3];

fs = h.Head.SampleRate(1);
rawBDF = d.Record([channels,mastoids,llEOG],:)';

%% Process
% Get and organize event triggers
triggers = mod(d.Record(end,:),256);
OnOff = [0,diff(triggers)];
eventList = OnOff(OnOff>0);

% Highpass filter continuous raw dataset
fprintf('Filtering raw BDF data (90Hz Highpass)...\n');
[data,Hd] = HPF(rawBDF,fs);

% Re-reference to mastoid
data = data-repmat(mean(data(:,mastoids),2),1,size(data,2));

% Remove eyeblinks using SSP
%  ref: Uusatilo & Ilmoniemi
if removeEyeBlink
    EOGthresh = 50;
    EOG = abs(data(:,1)-data(:,llEOG));
    EOGsamples = EOG>EOGthresh;
    blinkData = data(EOGsamples,:);
    blinkCov = cov(blinkData);
    [E,S] = eig(blinkCov);
    maxEigVec = find(diag(S)==max(diag(S)));
    C = E(:,maxEigVec);
    P = eye(size(data,2)) - C*C';

    data = (P*data')';
end

% Epoch
fprintf('Epoching data...\n');
tPre = 0.100;
tPost = 5.100;
preroll = fix(tPre*fs);
postroll = fix(tPost*fs);

TR1 = find(OnOff==32);   %  0 dB, +polarity
TR2 = find(OnOff==33);   %  0 dB, -polarity

TR3 = find(OnOff==34);   % -2 dB, +polarity
TR4 = find(OnOff==35);   % -2 dB, -polarity

TR5 = find(OnOff==36);   % -4 dB, +polarity
TR6 = find(OnOff==37);   % -4 dB, -polarity

numStims = min([length(TR1),length(TR2),length(TR3),length(TR4),length(TR5),length(TR6)]);
fprintf('Minimum number of trial presentations = %1.0f\n',numStims);

for k = 1:numStims
    mod0(:,:,k) = data(TR1(k)-preroll:TR1(k)+postroll,channels);
    mod0(:,:,k+numStims) = data(TR2(k)-preroll:TR2(k)+postroll,channels);
    
    mod2(:,:,k) = data(TR3(k)-preroll:TR3(k)+postroll,channels);
    mod2(:,:,k+numStims) = data(TR4(k)-preroll:TR4(k)+postroll,channels);
    
    mod4(:,:,k) = data(TR5(k)-preroll:TR5(k)+postroll,channels);
    mod4(:,:,k+numStims) = data(TR6(k)-preroll:TR6(k)+postroll,channels);
end

% Partition Long eFFR to shorter epochs containing 10 periods of the
% envelope
nPeriods = 10;
nSamp = 400;
nSamp = nPeriods*fs/102.4;
n = size(mod0,1)-mod(size(mod0,1),nSamp);
nCh = channels(end);

mod0 = reshape(permute(mod0(1:n,:,:),[1 3 2]),[nSamp 2*numStims*n/nSamp nCh]);
mod2 = reshape(permute(mod2(1:n,:,:),[1 3 2]),[nSamp 2*numStims*n/nSamp nCh]);
mod4 = reshape(permute(mod4(1:n,:,:),[1 3 2]),[nSamp 2*numStims*n/nSamp nCh]);


% Permute dimensions to fit with Complex PLV function cplxplv( )
mod0 = permute(mod0,[3 1 2]);
mod2 = permute(mod2,[3 1 2]);
mod4 = permute(mod4,[3 1 2]);


FFR.filename = bdfFilename;
FFR.erp(:,:,:,1) = mod0;
FFR.erp(:,:,:,2) = mod2;
FFR.erp(:,:,:,3) = mod4;
FFR.fs = fs;
FFR.t = (0:size(FFR.erp,2)-1)/fs;
FFR.f = ((0:size(FFR.erp,2)-1)*fs)/size(FFR.erp,2);
FFR.numEpochs = 2*numStims*n/nSamp;
FFR.numPeriodsAnalyzed = nPeriods;
FFR.envFreq = 102.4;
FFR.tfsFreq = 4096;
FFR.modDepth = [0 -2 -4];
FFR.preStimTime = tPre;
FFR.postStimTime = tPost;
FFR.HPF = Hd;
FFR.dimLabel = '(channel x time x epoch x modDepth)';

%% Save data files to new directory
if saveData
    save(strcat('processedDATA/FFR/',subjFileID,'_FFRs.mat'),'FFR');
end



