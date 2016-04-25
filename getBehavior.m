function SCORE = getBehavior(matFilename)
% Make sure EEGLab is not in the path
% edits:
%   2015-04-01: adjusted epoching scheme SCB

if nargin<1
    [fname,pname] = uigetfile('*.mat');
    matFilename = [pname,fname];
end
% Get root file name:
[subjDir,subjID,~] = fileparts(matFilename);
subjID = subjID(1:5); % get Subject ID (convention: first 5 characters)

load(matFilename); % load experimental data file

%% Behavioral Data
results.answ(isnan(results.answ))=0;
answ = results.answ;
resp = squeeze(results.resp(:,1,:));
resp(isnan(resp)) = 0;

%   attend : 1=Lead, 2=Passive, 3=Lag
% stimIsLR : 1=Left, 2=Right
attLdL = results.attend==1 & results.stimIsLR==1;
attLdR = results.attend==1 & results.stimIsLR==2;
attLgL = results.attend==3 & results.stimIsLR==1;
attLgR = results.attend==3 & results.stimIsLR==2;
attPs = results.attend==2;


SCORE.subjID = results.subjID;
SCORE.TotalPC = mean(resp(:)==answ(:));
SCORE.perblockPC = results.score;
SCORE.attend = results.attend;
SCORE.stimIsLR = results.stimIsLR;
SCORE.answ = answ;
SCORE.resp = resp;
SCORE.hits = (resp==answ);
SCORE.FA = (results.attend==2 & resp~=0);
SCORE.NR = (results.attend~=2 & resp==0);
SCORE.FArate = sum(SCORE.FA(:))/sum(SCORE.attend(:)==2);
SCORE.NRrate = sum(SCORE.NR(:))/sum(SCORE.attend(:)~=2);
SCORE.attLeadL = attLdL;
SCORE.attLeadR = attLdR;
SCORE.attLagL = attLgL;
SCORE.attLagR = attLgR;
SCORE.attPassive = attPs;
SCORE.startTime = results.startTime;
SCORE.stopTime = results.stopTime;

end



