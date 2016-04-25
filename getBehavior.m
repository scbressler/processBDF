function SCORE = getBehavior(matFilename)
% SCORE = getBehavior(matFilename)
%   This function is designed specifically for the "Auditory Selective
%   Attention Task" and "Visual/FFR--Visual Selective Attention Task with
%   simultaneous envelope FFR"
%
%   In both of these experiments, each experimental session saves a .mat
%   file with the experimental parameters and subject's behavioral data.
%   This function summarizes this data and returns a MATLAB structure
%   called 'SCORE'
%
%   The function only takes one input, which is the full path and filename
%   of the .mat data file.  If no input is given, the function will pull up
%   a GUI and ask the user to select which file to process.
%
%   DESCRIPTION OF BEHAVIORAL TASK:
%       Subject instructed to either attend to one of three simultaneously
%       occuring stimuli, or ignore the trial altogether.  Attended stimuli
%       are either LEFT or RIGHT, center stimulus is always ignored.
%       Attended stimuli are classified as LEADING or LAGGING.  Subject are
%       tasked with identifying the contour of the changing patterns as:
%           1 = RISING
%           2 = FALLING
%           3 = ZIG-ZAGGING
%       Responses (resp) are entered using the numeric keypad and are check
%       against the answers (answ).
%
%       In cases where subjects are instructed to ignore the trial,
%       subjects are asked to withold their response.  Trials where a
%       response is given are marked as a FALSE ALARM.
%
%       As of 25-Apr-2016, the current experimental design always presents
%       the LEADING stimulus on the LEFT, and the LAGGING stimulus on the
%       RIGHT.  (Note, this parameter (stimIsLR) can be changed on the
%       experimental script, which is why there are options to view 
%       Lag/Left and Lead/Right data.)
%
% Revision history:
%   2016-04-25: created by Scott Bressler (SCB) 

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



