function [diodes,tagIt] = cleanupDiodes(diodes,thresh)
%
% [diodes,tagIt] = cleanupDiodes(diodes,thresh)
% This function addresses the possibility of the light-sensing diodes
% stuttering causing multiple onset events.  It assumes that if a
% stuttering event occurs it is within a small time window defined by the
% variable 'thresh', which is the window span in samples
%
% INPUT VARIABLE
%   diodes : 2- or 3-dimension vector of on-off diode information
%            [samples x diode# x epochs]
%   thresh : number of samples to look for stuttering diode information
%
% OUTPUT VARIABLES
%   diodes : cleaned up diode information vector
%    tagIt : no idea what this is
%
% Created 2016-08-04 Scott Bressler

if nargin<2
    thresh = 1000; % default number of samples for stutter detection
end

nEpochs = size(diodes,3); % number of epochs
nTriggers = size(diodes,2); % number of diode channels
nSamples = size(diodes,1); % number of samples

for k = 1:nEpochs
    for m = 1:nTriggers
        
        wTrigs = find(diodes(:,m,k));
        n = find(diff(find(diodes(:,m,k)))<=thresh);
        
        if(sum(n)>0)
            diodes(wTrigs(n+1),m,k) = 0;
            tagIt(k,m) = 1;
        end
            
    end
end

        