function [diodes,tagIt] = cleanupDiodes(diodes)
%
% [diodes,tagIt] = cleanupDiodes(diodes)
% This function addresses the possibility of the light-sensing diodes
% stuttering causing multiple onset events.
%
% INPUT VARIABLE
%   diodes : 
%
% OUTPUT VARIABLES
%   diodes :
%    tagIt :
%
% Created 2016-08-04 Scott Bressler

nEpochs = size(diodes,3);
nTriggers = size(diodes,2);
nSamples = size(diodes,1);

thresh = 1000;

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

        