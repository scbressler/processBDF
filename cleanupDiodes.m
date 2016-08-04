function [diodes,tagIt] = cleanupDiodes(diodes)
%
% [diodes,tagIt] = cleanupDiodes(diodes)
%
% Addresses the possibility of unreliable diode information.
% Occassionally, the light-sensing diodes will "stutter" during the course
% of a pulsed light-on event.  This function will search for possible
% stutter events and fix.
%
% INPUT VARIABLE
%   diodes : 2- or 3-D matrix of diode information
%            [samples x diode# x epoch]
%
% OUTPUT VARIABLE
%   diodes : cleaned up diode information matrix
%    tagIt : no idea what this is...yet
%
% Created: 2016-08-04 Scott Bressler

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

end

        