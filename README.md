# processBDF
Pre-processing script for raw multi-channel .bdf EEG data.  Outputs filtered, trial-epoched data.

Opens and reads .bdf files using BioSemi-derived scripts
Filters continuous multi-channel EEG data using zero-phase FIR bandpass filter
Eyeblink removal using SSP [optional]
Re-reference to average of left and right mastoid/earlobe signal
Baseline correction [optional]
Trial epoching to event trigger(s)

Outputs MATLAB structure:
  ERP.subjID = subjID;
  ERP.erp = erp; % 3-D data matrix (samples x channels x epochs)
  ERP.triggers = epTrigs; % event trigger vector
  ERP.diodes = epDiodes;  % diode event data 2-D matrix (samples x diode)
  ERP.t = t; % epoch-related time vector (seconds)
  ERP.preStimTime = tPre; % pre-stimulus period
  ERP.postStimTime = tPost; % post-stimulus window
  ERP.fs = fs; % data sampling frequency
  ERP.eyeblinkRej = 0 or 1;
  ERP.eyeblinkRejThreshold = EOGthresh;
  ERP.filterObj = Hd; % zero-phase FIR filter object
  ERP.filterCutOffs = [cf1,cf2]; % filter cutoffs
  ERP.params = results.parameters; % summary of experimental parameters
  
