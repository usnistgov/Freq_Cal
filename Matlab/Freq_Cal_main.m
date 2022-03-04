 FC = Freq_Cal('F0',50,'Fs',50);
%FC = FC.calcDelayTime;
%FC.calcEffResolution(1);
%FC.calcFreqSettlingTime; % use this for the positive step
%FC.calcFreqSettlingTime; % use this for the negative step
%FC.calcRocofSettlingTime; % use this for the positive step
%FC.calcRocofSettlingTime; % use this for the negative step
%FC.calcFreqRng
%FC.calcInterHarm
%FC.calcHarm
%FC.calcDynMeasRange
FC.calcStep;
%FC.calcCombinedStep
