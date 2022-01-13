FC = Freq_Cal('F0',50,'Fs',50);
FC = FC.calcDelayTime;
FC.calcEffResolution(1);
%FC.calcFreqRng
%FC.calcInterHarm
%FC.calcDynMeasRange