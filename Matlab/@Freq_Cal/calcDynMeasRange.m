function calcDynMeasRange(obj)
% calculates the accuracy of the dynamic measuring range test.  This test
% inputs a series of frequency modulated signals with modulation frequency
% (Fm) held constant at the measurement range limit (e.g. for +/- 5 Hz
% range, Fm is held constant at 5 Hz.  Then the modulation index is varied
% from 0.1 to 1 in increments of 0.1.
%
% Get the measuring range data files
prompt = sprintf('Choose the folder containing the Dynamic Measuring Range data');
obj.getResultsFileList(prompt)

[modDepth,maxAbsFE,peakROCOF,maxAbsRFE] = RngData(obj);

figure(obj.fig); obj.fig=obj.fig+1;
subplot(2,1,1)
plot(modDepth,maxAbsFE,'b*')
title ('Dynamic Frequency Range')
xlabel('Delta Frequency (Hz)')
set(gca,'YScale','log');
ylabel('Max Absolute Frequency Error (Hz)')
set(gca,'FontSize',12);
xline(5,'b')
yline(obj.MaxAbsFreqErrorDyn,'--r')
hold on


subplot(2,1,2)
plot(peakROCOF,maxAbsRFE,'b*')
title ('Dynamic ROCOF Range')
xlabel('ROCOF (Hz/s)')
set(gca,'YScale','log')
ylabel('Max Absolute ROCOF Error (Hz/s)')
set(gca,'FontSize',12);
xline(60,'b')
yline(obj.MaxAbsRocofErrorDyn,'--r')

hold on

% Get the measuring range data files
prompt = sprintf('Choose the folder containing the Dynamic Operating Range data');
obj.getResultsFileList(prompt)
[modDepth,maxAbsFE,peakROCOF,maxAbsRFE] = RngData(obj);

%figure(obj.fig); obj.fig=obj.fig+1;
subplot(2,1,1)
plot(modDepth,maxAbsFE,'r*')
title ('Dynamic Frequency Range')
xlabel('Delta Frequency (Hz)')
set(gca,'YScale','log');
ylabel('Max Absolute Frequency Error (Hz)')
set(gca,'FontSize',12);

hold off


subplot(2,1,2)
plot(peakROCOF,maxAbsRFE,'r*')
title ('Dynamic ROCOF Range')
xlabel('ROCOF (Hz/s)')
set(gca,'YScale','log')
ylabel('Max Absolute ROCOF Error (Hz/s)')
set(gca,'FontSize',12);
hold off

end

function [modDepth,maxAbsFE,peakROCOF,maxAbsRFE] = RngData(obj)
modDepth = zeros(numel(obj.dataFiles),1);
peakROCOF = modDepth;
maxAbsFE = zeros(numel(obj.dataFiles),1);
maxAbsRFE = maxAbsFE;

for i = 1:numel(obj.dataFiles)
    C = readcell(cell2mat(obj.dataFiles(i)));
    P = readcell(cell2mat(obj.paramFiles{i}));
    
    % get the Fe and RFe for this data file
    hdr = string(C(1,:));
    %col = hdr=='FE';
    col = find(hdr=='FE');
    maxAbsFE(i) = max(abs(cell2mat(C(2:end,col))));
    col = find(hdr=='RFE');
    maxAbsRFE(i) = max(abs(cell2mat(C(2:end,col))));
    
    
    % get the depth of modulation for this parameter file
    hdr = string(P(:,1));
    row = find(hdr=='Fa');
    Fa = cell2mat(P(row,2));
    row = find(hdr=='Ka');
    Ka = cell2mat(P(row,2));
    modDepth(i) = Fa*Ka;  
    peakROCOF(i) = (Ka*Fa^2)*(2*pi);  
    
end
% sort the files in the order of the modulation depth
[modDepth, idx] = sort(modDepth);
peakROCOF = peakROCOF(idx);
maxAbsFE  = maxAbsFE(idx);
maxAbsRFE = maxAbsRFE(idx);
    
end