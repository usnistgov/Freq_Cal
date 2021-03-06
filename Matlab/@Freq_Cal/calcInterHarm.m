function calcInterHarm(obj)
% calculates accuracy (%) of frequency range test data
% accuracy is based on the maximum absolute error for each data file

% this function will create a table with one row per data file, columns will be:
% InterHarmFreq, Freq_accuracy, ROCOF_accuracy
T = table('Size',[numel(obj.dataFiles),3],...
    'VariableTypes',{'double','double','double'},...
    'VariableNames',{'InterHarmFreq','Freq_Accuracy','ROCOF_Accuracy'});

for i = 1:numel(obj.dataFiles)
    C = readcell(cell2mat(obj.dataFiles(i)));
    P = readcell(cell2mat(obj.paramFiles{i}));
    
    % data file frequency
    hdr = string(P(:,1));
    row = find(hdr=='Fh');
    f = cell2mat(P(row,2));
    
    % Calculate the frequency accuracy
    hdr = string(C(1,:));
    col = find(hdr=='PMU_FREQ');
    inst = cell2mat(C(2:end,col));    
    col = find(hdr=='REF_FREQ');
    ref = cell2mat(C(2:end,col));
    fAcc = abs(inst-ref);
    
    % Calculate the ROCOF accuracy
    col = find(hdr=='PMU_ROCOF'); 
    inst = cell2mat(C(2:end,col));    
    col = find(hdr=='REF_ROCOF'); 
    ref = cell2mat(C(2:end,col));    
    rAcc = abs(inst-ref);    
    
    T(i,:) = {f(1), min(fAcc),min(rAcc)};
    
    
end
T = sortrows(T);    % sort the table in order of frequency
writetable(T,fullfile(obj.resultPath,'FreqRng_Results.csv'))

% plot the frequency and ROCOF accuracy
obj.fig = obj.fig+1;
figure(obj.fig)
subplot(2,1,1)
xData = T{:,'InterHarmFreq'}; yData = T{:,'Freq_Accuracy'};
plot(xData,yData,'*k')
title('Frequency Absolute Error (Hz)')
xlabel('Interfering signal frequency (Hz)')
ylabel('Absolute Error (Hz)')
yline(obj.MaxAbsFreqError,'--r')
yl = ylim;
ylim([-0.0005,yl(2)+0.001])
xlim([0,110])
xline(25,'--r')
xline(75,'--r')
xline(45,'--b')
xline(55,'--b')
%yline(99.9900,'r')

subplot(2,1,2)
yData = T{:,'ROCOF_Accuracy'};
plot(xData,yData,'*k')
title('ROCOF Absolute Error (Hz/s)')
xlabel('Interfering signal frequency (Hz)')
ylabel('Absolute Error (Hz/s)')
yline(obj.MaxAbsRocofError,'--r')
yl = ylim;
if yl(2) < obj.MaxAbsRocofError, yl(2) = obj.MaxAbsRocofError; end
ylim([-0.005,yl(2)+0.01])
xlim([0,110])
xline(25,'--r')
xline(75,'--r')
xline(45,'--b')
xline(55,'--b')

end
