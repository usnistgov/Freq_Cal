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
    fAcc = calcAcc(inst,ref);
    
    % Calculate the ROCOF accuracy
    col = find(hdr=='PMU_ROCOF'); 
    inst = cell2mat(C(2:end,col));    
    col = find(hdr=='REF_ROCOF'); 
    ref = cell2mat(C(2:end,col));    
    rAcc = (1-abs(inst-ref))*100;    
    
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
title('Frequency Accuracy (%)')
xlabel('Interfering signal frequency (Hz)')
ylabel('Accuracy (%)')
ylim([min(yData)-0.01,100+0.01])
xlim([0,110])
xline(25,'--r')
xline(75,'--r')
xline(45,'--b')
xline(55,'--b')
%yline(99.9900,'r')

subplot(2,1,2)
yData = T{:,'ROCOF_Accuracy'};
plot(xData,yData,'*k')
title('ROCOF Accuracy (%)')
xlabel('Interfering signal frequency (Hz)')
ylabel('Accuracy (%)')
ylim([min(yData)-0.01,100+0.01])
xlim([0,110])
xline(25,'--r')
xline(75,'--r')
xline(45,'--b')
xline(55,'--b')

end


%%-------------------------------------------------------------------------
% local functions
function [accuracy] = calcAcc(inst,ref)
% calculates % accuracy given the instrument value and the reference value
err = inst-ref;
accuracy = ((ref-abs(err))./ref)*100;
end