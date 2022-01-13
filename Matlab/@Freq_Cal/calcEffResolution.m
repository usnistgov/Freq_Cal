function calcEffResolution(obj,idx)

    obj.getResultsFileList('Choose the folder containing the Effective Resolution Test Data')

    % get the ata and parameter info from the indexed files
    C = readcell(cell2mat(obj.dataFiles(idx)));
    P = readcell(cell2mat(obj.paramFiles{idx}));

    % frequency
    hdr = string(C(1,:));
    dt = 1/obj.Fs;
    t = dt;
    col = find(hdr=='FREQ');    % column number
    Y(:,1) = cell2mat(C(2:end,col))';
    col = find(hdr=='ROCOF');
    Y(:,2) = cell2mat(C(2:end,col))';
    
    %---------------- debug: plot the data ------------------------------
    %tVec = (0:1:length(Y)-1)/dt;
    %plot(tVec,Y(:,1),tVec,Y(:,2));
    
    % initial guess of frequency from the parameter file
    w = 2*pi*f;
    
    row = find(hdr=='Ka');
    Xm = cell2mat(P(row,2));
    
    
    % Instantiate a SineFit object 
    %  Use the NIST sineFit class found on github at https://github.com/usnistgov/SineFit
    SF = SineFit('Y',Y,'t',t,'w',w);
    
    % calculate the 4-P fit
    SF.fitter(1);
    
    % plot the fits
    close all
    figure(1)
    subplot(2,1,1)
    hold on
    SF.plot(1,'Y','Y_est','xlim',[0 5])
    hold off
    title(sprintf('Frequency modulation fit, Xm = %.2f, fm = %.2f',Xm,f))
    subplot(2,1,2)
    SF.plot(1,'Y_res','xlim',[0 5]) 
    title(sprintf('Frequency Eff Resolution = %.4f Hz',SF.rmserr(1)/max(SF.Y_est(:,1))))
    
    figure(2)
    subplot(2,1,1)
    hold on
    SF.plot(2,'Y','Y_est','xlim',[0 5])
    hold off
    title(sprintf('ROCOF modulation fit, Xm = %.2f, fm = %.2f',Xm,f))
    subplot(2,1,2)
    SF.plot(2,'Y_res','xlim',[0 5]) 
    title(sprintf('ROCOF Eff. Resolution = %f Hz/s',SF.rmserr(2)/max(SF.Y_est(:,2))))


end