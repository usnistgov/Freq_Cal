function obj = calcHarm(obj)

        obj.getResultsFileList('Choose the folder containing the Harmonic Test Data')
        
        % Allocate a Cell Array for tables
        T_array = cell(3,numel(obj.dataFiles));
        
        for i = 1:numel(obj.dataFiles)

            
            C = readcell(cell2mat(obj.dataFiles(i)));
            P = readcell(cell2mat(obj.paramFiles{i}));
            
            % get the frequency from the Params file
            hdr = string(P(:,1));
            row = find(hdr == 'Fin');
            T_array{1,i} = cell2mat(P(row,2));

            % get the harmonic phase
            row = 8;
            T_array{2,i} = cell2mat(P(row,2));

            % Allocate a table to be stored in the cell array
            T = table('Size',[size(C,1)-1,3],...
                'VariableTypes',{'double','double','double'},...
                'VariableNames',{'Time','Fe','RFe'});
            
            % get the timestamps
            hdr = string(C(1,:));
            col = find(hdr == 'Timestamp',1,'first');
            t = cell2mat(C(2:end,col));
            T.Time = t-t(1);
             
            % get the Freequency Error
            col = find(hdr == "FE",1,'first');
            T.Fe = cell2mat(C(2:end,col));

            % get the ROCOF Error
            col = find(hdr == "RFE",1,'first');
            T.RFe = cell2mat(C(2:end,col));

            T_array{3,i} = T;
        end
        
        % sort the array by phase then frequency
        T_array = sortrows(T_array',[2,1])';
        
        for k = 1:2
            % Plotting
            fFe = figure(obj.fig); obj.fig = obj.fig+1;
            fRFe = figure(obj.fig); obj.fig = obj.fig+1;
            
            % get the results
            for i = 1:3
                n = (k-1)*3+i;
                f = T_array{1,i*k};
                p = T_array{2,n};
                t = T_array{3,n}.Time;
                Fe = T_array{3,n}.Fe;
                RFe = T_array{3,n}.RFe;
                
                
                % frequency plots
                set(0, 'currentfigure', fFe)
                subplot(3,1,i)
                plot(t,Fe,'*k')
                title(sprintf('Frequency Error, F = %d, H. Phase = %d',f,p))
                ylabel('Frequency Error (Hz)')
                xlabel('Time (s)')
                yline(obj.MaxAbsFreqError,'--r')
                ylim([-0.001,1.1*obj.MaxAbsFreqError])
                xlim([0,5])
                set(gca,'FontSize',12)
                
                set(0, 'currentfigure', fRFe)
                subplot(3,1,i)
                plot(t,RFe,'*k')
                title(sprintf('ROCOF Error, F = %d, H. Phase = %d',f,p))
                xlabel('ROCOF Error (Hz/s)')
                ylabel('Time (s)')
                yline(obj.MaxAbsRocofError,'--r')
                ylim([-0.001,1.1*obj.MaxAbsRocofError])
                xlim([0,5])
                set(gca,'FontSize',12)

                
            end
        end
end