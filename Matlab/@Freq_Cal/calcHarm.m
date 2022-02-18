function obj = calcHarm(obj)

        obj.getResultsFileList('Choose the folder containing the Harmonic Test Data')
        
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

            T = table('Size',[size(C,1),3],...
                'VariableTypes',{'double','double','double'},...
                'VariableNames',{'Time','Fe','RFe'});
            % get the timestamps
            hdr = string(C(1,:));
            col = find(hdr == 'Timestamp',1,'first');
            t = cell2mat(C(2:end,col));
            t = t-t(1);



end