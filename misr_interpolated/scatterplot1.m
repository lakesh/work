%Collocate MISR and AERONET readings for day 2 of year 2007
%load('/work/misr_interpolated/2007/2/2.mat');
%load('/work/aeronet/2007/2.mat');

misrFolderPath = '/work/misr_interpolated/2007/';
aeronetFolderPath = '/work/aeronet/2007/';
folders = dir(misrFolderPath);
collocatedData = zeros(10000,2);
index = 1;

for counter=1:length(folders)
    folderName = folders(counter).name;
    try
        if folderName(1:1) ~= '.'
            %cd(folderName);
            try 
                load([folderName '/interpolatedData_' folderName '.mat']);
                load([aeronetFolderPath folderName '.mat']);
                [row1 column1, layer1] = size(data);
                [row2 column2] = size(aeronetdata);
                
                for i=1:row1
                    for k=1:column1
                        dataPointMISR = data(i,k,:);
                        for j=1:row2
                            dataPointAERONET = aeronetdata(j,:);
                            distance = calculateDistance(dataPointMISR(1,1,6),dataPointMISR(1,1,7),dataPointAERONET(1,17),dataPointAERONET(1,18));
                            disp(distance);
                            if distance <= 30
                                
                                hourMISR = dataPointMISR(1,1,4);
                                minuteMISR = dataPointMISR(1,1,5);
                                totalMinutesMISR = hourMISR*60 + minuteMISR;
                                hourAERONET = dataPointAERONET(1,4);
                                minuteAERONET = dataPointAERONET(1,5);
                                totalMinutesAERONET = hourAERONET*60 + minuteAERONET;
                                if abs(totalMinutesAERONET - totalMinutesMISR) <= 30
                                    collocatedData(index,1) = dataPointMISR(1,8);
                                    collocatedData(index,2) = dataPointAERONET(1,10);
                                    index = index+1;
                                end
                            end
                        end
                    end

                end
            catch err1
                
            end
        end
    catch err
    end
    
end
collocatedData = collocatedData(1:index,:);
scatter(collocatedData(:,1),collocatedData(:,2));

index = find(collocatedData(:,1) == -9999);
collocatedData(index,:) = [];