%Collocate MISR and AERONET readings for day 2 of year 2007
%load('/work/misr_interpolated/2007/2/2.mat');
%load('/work/aeronet/2007/2.mat');

misrFolderPath = '/work/misr_interpolated/2006/';
aeronetFolderPath = '/work/aeronet/2006/';
folders = dir(misrFolderPath);
collocatedData = zeros(10000,9);
index = 1;

deltaY = 0.1349;
deltaX = 0.179562; 

actualLeftBoundary = -72.18800;
actualRightBoundary = -70.55000;
actualTopBoundary = 43.1089;
actualBottomBoundary = 41.30000;

for counter=1:length(folders)
    folderName = folders(counter).name;
    try
        if folderName(1:1) ~= '.'
            %cd(folderName);
            try 
                %load([misrFolderPath folderName '/interpolatedData_' folderName '.mat']);
                load([misrFolderPath folderName '/interpolated_' folderName '.mat']);
                load([aeronetFolderPath folderName '.mat']);
                [row1 column1, layer1] = size(data);
                [row2 column2] = size(aeronetdata);
                
                for i=1:row1
                    for k=1:column1
                        dataPointMISR = data(i,k,:);
                        for j=1:row2
                            dataPointAERONET = aeronetdata(j,:);
                            distance = calculateDistance(dataPointMISR(1,1,6),dataPointMISR(1,1,7),dataPointAERONET(1,17),dataPointAERONET(1,18));
                            leftBoundary = dataPointMISR(1,1,7) - deltaX/2;
                            bottomBoundary = dataPointMISR(1,1,6) - deltaY/2;
                            
                            
                            if k == column1
                                topBoundary = actualTopBoundary;
                            else
                                topBoundary = dataPointMISR(1,1,6) + deltaY/2;
                            end
                            
                            if i == row1
                                rightBoundary = actualRightBoundary;
                            else
                                rightBoundary = dataPointMISR(1,1,7) + deltaX/2;
                            end
                            
                            %Apply distance threshold for MISR and AERONET
                            %data point
                            %if distance <= 15
                            if dataPointAERONET(1,17) >= bottomBoundary && dataPointAERONET(1,17) <= topBoundary && dataPointAERONET(1,18) >= leftBoundary && dataPointAERONET(1,18) <= rightBoundary
                                
                                hourMISR = dataPointMISR(1,1,4);
                                minuteMISR = dataPointMISR(1,1,5);
                                totalMinutesMISR = hourMISR*60 + minuteMISR;
                                hourAERONET = dataPointAERONET(1,4);
                                minuteAERONET = dataPointAERONET(1,5);
                                totalMinutesAERONET = hourAERONET*60 + minuteAERONET;
                                
                                %Apply time threshold for MISR and AERONET
                                %data point
                                if abs(totalMinutesAERONET - totalMinutesMISR) <= 30
                                    %MISR AOD
                                    %446
                                    %collocatedData(index,1) = dataPointMISR(1,8);
                                    %558
                                    collocatedData(index,1) = dataPointMISR(1,10);
                                    %AERONET AOD
                                    collocatedData(index,2) = dataPointAERONET(1,10);
                                    
                                    
                                    %Cloud screening(IMisr(xi))
                                    if dataPointMISR(1,10) <= 0 || isnan(dataPointMISR(1,10));
                                        collocatedData(index,3) = 0;
                                        %Convert the negative, nan AOD
                                        %values to 0
                                        collocatedData(index,1) = 0;
                                    else
                                        collocatedData(index,3) = 1;
                                    end
                                    
                                    
                                    %Latitude
                                    collocatedData(index,4) = dataPointMISR(1,6);
                                    %Longitude
                                    collocatedData(index,5) = dataPointMISR(1,7);
                                    %Day of year
                                    collocatedData(index,6) = dataPointMISR(1,26);
                                    %hour
                                    collocatedData(index,7) = dataPointMISR(1,4);
                                    %minute
                                    collocatedData(index,8) = dataPointMISR(1,5);
                                    %row
                                    collocatedData(index,9) = i;
                                    %column
                                    collocatedData(index,10) = k;
                                    index = index+1;
                                    break;
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
collocatedData = collocatedData(1:index-1,:);

%Remove the -9999 values
index = find(collocatedData(:,1) == -9999);
collocatedData(index,:) = [];

%Scatter plot
scatter(collocatedData(:,1),collocatedData(:,2));