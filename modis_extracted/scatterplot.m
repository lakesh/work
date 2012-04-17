modisFolderPath = '/work/modis_extracted/2003/';
aeronetFolderPath = '/work/aeronet/2003/';
folders = dir(modisFolderPath);
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
            disp(folderName);
            try 
                %load([modisFolderPath folderName '/interpolatedData_' folderName '.mat']);
                load([modisFolderPath folderName '/interpolated_' folderName '.mat']);
                load([aeronetFolderPath folderName '.mat']);
                [row1 column1, layer1] = size(data);
                [row2 column2] = size(aeronetdata);
                
                for i=1:row1
                    for k=1:column1
                        dataPointMODIS = data(i,k,:);
                        for j=1:row2
                            dataPointAERONET = aeronetdata(j,:);
                            distance = calculateDistance(dataPointMODIS(1,1,5),dataPointMODIS(1,1,6),dataPointAERONET(1,17),dataPointAERONET(1,18));
                            leftBoundary = dataPointMODIS(1,1,6) - deltaX/2;
                            %rightBoundary = dataPointMODIS(1,1,6) + deltaX/2;
                            bottomBoundary = dataPointMODIS(1,1,5) - deltaY/2;
                            %topBoundary = dataPointMODIS(1,1,5) + deltaY/2;
                            
                            if k == column1
                                topBoundary = actualTopBoundary;
                            else
                                topBoundary = dataPointMODIS(1,1,5) + deltaY/2;
                            end
                            
                            if i == row1
                                rightBoundary = actualRightBoundary;
                            else
                                rightBoundary = dataPointMODIS(1,1,6) + deltaX/2;
                            end
                            
                            %Apply distance threshold for MODIS and AERONET
                            %data point
                            %if distance <= 30
                            if dataPointAERONET(1,17) >= bottomBoundary && dataPointAERONET(1,17) <= topBoundary && dataPointAERONET(1,18) >= leftBoundary && dataPointAERONET(1,18) <= rightBoundary
                                hourMODIS = dataPointMODIS(1,1,3);
                                minuteMODIS = dataPointMODIS(1,1,4);
                                totalMinutesMODIS = hourMODIS*60 + minuteMODIS;
                                hourAERONET = dataPointAERONET(1,4);
                                minuteAERONET = dataPointAERONET(1,5);
                                totalMinutesAERONET = hourAERONET*60 + minuteAERONET;
                                
                                %Apply time threshold for MISR and AERONET
                                %data point
                                if abs(totalMinutesAERONET - totalMinutesMODIS) <= 30
                                    %collocatedData(index,1) = dataPointMODIS(1,1,23);
                                    %MODIS AOD value
                                    collocatedData(index,1) = dataPointMODIS(1,1,23);
                                    %AERONET AOD value
                                    collocatedData(index,2) = dataPointAERONET(1,10);
                                    
                                    %Cloud screening(IModis(xi))
                                    %dataPointMODIS(1,8) > 0 means cloud
                                    %masked
                                    if dataPointMODIS(1,8) > 0 || dataPointMODIS(1,1,23) <= 0 || isnan(dataPointMODIS(1,1,23));
                                        collocatedData(index,3) = 0;
                                        collocatedData(index,1) = 0;
                                    else
                                        collocatedData(index,3) = 1;
                                    end
                                    
                                    %Latitude
                                    collocatedData(index,4) = dataPointMODIS(1,1,5);
                                    %Longitude
                                    collocatedData(index,5) = dataPointMODIS(1,1,6);
                                    %day of  year
                                    collocatedData(index,6) = str2double(folderName);
                                    %hour 
                                    collocatedData(index,7) = dataPointMODIS(1,1,3);
                                    %minute
                                    collocatedData(index,8) = dataPointMODIS(1,1,4);
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

%Remove the -9999 values and nan values for MODIS
index = find(collocatedData(:,1) == -9999);
collocatedData(index,:) = [];

index = find(isnan(collocatedData(:,1)));
collocatedData(index,:) = [];

%Remove the -9999 values and nan values for AERONET
index = find(isnan(collocatedData(:,2)));
collocatedData(index,:) = [];

%Scatter plot
scatter(collocatedData(:,1),collocatedData(:,2));