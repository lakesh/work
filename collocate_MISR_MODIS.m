%MODIS = load('collocated_MODIS_AERONET.mat');
%MISR = load('collocated_MISR_AERONET.mat');

MODIS = load('collocated_MODIS_AERONET_2006.mat');
MISR = load('collocated_MISR_AERONET_2006.mat');

MODIS = MODIS.collocatedData;
MISR = MISR.collocatedData;

[row column] = size(MISR);
[row1 column1] = size(MODIS);

collocatedData = zeros(row,10);
index = 1;


%TODO
%In some cases the AERONET AOD values for the collocated MISR and MODIS
% points are differnet. Right now I am using the aeronet data from
% collocated MODIS_AERONET data. Resolve this issue


for i=1:row
    MISRdataPoint = MISR(i,:);
    for j=1:row1
        MODISdataPoint = MODIS(j,:);
            % If the location of the points are same and they are of the same
            % day, then they are collocated
            if MISRdataPoint(1,6) == MODISdataPoint(1,6) && MISRdataPoint(1,9) == MODISdataPoint(1,9) && MISRdataPoint(1,10) == MODISdataPoint(1,10)
                %MISR AOD value
                collocatedData(index,1) = MISRdataPoint(1,1);
                
                %IMisr(xi)
                collocatedData(index,2) = MISRdataPoint(1,3);
                
                
                %MODIS AOD value
                collocatedData(index,3) = MODISdataPoint(1,1);
                
                %IModis(xi)
                collocatedData(index,4) = MODISdataPoint(1,3);
                
                %AERONET AOD value
                collocatedData(index,5) = MODISdataPoint(1,2);
                %Latitude
                collocatedData(index,6) = MODISdataPoint(1,4);
                %Longitude
                collocatedData(index,7) = MODISdataPoint(1,5);
                %Day of year
                collocatedData(index,8) = MODISdataPoint(1,6);
                %Row
                collocatedData(index,9) = MODISdataPoint(1,9);
                %Column
                collocatedData(index,10) = MODISdataPoint(1,10);
                index = index + 1;
            end
    end
end

nanValues = find(isnan(collocatedData(:,5)));
collocatedData(nanValues,:) = [];