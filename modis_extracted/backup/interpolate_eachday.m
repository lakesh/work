leftBoundary = -72.1889;
rightBoundary = -70.550;
topBoundary = 43.10889;
bottomBoundary = 41.3;

deltaX = 0.1349;
deltaY = 0.179562;

missingValue = -9999;

[xq yq] = meshgrid(41.3 + deltaX/2:0.1349:43.10889, -72.1889+deltaY/2:0.179562:-70.550);
[row1 column1]=size(xq);

[b m n] = unique(data(:,4));
[row column] = size(b);

for i=1:row
    cd(num2str(b(i,1)));
    load([num2str(b(i,1)) '.mat']);
    interpolatedData = zeros(row1,column1,40);
    
    try 
        dataWithNormalization = data;
        index = find(dataWithNormalization == -9999);
        dataWithNormalization(index) = 1.0;

        %year,dayofyear,hour,minute
        interpolatedData(:,:,1) = data(1,3);
        interpolatedData(:,:,2) = data(1,4);
        interpolatedData(:,:,3) = griddata(data(:,1),data(:,2),data(:,5),xq,yq,'nearest');
        interpolatedData(:,:,4) = griddata(data(:,1),data(:,2),data(:,6),xq,yq,'nearest');

        %latitude,longitude
        interpolatedData(:,:,5) = xq;
        interpolatedData(:,:,6) = yq;

        %Land_Type, CloudMask_Quality, CloudMask_Summary, SnowIce
        interpolatedData(:,:,7) = griddata(data(:,1),data(:,2),data(:,7),xq,yq,'nearest');
        interpolatedData(:,:,8) = griddata(data(:,1),data(:,2),data(:,8),xq,yq,'nearest');
        interpolatedData(:,:,9) = griddata(data(:,1),data(:,2),data(:,9),xq,yq,'nearest');
        interpolatedData(:,:,10) = griddata(data(:,1),data(:,2),data(:,10),xq,yq,'nearest');


        %Aerosol_Type_Land
        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,11) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,11) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,11),xq,yq);
        interpolatedData(:,:,12) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,11),xq,yq);

        %CloudFractionLand
        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,12) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,13) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,12),xq,yq);
        interpolatedData(:,:,14) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,12),xq,yq);

        %Quality_Assurance_Land
        %DB_QA_Usefulness, DB_QA_Conf, DB_Aerosol_Type, DB_Retrieving_Condition
        interpolatedData(:,:,15) = griddata(data(:,1),data(:,2),data(:,13),xq,yq,'nearest');
        interpolatedData(:,:,16) = griddata(data(:,1),data(:,2),data(:,14),xq,yq,'nearest');
        interpolatedData(:,:,17) = griddata(data(:,1),data(:,2),data(:,15),xq,yq,'nearest');
        interpolatedData(:,:,18) = griddata(data(:,1),data(:,2),data(:,16),xq,yq,'nearest');


        %%%%%%%%%%%%%%Corrected_Optical_Depth_Land%%%%%%%%%%%%%
        %AOD_213,AOD_470,AOD_550,AOD_660
        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,17) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,19) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,17),xq,yq);
        interpolatedData(:,:,20) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,17),xq,yq);

        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,18) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,21) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,18),xq,yq);
        interpolatedData(:,:,22) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,18),xq,yq);

        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,19) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,23) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,19),xq,yq);
        interpolatedData(:,:,24) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,19),xq,yq);


        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,20) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,25) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,20),xq,yq);
        interpolatedData(:,:,26) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,20),xq,yq);


        %%%%%%%%%%%%%Angstrom_Exponent_Land %%%%%%%%%%%%%%%%
        %Angstrom_C005
        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,21) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,27) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,21),xq,yq);
        interpolatedData(:,:,28) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,21),xq,yq);

        %%%%%%%%Number_Pixels_Used_Land%%%%%%%%%%%%%%%%
        %PixelsUsed470,PixelsUsed660
        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,22) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,29) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,22),xq,yq);
        interpolatedData(:,:,30) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,22),xq,yq);

        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,23) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,31) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,23),xq,yq);
        interpolatedData(:,:,32) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,23),xq,yq);


         %%%%%%%%%%%%%% Solar Zenith, Solar Azimuth, Sensor Azimuth, Scattering
         %%%%%%%%%%%%%% Angle %%%%%%%%%%%%%%%%
        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,24) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,33) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,24),xq,yq);
        interpolatedData(:,:,34) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,24),xq,yq);

        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,25) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,35) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,25),xq,yq);
        interpolatedData(:,:,36) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,25),xq,yq);

        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,26) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,37) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,26),xq,yq);
        interpolatedData(:,:,38) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,26),xq,yq);

        dataWithoutMissingAttributes = data;
        index = find(dataWithoutMissingAttributes(:,27) == -9999);
        dataWithoutMissingAttributes(index,:) = [];
        interpolatedData(:,:,39) = griddata(dataWithoutMissingAttributes(:,1),dataWithoutMissingAttributes(:,2),dataWithoutMissingAttributes(:,27),xq,yq);
        interpolatedData(:,:,40) = griddata(dataWithNormalization(:,1),dataWithNormalization(:,2),dataWithNormalization(:,27),xq,yq);
        
        data = interpolatedData;
        
        save(['interpolatedData_' num2str(b(i,1)) '.mat'],'data');
        clear interpolatedData data
    catch err
        
    end
    
    cd ..
end