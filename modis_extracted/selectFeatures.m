year = 2003;

cd(num2str(year));
AODFolderList = dir('./');
curdir=pwd;

for counter=1:length(AODFolderList)
    AODFolderName = AODFolderList(counter).name;
    cd(AODFolderName);
    try 
        load([AODFolderName '.mat']);

        [row column] = size(data);

        newData = zeros(row,27);

        % Latitude, Longitude, year, dayofyear, hour, minute
        newData(:,1:6) = data(:,1:6);


        %%%%%%%Cloud_Mask_QA%%%%%%%%
        %Land_Type, CloudMask_Quality, CloudMask_Summary, SnowIce(nearest neighbor
        %interpolation)
        newData(:,7) = data(:,37);
        newData(:,8) = data(:,23);
        newData(:,9) = data(:,24);
        newData(:,10) = data(:,76);

        %%%%%%%%Aerosol_Type_Land %%%%%%%%%%%
        newData(:,11) = data(:,19);

        %%%%%%%%%%% Mean_Reflectance_Land_All%%%%%%%%
        % Reflectance_After_Cloud_Screening_470,
        % Reflectance_After_Cloud_Screening_660,
        % Reflectance_After_Cloud_Screening_2100
        %newData(:,12) = data(:,56);
        %newData(:,13) = data(:,57);
        %newData(:,14) = data(:,55);


        %%%%%%%%%%%% Cloud_Fraction_Land %%%%%%%%%%%%%%
        newData(:,12) = data(:,22);

        %%%%%%%%%%%% Quality_Assurance_Land %%%%%%%%%%%
        %DB_QA_Usefulness, DB_QA_Conf, DB_Aerosol_Type,
        %DB_Retrieving_Condition(nearest neighbor interpolation)
        newData(:,13) = data(:,29);
        newData(:,14) = data(:,28);
        newData(:,15) = data(:,27);
        newData(:,16) = data(:,30);


        %%%%%%%%%%%%%%Corrected_Optical_Depth_Land%%%%%%%%%%%%%
        %AOD_213,AOD_470,AOD_550,AOD_660
        newData(:,17) = data(:,11);
        newData(:,18) = data(:,12);
        newData(:,19) = data(:,13);
        newData(:,20) = data(:,14);

        %%%%%%%%%%%%%Angstrom_Exponent_Land %%%%%%%%%%%%%%%%
        %Angstrom_C005
        newData(:,21) = data(:,20);

        %%%%%%%%Number_Pixels_Used_Land%%%%%%%%%%%%%%%%
        %PixelsUsed470,PixelsUsed660
        newData(:,22) = data(:,45);
        newData(:,23) = data(:,46);



         %%%%%%%%%%%%%% Solar Zenith, Solar Azimuth, Sensor Azimuth, Scattering
         %%%%%%%%%%%%%% Angle %%%%%%%%%%%%%%%%%
         newData(:,24) = data(:,78);
         newData(:,25) = data(:,77);
         newData(:,26) = data(:,73);
         newData(:,27) = data(:,72);
         
         data = newData;
         save([AODFolderName '.mat'],'data');
         cd(curdir);
    catch
        cd(curdir);
    end
end