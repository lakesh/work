[row column] = size(data);

newData = zeros(row,40);

% Latitude, Longitude, year, dayofyear, hour, minute
newData(:,1:6) = data(:,1:6);


%%%%%%%Cloud_Mask_QA%%%%%%%%
%Land_Type, CloudMask_Quality, CloudMask_Summary, SnowIce
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
newData(:,12) = data(:,56);
newData(:,13) = data(:,57);
newData(:,14) = data(:,55);


%%%%%%%%%%%% Cloud_Fraction_Land %%%%%%%%%%%%%%
newData(:,15) = data(:,22);

%%%%%%%%%%%% Quality_Assurance_Land %%%%%%%%%%%
%DB_QA_Usefulness, DB_QA_Conf, DB_Aerosol_Type, DB_Retrieving_Condition
newData(:,16) = data(:,29);
newData(:,17) = data(:,28);
newData(:,18) = data(:,27);
newData(:,19) = data(:,30);


%%%%%%%%%%%  DeepBlueMeanReflectanceLand %%%%%%%
% ReflectanceDB_412, ReflectanceDB_470, ReflectanceDB_660
newData(:,20) = data(:,52);
newData(:,21) = data(:,53);
newData(:,22) = data(:,54);


%%%%%%%%%%%  DeepBlueSurfaceReflectanceLand %%%%%%%
% SurfaceReflectance_DB_412, SurfaceReflectanceDB_470,
% SurfaceReflectanceDB_660
newData(:,23) = data(:,79);
newData(:,24) = data(:,80);
newData(:,25) = data(:,81);

%%%%%%%%%%%% DeepBlueNumberOfPixelsUsed %%%%%%%%
% NumberPixelsUsedDB_412,NumberPixelsUsedDB_470,NumberPixelsUsedDB_660
newData(:,26) = data(:,39);
newData(:,27) = data(:,40);
newData(:,28) = data(:,41);

 %%%%%%%%%%%% Deep_Blue_Aerosol_Optical_Depth_550 %%%%%%
 %AOD_DB_550
 newData(:,29) = data(:,17);
 
 
 %%%%%%%%%%%% DeepBlue_Aerosol_Optical_Depth_Land %%%%%%%
 %AOD_DB_412, AOD_DB_470, AOD_DB_660
 newData(:,30) = data(:,15);
 newData(:,31) = data(:,16);
 newData(:,32) = data(:,18);
 
 %%%%%%%%%%%% DeepBlue_Angstrom_Exponent_Land %%%%%%%%
 %Angstrom_DB
 newData(:,33) = data(:,21);
 
 %%%%%%%%%%% Deep_Blue_Scattering_Albedo_Land %%%%%%%%%%%
 %SSA_DB_412, SSA_DB_470, SSA_DB_660 
 
 newData(:,34) = data(:,61);
 newData(:,35) = data(:,62);
 newData(:,36) = data(:,63);
 
 
 %%%%%%%%%%%%%% Solar Zenith, Solar Azimuth, Sensor Azimuth, Scattering
 %%%%%%%%%%%%%% Angle %%%%%%%%%%%%%%%%%
 newData(:,37) = data(:,78);
 newData(:,38) = data(:,77);
 newData(:,39) = data(:,73);
 newData(:,40) = data(:,72);
 
 data=newData;