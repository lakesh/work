%%%%%%%%%%%%%%%%%%% This file is for extracting the MODIS data for our
%%%%%%%%%%%%%%%%%%% selected region from Kosta's data
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%




path='/home/lakesh/Desktop/';
years=[2006;2007];
%years=[2006];

leftBoundary=-72.1880;
rightBoundary=-70.55000;
%rightBoundary=-72.88917;
bottomBoundary=41.3000;
topBoundary=43.108889;
%topBoundary=40.95452;

dataPoints=zeros(1425060,87);
curdir=pwd;

width=203;
height=135;
index = 1;

for y=1:length(years)
    AODFileList = dir([path num2str(years(y)) '/*.mat']);

    for counter=1:length(AODFileList)
        AODFileName = AODFileList(counter).name;
        data = load([path num2str(years(y)) '/' AODFileName]);
        latitudeData=data.Latitude;
        longitudeData=data.Longitude;
        disp(AODFileName);
        name=regexp(AODFileName,'\.','split');
        for i=1:width
            for j=1:height
                try 
                    %If the point is inside our region include it
                    if latitudeData(i,j) > bottomBoundary && latitudeData(i,j) < topBoundary && longitudeData(i,j) < rightBoundary && longitudeData(i,j) > leftBoundary
                        dataPoints(index,1) = data.Latitude(i,j);
                        dataPoints(index,2) = data.Longitude(i,j);
                        dataPoints(index,3) = data.AODSmall_213(i,j); %7 // This reference is for the new data file that is created 
                        %%%% afterwards where the columns are arranged such that the
                        %%%% first 6 columns give
                        %%%% latitude,longitude,year,day of year,hour and
                        %%%% minute information. So the 7th column gives
                        %%%% AODSmall_213 and so on
                        dataPoints(index,4) = data.AODSmall_470(i,j); %8
                        dataPoints(index,5) = data.AODSmall_550(i,j); %9 
                        dataPoints(index,6) = data.AODSmall_660(i,j); %10
                        dataPoints(index,7)= data.AOD_213(i,j); %11
                        dataPoints(index,8)= data.AOD_470(i,j); %12
                        dataPoints(index,9)= data.AOD_550(i,j); %13
                        dataPoints(index,10)= data.AOD_660(i,j); %14
                        dataPoints(index,11)= data.AOD_DB_412(i,j); %15
                        dataPoints(index,12)= data.AOD_DB_470(i,j); %16
                        dataPoints(index,13)= data.AOD_DB_550(i,j); %17
                        dataPoints(index,14)= data.AOD_DB_660(i,j); %18
                        dataPoints(index,15)= data.AerosolType(i,j); %19
                        dataPoints(index,16)= data.Angstrom_C005(i,j); %20
                        dataPoints(index,17)= data.Angstrom_DB(i,j); %21
                        dataPoints(index,18)= data.CloudFraction(i,j); %22
                        dataPoints(index,19)= data.CloudMask_Quality(i,j); %23
                        dataPoints(index,20)= data.CloudMask_Summary(i,j); %24
                        dataPoints(index,21)= data.CriticalReflectance_470(i,j); %25
                        dataPoints(index,22)= data.CriticalReflectance_660(i , j); %26
                        dataPoints(index,23)= data.DB_Aerosol_Type(i,j); %27
                        dataPoints(index,24)= data.DB_QA_Conf(i,j); %28
                        dataPoints(index,25)= data.DB_QA_Usefulness(i,j); %29
                        dataPoints(index,26)= data.DB_Retrieving_Condition(i,j); %30
                        dataPoints(index,27)= data.Elevation(i,j); %31
                        dataPoints(index,28)= data.Error_CriticalReflectance_470(i,j); %32
                        dataPoints(index,29)= data.Error_CriticalReflectance_660(i,j); %33
                        dataPoints(index,30)= data.Error_PathRadiance_470(i,j); %34
                        dataPoints(index,31)= data.Error_PathRadiance_660(i,j); %35
                        dataPoints(index,32)= data.FittingError(i,j); %36
                        dataPoints(index,33)= data.Land_Type(i,j); %37
                        dataPoints(index,34)= data.MassConc(i,j); %38
                        dataPoints(index,35)= data.NumberPixelsUsedDB_412(i,j); %39
                        dataPoints(index,36)= data.NumberPixelsUsedDB_470(i,j); %40
                        dataPoints(index,37)= data.NumberPixelsUsedDB_660(i,j); %41
                        dataPoints(index,38)= data.OpticalDepthRatio(i,j); %42
                        dataPoints(index,39)= data.PathRadiance_470(i,j); %43
                        dataPoints(index,40)= data.PathRadiance_660(i,j); %44
                        dataPoints(index,41)= data.PixelsUsed_470(i,j); %45
                        dataPoints(index,42)= data.PixelsUsed_660(i,j); %46
                        dataPoints(index,43)= data.ProcPathPart_1(i,j); %47
                        dataPoints(index,44)= data.ProcPathPart_2(i,j); %48
                        dataPoints(index,45)= data.QA_Conf_470(i,j); %49
                        dataPoints(index,46)= data.QA_Conf_660(i,j); %50
                        dataPoints(index,47)= data.QualityWeight_CriticalReflectance_470(i,j); %51
                        dataPoints(index,48)= data.QualityWeight_CriticalReflectance_660(i,j); %52
                        dataPoints(index,49)= data.QualityWeight_PathRadiance_470(i,j); %53
                        dataPoints(index,50)= data.QualityWeight_PathRadiance_660(i,j); %54
                        %dataPoints(index,51)= data.Reflectance(i,j);
                        dataPoints(index,51) = years(y);



                        dataPoints(index,52)= data.ReflectanceDB_412(i,j); %55
                        dataPoints(index,53)= data.ReflectanceDB_470(i,j); %56
                        dataPoints(index,54)= data.ReflectanceDB_660(i,j); %57
                        dataPoints(index,55)= data.Reflectance_After_Cloud_Sceening_2100(i,j); %58
                        dataPoints(index,56)= data.Reflectance_After_Cloud_Screening_470(i,j); %59
                        dataPoints(index,57)= data.Reflectance_After_Cloud_Screening_660(i,j); %60
                        dataPoints(index,58)= data.SSA_DB_412(i,j); %61
                        dataPoints(index,59)= data.SSA_DB_470(i,j); %62
                        dataPoints(index,60)= data.SSA_DB_660(i,j); %63
                        dataPoints(index,61)= data.STD_AOD_DB_412(i,j); %64
                        dataPoints(index,62)= data.STD_AOD_DB_470(i,j); %65
                        dataPoints(index,63)= data.STD_AOD_DB_550(i,j); %66
                        dataPoints(index,64)= data.STD_AOD_DB_660(i,j); %67
                        %dataPoints(index,65)=
                        %data.STD_Reflectance(i,j) %68
                        dataPoints(index,66)= data.STD_Reflectance_After_Cloud_Screening_2100(i,j); %69
                        dataPoints(index,67)= data.STD_Reflectance_After_Cloud_Screening_470(i,j); %70
                        dataPoints(index,68)= data.STD_Reflectance_After_Cloud_Screening_660(i,j); %71
                        dataPoints(index,69)= data.ScatteringAngle(i,j); %72
                        dataPoints(index,70)= data.SensorAzimuth(i,j); %73
                        dataPoints(index,71)= data.SensorZenith(i,j); %74
                        dataPoints(index,72)= data.SnowCover(i,j); %75
                        dataPoints(index,73)= data.SnowIce(i,j); %76
                        dataPoints(index,74)= data.SolarAzimuth(i,j); %77
                        dataPoints(index,75)= data.SolarZenith(i,j); %78
                        dataPoints(index,76)= data.SurfaceReflectanceDB_412(i,j); %79
                        dataPoints(index,77)= data.SurfaceReflectanceDB_470(i,j); %80
                        dataPoints(index,78)= data.SurfaceReflectanceDB_660(i,j); %81
                        dataPoints(index,79)= data.SurfaceReflectance_2130(i,j); %82
                        dataPoints(index,80)= data.SurfaceReflectance_470(i,j); %83
                        dataPoints(index,81)= data.SurfaceReflectance_660(i,j); %84
                        dataPoints(index,82)= data.TotalOzone(i,j); %85
                        dataPoints(index,83)= data.TotalPrecWater(i,j); %86
                        dataPoints(index,84)= data.QA_Usefulness_470(i,j); %87
                        dataPoints(index,85)= data.QA_Usefulness_660(i,j); %88
                        dataPoints(index,86) = str2double(name{3});
                        dataPoints(index,87) = str2double(name{2}(5:7));
                        index=index+1;
                    end
                catch err

                end
            end
        end
            
        clear data
    end
    

end

