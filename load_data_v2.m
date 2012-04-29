function [inputData]= loaddata_v2()

    %For year 2006%
    interpolatedMISRPath = '/work/testset/MISR/';
    interpolatedMODISPath = '/work/testset/MODIS/';

    inputData = zeros(42705,12);

    %filePrefix = 'interpolatedData_';
    filePrefix = 'interpolatedData_';

    row = 9;
    column = 13;
    N = row*column;
    index = 1;


    %Load the MISR data
    for day=1:365
        try 
            load([interpolatedMISRPath filePrefix num2str(day) '.mat']);
            for i=1:row
                for j=1:column
                    %year
                    inputData(index,1) = data(i,j,1);
                    %dayofyear
                    inputData(index,2) = data(i,j,26);
                    %latitude
                    inputData(index,3) = data(i,j,6);
                    %longitude
                    inputData(index,4) = data(i,j,7);
                    
                    %AOD value
                    inputData(index,5) = data(i,j,10);
                    %Quality
                    if data(i,j,10) == 0 || isnan(data(i,j,10))
                        inputData(index,6) = 0;
                        %Replace the Nan values by 0 as well
                        inputData(index,5) = 0;
                    else
                        inputData(index,6) = 1;
                    end
                    inputData(index,7) = data(i,j,4);
                    inputData(index,8) = data(i,j,5);
                    index = index + 1;
                end
            end
        catch err
            %In case we don't have MISR data for a particular day use 0
            index = index + N;
        end

    end

    index = 1;

    %Load the MODIS data
    for day=1:365
        try 
            load([interpolatedMODISPath filePrefix num2str(day) '.mat']);
            for i=1:row
                for j=1:column
                    if data(i,j,23) < 0 || isnan(data(i,j,23))
                        inputData(index,9) = 0;
                    else
                        inputData(index,9) = data(i,j,23);
                    end
                    if inputData(index,9) == 0 || isnan(inputData(index,9)) || data(i,j,8) > 0 
                        %Replace the Nan values by 0 as well
                        inputData(index,10) = 0;
                    else
                        inputData(index,10) = 1;
                    end
                    inputData(index,11) = data(i,j,3);
                    inputData(index,12) = data(i,j,4);
                    index = index + 1;
                end
            end
        catch err
            %In case we don't have MODIS data for a particular day use 0
            index = index + N;
        end

    end
   
    
    clear N column data day err filePrefix i index interpolatedMISRPath interpolatedMODISPath j row;