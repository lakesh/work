function [inputData missingValuesIndex] = loaddata()

    %For year 2006%
    interpolatedMISRPath = '/work/testset/MISR/';
    interpolatedMODISPath = '/work/testset/MODIS/';

    inputData = zeros(42705,4);

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
                    inputData(index,1) = data(i,j,8);
                    if data(i,j,8) == 0 || isnan(data(i,j,8))
                        inputData(index,2) = 0;
                        %Replace the Nan values by 0 as well
                        inputData(index,1) = 0;
                    else
                        inputData(index,2) = 1;
                    end
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
                    if data(i,j,23) < 0
                        inputData(index,3) = 0;
                    else
                        inputData(index,3) = data(i,j,23);
                    end
                    if inputData(index,3) == 0 || isnan(inputData(index,3)) || data(i,j,8) > 0 
                        %Replace the Nan values by 0 as well
                        inputData(index,3) = 0;
                        inputData(index,4) = 0;
                    else
                        inputData(index,4) = 1;
                    end
                    index = index + 1;
                end
            end
        catch err
            %In case we don't have MODIS data for a particular day use 0
            index = index + N;
        end

    end
    
    
    zeroMISR =  find(inputData(:,1)==0);
    zeroMODIS = find(inputData(:,3)==0);
    missingValuesIndex = intersect(zeroMISR,zeroMODIS);
    
    clear N column data day err filePrefix i index interpolatedMISRPath interpolatedMODISPath j row;