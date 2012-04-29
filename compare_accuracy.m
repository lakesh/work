%load('collocated_MODIS_AERONET');
%load('collocated_MODIS_AERONET_2004');
load('collocated_MISR_MODIS_AERONET_2005');
load('u');
%load('missingValuesIndex');

row = 9;
column = 13;

[m n] = size(collocatedData);
%[m1 n1] = size(missingValuesIndex);

resultData = zeros(m,3);

for i=1:m
    resultData(i,1) = collocatedData(i,1);
    resultData(i,2) = collocatedData(i,3);
    
    %index = collocatedData(i,6)-1;
    index = collocatedData(i,8)-1;
    %since we are removing the data points with missing attributes, get the
    %new index of the elements in the new array
    index = index*row*column + (collocatedData(i,9) - 1) * column + collocatedData(i,10);
    resultData(i,3) = u(index,1);
    resultData(i,4) = collocatedData(i,5);
end

errorMODIS = 0;
errorMISR = 0;
errorCRF = 0;

for i=1:m
    errorMISR = errorMISR + (resultData(i,1) - resultData(i,4))^2;
    errorMODIS = errorMODIS + (resultData(i,2) - resultData(i,4))^2;
    errorCRF = errorCRF + (resultData(i,3) - resultData(i,4))^2;
end

errorMISR = errorMISR/m;
errorMODIS = errorMODIS/m;
errorCRF = errorCRF/m;