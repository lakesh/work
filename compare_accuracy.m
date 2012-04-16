%load('collocated_MODIS_AERONET');
load('collocated_MODIS_AERONET_2006');
load('u');
%load('missingValuesIndex');

row = 9;
column = 13;

[m n] = size(collocatedData);
%[m1 n1] = size(missingValuesIndex);

resultData = zeros(m,3);

for i=1:m
    resultData(i,1) = collocatedData(i,1);
    resultData(i,2) = collocatedData(i,2);
    index = collocatedData(i,6)-1;
    %since we are removing the data points with missing attributes, get the
    %new index of the elements in the new array
    index = index*row*column + (collocatedData(i,9) - 1) * column + collocatedData(i,10);
    %a = find(missingValuesIndex < index);
    %[m2 n2] = size(a);
    %index = index - m2;
    resultData(i,3) = u(index,1);
end

errorMODIS = 0;
errorCRF = 0;

for i=1:m
    errorMODIS = errorMODIS + (resultData(i,1) - resultData(i,2))^2;
    errorCRF = errorCRF + (resultData(i,2) - resultData(i,3))^2;
end

errorMODIS = errorMODIS/m;
errorCRF = errorCRF/m;