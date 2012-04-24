errorMISR = 0;
errorMODIS = 0;

missingData = find(collocatedData(:,1) == 0 | collocatedData(:,3) == 0);
collocatedData(missingData,:) = [];

[row column] = size(collocatedData);

for i=1:row
  
   [(collocatedData(i,1) - collocatedData(i,3))^2]
   [(collocatedData(i,3) - collocatedData(i,3))^2]
   errorMISR = errorMISR + (collocatedData(i,1) - collocatedData(i,3))^2;
   errorMODIS = errorMODIS + (collocatedData(i,2) - collocatedData(i,3))^2;
end

errorMISR = errorMISR/row;
errorMODIS = errorMODIS/row;