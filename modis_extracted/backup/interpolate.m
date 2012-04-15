leftBoundary = -72.1889;
rightBoundary = -70.550;
topBoundary = 43.10889;
bottomBoundary = 41.3;

deltaX = 0.1349;
deltaY = 0.179562;

missingValue = -9999;

[xq yq] = meshgrid(41.3 + deltaX/2:0.1349:43.10889, -72.1889+deltaY/2:0.179562:-70.550);
[row column]=size(xq);

interpolatedData = zeros(row,column,26);
%dataWithoutMissingAttributes = data;
dataWithNormalization = data;

index = find(dataWithNormalization == -9999);
dataWithNormalization(index) = 1.0;

interpolatedData(:,:,1) = data(1,1);
interpolatedData(:,:,2) = data(1,2);
interpolatedData(:,:,3) = data(1,3);
interpolatedData(:,:,4) = griddata(data(:,6),data(:,7),data(:,4),xq,yq,'nearest');
interpolatedData(:,:,5) = griddata(data(:,6),data(:,7),data(:,5),xq,yq,'nearest');
interpolatedData(:,:,6) = xq;
interpolatedData(:,:,7) = yq;

dataWithoutMissingAttributes = data;
index = find(dataWithoutMissingAttributes(:,8) == -9999);
dataWithoutMissingAttributes(index,:) = [];
interpolatedData(:,:,8) = griddata(dataWithoutMissingAttributes(:,6),dataWithoutMissingAttributes(:,7),dataWithoutMissingAttributes(:,8),xq,yq);
interpolatedData(:,:,9) = griddata(dataWithNormalization(:,6),dataWithNormalization(:,7),dataWithNormalization(:,8),xq,yq);

dataWithoutMissingAttributes = data;
index = find(dataWithoutMissingAttributes(:,9) == -9999);
dataWithoutMissingAttributes(index,:) = [];
interpolatedData(:,:,10) = griddata(dataWithoutMissingAttributes(:,6),dataWithoutMissingAttributes(:,7),dataWithoutMissingAttributes(:,9),xq,yq);
interpolatedData(:,:,11) = griddata(dataWithNormalization(:,6),dataWithNormalization(:,7),dataWithNormalization(:,9),xq,yq);


dataWithoutMissingAttributes = data;
index = find(dataWithoutMissingAttributes(:,10) == -9999);
dataWithoutMissingAttributes(index,:) = [];
interpolatedData(:,:,12) = griddata(dataWithoutMissingAttributes(:,6),dataWithoutMissingAttributes(:,7),dataWithoutMissingAttributes(:,10),xq,yq);
interpolatedData(:,:,13) = griddata(dataWithNormalization(:,6),dataWithNormalization(:,7),dataWithNormalization(:,10),xq,yq);

dataWithoutMissingAttributes = data;
index = find(dataWithoutMissingAttributes(:,11) == -9999);
dataWithoutMissingAttributes(index,:) = [];
interpolatedData(:,:,14) = griddata(dataWithoutMissingAttributes(:,6),dataWithoutMissingAttributes(:,7),dataWithoutMissingAttributes(:,11),xq,yq);
interpolatedData(:,:,15) = griddata(dataWithNormalization(:,6),dataWithNormalization(:,7),dataWithNormalization(:,11),xq,yq);

dataWithoutMissingAttributes = data;
index = find(dataWithoutMissingAttributes(:,12) == -9999);
dataWithoutMissingAttributes(index,:) = [];
interpolatedData(:,:,16) = griddata(dataWithoutMissingAttributes(:,6),dataWithoutMissingAttributes(:,7),dataWithoutMissingAttributes(:,12),xq,yq);
interpolatedData(:,:,17) = griddata(dataWithNormalization(:,6),dataWithNormalization(:,7),dataWithNormalization(:,12),xq,yq);



dataWithoutMissingAttributes = data;
index = find(dataWithoutMissingAttributes(:,13) == -9999);
dataWithoutMissingAttributes(index,:) = [];
interpolatedData(:,:,18) = griddata(dataWithoutMissingAttributes(:,6),dataWithoutMissingAttributes(:,7),dataWithoutMissingAttributes(:,13),xq,yq);
interpolatedData(:,:,19) = griddata(dataWithNormalization(:,6),dataWithNormalization(:,7),dataWithNormalization(:,13),xq,yq);


dataWithoutMissingAttributes = data;
index = find(dataWithoutMissingAttributes(:,14) == -9999);
dataWithoutMissingAttributes(index,:) = [];
interpolatedData(:,:,20) = griddata(dataWithoutMissingAttributes(:,6),dataWithoutMissingAttributes(:,7),dataWithoutMissingAttributes(:,14),xq,yq);
interpolatedData(:,:,21) = griddata(dataWithNormalization(:,6),dataWithNormalization(:,7),dataWithNormalization(:,14),xq,yq);


dataWithoutMissingAttributes = data;
index = find(dataWithoutMissingAttributes(:,15) == -9999);
dataWithoutMissingAttributes(index,:) = [];
interpolatedData(:,:,22) = griddata(dataWithoutMissingAttributes(:,6),dataWithoutMissingAttributes(:,7),dataWithoutMissingAttributes(:,15),xq,yq);
interpolatedData(:,:,23) = griddata(dataWithNormalization(:,6),dataWithNormalization(:,7),dataWithNormalization(:,15),xq,yq);

dataWithoutMissingAttributes = data;
index = find(dataWithoutMissingAttributes(:,16) == -9999);
dataWithoutMissingAttributes(index,:) = [];
interpolatedData(:,:,24) = griddata(dataWithoutMissingAttributes(:,6),dataWithoutMissingAttributes(:,7),dataWithoutMissingAttributes(:,16),xq,yq);
interpolatedData(:,:,25) = griddata(dataWithNormalization(:,6),dataWithNormalization(:,7),dataWithNormalization(:,16),xq,yq);


interpolatedData(:,:,26) = data(1,17);
z1 = griddata(data(:,6),data(:,7),data(:,8),xq,yq);