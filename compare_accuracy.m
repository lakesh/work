load('collocated_MISR_MODIS_AERONET_NEW_2004.mat');
load('u.mat');

colloc_data=collocated_misr_modis_aeronet;

errorMISR=0;
errorMODIS=0;
errorCRF=0;

a=find(colloc_data(:,1) ~= 0 & colloc_data(:,3) ~= 0);
b=find(colloc_data(:,2) ~= 0 & colloc_data(:,3) ~= 0);

colloc_misr = colloc_data(a,:);
colloc_modis = colloc_data(b,:);

u_misr=u(a,:);
u_modis=u(b,:);

[row column] = size(colloc_misr);
for i=1:row
    errorMISR = errorMISR + (colloc_misr(i,1) - colloc_misr(i,3)) * (colloc_misr(i,1) - colloc_misr(i,3));
    errorCRF = errorCRF + (u_misr(i,1) - colloc_misr(i,3)) * (u_misr(i,1) - colloc_misr(i,3));
end

errorCRF = 0;

[row column] = size(colloc_modis);
for i=1:row
    errorMODIS = errorMODIS + (colloc_modis(i,2) - colloc_modis(i,3)) * (colloc_modis(i,2) - colloc_modis(i,3));
    errorCRF = errorCRF + (u_modis(i,1) - colloc_modis(i,3)) * (u_modis(i,1) - colloc_modis(i,3));
end

