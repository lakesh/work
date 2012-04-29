rows = 42705;
collocated_misr_modis_aeronet = zeros(rows,3);
count = 1;

for i=1:rows
    if collocated_misr(i,1) == 0 && collocated_modis(i,1) == 0
        %Just take the label(if available) and the remaining are missing attributes
        collocated_misr_modis_aeronet(i,3) = collocated_modis(i,2);        
    elseif collocated_misr(i,1) == 0 && collocated_modis(i,1) ~= 0
        %AOD MODIS
        collocated_misr_modis_aeronet(i,2) = collocated_modis(i,1);
        %AOD AERONET
        collocated_misr_modis_aeronet(i,3) = collocated_modis(i,2);        
    elseif collocated_misr(i,1) ~= 0 && collocated_modis(i,1) == 0
        %AOD MISR
        collocated_misr_modis_aeronet(i,1) = collocated_misr(i,1);
        %AOD AERONET
        collocated_misr_modis_aeronet(i,3) = collocated_misr(i,2);        
    elseif collocated_misr(i,1) ~= 0 && collocated_modis(i,1) ~= 0
        %misr_hour = collocated_misr(i,3);
        %misr_minute = collocated_misr(i,4);
        %modis_hour = collocated_modis(i,3);
        %modis_minute = collocated_modis(i,4);
        %misr_total_time = misr_hour*60 + misr_minute;
        %modis_total_time = modis_hour*60 + modis_minute;
        %AOD MISR
        collocated_misr_modis_aeronet(i,1) = collocated_misr(i,1);
        %AOD MODIS
        collocated_misr_modis_aeronet(i,2) = collocated_modis(i,1);
        %AOD AERONET
        collocated_misr_modis_aeronet(i,3) = collocated_misr(i,2);        
    end
end