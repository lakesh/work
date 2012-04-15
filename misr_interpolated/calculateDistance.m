function [ distance ] = calculateDistance( latitude1,longitude1,latitude2,longitude2 )
    radius = 6371;
    dLat = degtorad(latitude2-latitude1); 
    %disp(dLat);
    dLon = degtorad(longitude2-longitude1); 
    %disp(dLon);
    a = sin(dLat/2) * sin(dLat/2) + cos(degtorad(latitude1)) * cos(degtorad(latitude2)) * sin(dLon/2) * sin(dLon/2); 
    %disp(a);
    c = 2 * atan2(sqrt(a), sqrt(1-a)); 
    %disp(c);
    distance = radius * c; 
end

