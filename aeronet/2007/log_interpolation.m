[row column]=size(aeronetdata);
AOD_558 = zeros(row,1);
for i=1:row
    alpha = -log(aeronetdata(i,9)/aeronetdata(i,11)) / log(675/440);
    AOD_Value = aeronetdata(i,9) * (675/558)^alpha;
    AOD_558(i,1) = AOD_Value;
end
newAeronetData = [aeronetdata(:,1:10), AOD_558, aeronetdata(:,11:17) ];
aeronetdata=newAeronetData;
