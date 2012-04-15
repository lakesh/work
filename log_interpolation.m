files=dir('./*.mat');
for counter=1:length(files)
    AODFileName = files(counter).name;
    load(AODFileName);
    [row column]=size(aeronetdata);
    AOD_558 = zeros(row,1);
    for i=1:row
        alpha = -log(aeronetdata(i,9)/aeronetdata(i,11)) / log(675/440);
        AOD_Value = aeronetdata(i,9) * (675/558)^alpha;
        AOD_558(i,1) = AOD_Value;
    end
    newAeronetData = [aeronetdata(:,1:9), AOD_558, aeronetdata(:,10:17) ];
    aeronetdata=newAeronetData;
    save(AODFileName,'aeronetdata');
    clear newAeronetData row column alpha AOD_Value AOD_558
end
