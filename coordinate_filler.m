

files=dir('./*.mat');
for counter=1:length(files)
    AODFileName = files(counter).name;
    load(AODFileName);
    [row column] = size(aeronetdata);
    newaeronetdata = zeros(row,17);
    newaeronetdata(:,1:15)=aeronetdata(:,:);
    for i=1:row
        if newaeronetdata(i,14) == 1
            newaeronetdata(i,16) = 42.527788;
            newaeronetdata(i,17) = -71.26889;
        elseif newaeronetdata(i,14) == 2    
            newaeronetdata(i,16) = 42.53200;
            newaeronetdata(i,17) = -72.1880;
        elseif newaeronetdata(i,14) == 3
            newaeronetdata(i,16) = 41.30;
            newaeronetdata(i,17) = -70.550;
        elseif newaeronetdata(i,14) == 4
            newaeronetdata(i,16) = 43.10889;
            newaeronetdata(i,17) = -70.947778;
        end
    end
    aeronetdata = newaeronetdata;
    save(AODFileName,'aeronetdata');
end