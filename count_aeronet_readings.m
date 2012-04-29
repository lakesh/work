count = 0;

for day=1:365
    try 
        load([num2str(day) '.mat']);
        for site=1:4
            a=find(aeronetdata(:,15) == site);
            if isempty(a)
            else
                count = count + 1;
            end
        end
    catch
    end
end
        