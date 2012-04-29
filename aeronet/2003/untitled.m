count = 0;
for i =1:365
    try 
        load([num2str(i) '.mat']);
        for j=1:4
            a=find(aeronetdata(:,15) == j);
            if isempty(a)

            else
                count = count + 1;
            end
        end
    catch
    end
end