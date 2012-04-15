count = 0;
for i=1:365
    try 
        load([num2str(i) '.mat']);
        a= find(aeronetdata(:,15) == 2);
        [row column] = size(a);
        count = count + row;
    catch
         
    end
end 
    