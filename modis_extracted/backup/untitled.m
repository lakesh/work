[b m n] = unique(data(:,4));
[row column] = size(b);
newdata = data;
clear data;

for i=1:row
    day = b(i,1);
    index = find(newdata(:,4) == day);
    data = newdata(index,:);
    save([num2str(day) '.mat'],'data');
end