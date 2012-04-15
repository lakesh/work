[b m n] = unique(data(:,4));
newData = data;
[row column] = size(b);

for i=1:row
    index = find(newData(:,4) == b(i,1));
    data = newData(index,:);
    mkdir(num2str(b(i,1)));
    save([num2str(b(i,1)) '/' num2str(b(i,1)) '.mat'],'data');
end