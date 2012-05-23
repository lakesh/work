alpha1=1;
alpha2=1;
beta1=1;

global len learning_rate alpha3 row column days numberOfDays N linearMatrix index sparseIndexX sparseIndexY sparseIndexValue sparseIndexX_WholeYear sparseIndexY_WholeYear sparseIndexValue_WholeYear colloc_datas QSpatial labelled_indexes unlabelled_indexes
learning_rate = 0.0001;

alpha3 = 0.001;

row = 9;
column = 13;
days = 365;
numberOfDays = 365;

N = row*column;


% Matrix to hold the index for each data point 
linearMatrix = zeros(row*column*numberOfDays,4);
index = 1;

for day=1:numberOfDays
    for i=1:row
        for j=1:column
            linearMatrix(index,1) = i;
            linearMatrix(index,2) = j;
            linearMatrix(index,3) = day;
            linearMatrix(index,4) = index;
            index = index+1;
        end
    end
end
%%%%%%%%%%%%%%%%Start calculating QSpatial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We calculate the precision matrix for a single day only
%Then replicate the same matrix through out the year

%Save the row, column and the values of those elements who have non zero
%values which can be used later for sparse matrix creation
sparseIndexX = zeros(N,1);
sparseIndexY = zeros(N,1);
sparseIndexValue = zeros(N,1);

row = 9;
column = 13;

index = 1;
for j=1:N
    x = linearMatrix(j,1);
    y = linearMatrix(j,2);
    for k=1:N
        x1 = linearMatrix(k,1);
        y1 = linearMatrix(k,2);

        if(j == k)
            %diagonal element
            sparseIndexX(index,1) = j;
            sparseIndexY(index,1) = k;


            % For the four corner elements number of neighbours = 2
            if j == 1 || j == N || (x == 1 && y == column) || (x == row && y == 1)
                sparseIndexValue(index,1) = 2;
            elseif x >=2 && y >= 2 && x < row && y < column
                % For the elements in the middle number of neighbours = 4
                sparseIndexValue(index,1) = 4;
            else
                % For the elements at the edge but not the corner
                % number of neighbours = 3
                sparseIndexValue(index,1) = 3;
            end
            index = index + 1;
        else
            %if the elements are adjacent 
            if (abs(x-x1) == 1 && abs(y-y1) == 0) || (abs(x-x1) == 0 && abs(y-y1) == 1)
                sparseIndexX(index,1) = j;
                sparseIndexY(index,1) = k;
                sparseIndexValue(index,1) = -1;
                index = index + 1;
            end
        end
    end
end


[m n] = size(sparseIndexX);

sparseIndexX_WholeYear = zeros(m*days,1);
sparseIndexY_WholeYear = zeros(m*days,1);
sparseIndexValue_WholeYear = zeros(m*days,1);

sparseIndexX_WholeYear(1:m,1)=sparseIndexX(1:m,1);
sparseIndexY_WholeYear(1:m,1)=sparseIndexY(1:m,1);
sparseIndexValue_WholeYear(1:m,1)=sparseIndexValue(1:m,1);

for day=1:days-1
    sparseIndexX_WholeYear(day*m+1:day*m+m,1)=sparseIndexX(1:m,1)+day*N;
    sparseIndexY_WholeYear(day*m+1:day*m+m,1)=sparseIndexY(1:m,1)+day*N;
    sparseIndexValue_WholeYear(day*m+1:day*m+m,1)=sparseIndexValue(1:m,1);
end


QSpatial = sparse(sparseIndexX_WholeYear', sparseIndexY_WholeYear', sparseIndexValue_WholeYear', N*days, N*days);

clear sparseIndexValue sparseIndexValue_WholeYear sparseIndexX sparseIndexX_WholeYear sparseIndexY sparseIndexY_WholeYear
%%%%%%%%%%%%%%%% End of calculating QSpatial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%years=[2004 2005 2006 2007];
%years=[2003 2004 2005 2007];
%years=[2003 2004 2005 2006];
years=[2004 2005 2006 2007];

total_years = 1;
%total_years = 1;

colloc_datas = cell(total_years,1);
%colloc_datas{1} = load('/work/collocated_MISR_MODIS_AERONET_NEW_2005.mat');
%colloc_datas{2} = load('/work/collocated_MISR_MODIS_AERONET_NEW_2006.mat');
%colloc_datas{3} = load('/work/collocated_MISR_MODIS_AERONET_NEW_2007.mat');
%colloc_datas{4} = load('/work/collocated_MISR_MODIS_AERONET_NEW_2004.mat');
colloc_datas{1}.collocated_misr_modis_aeronet = colloc_data;
%colloc_datas{1} = load('/work/collocated_MISR_MODIS_AERONET_NEW_2004.mat');
%colloc_datas{2} = load('/work/collocated_MISR_MODIS_AERONET_NEW_2005.mat');
%colloc_datas{3} = load('/work/collocated_MISR_MODIS_AERONET_NEW_2006.mat');
%colloc_datas{4} = load('/work/collocated_MISR_MODIS_AERONET_NEW_2007.mat');
%colloc_datas{1}.collocated_misr_modis_aeronet = collocated_misr_modis_aeronet;
%colloc_datas{1}.collocated_misr_modis_aeronet = colloc_data;
labelled_indexes = cell(total_years,1);
unlabelled_indexes = cell(total_years,1);



for i=1:total_years
    data = colloc_datas{i}.collocated_misr_modis_aeronet;
    labelled_indexes{i}.labelled_index = find(data(:,3) ~= 0);
    unlabelled_indexes{i}.unlabelled_index = linearMatrix(:,4);
    unlabelled_indexes{i}.unlabelled_index(labelled_indexes{i}.labelled_index,:)=[];
end
    
len = days*N;
