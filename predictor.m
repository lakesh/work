row=9;
column=13;
days = 365;

alpha1=124.0082;
alpha2=101.4777;
%alpha1=1;
%alpha2=0;

beta1=29.6218;
beta2=0;
alpha3=0.001;


%Default AOD value for the dummy predictor
default_predicted = 0;


%Load the data
load('collocated_MISR_MODIS_AERONET_NEW_2004.mat');
data=collocated_misr_modis_aeronet;

%Initialize the b matrix
% alpha*indicatorfunction*AODvalue
b = zeros(days*column*row,1);
for i=1:days*row*column
    b(i,1) = 2*(data(i,1)* alpha1 + data(i,2) * alpha2);
    %b(i,1) = 2*(data(i,1) * alpha2 + data(i,3) * alpha1);
end

%Replace the NaN values by 0
%b(isnan(b(:,1)),:) = 0;

% Matrix to hold the index for each data point 
linearMatrix = zeros(row*column*days,4);
index = 1;

for day=1:days
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

N = row*column;

index = 1;

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
Q2=beta1*QSpatial;
clear sparseIndexValue sparseIndexValue_WholeYear sparseIndexX sparseIndexX_WholeYear sparseIndexY sparseIndexY_WholeYear



truth_aeronet = data(:, 3);
forecast_modis = data(:, 2);
forecast_misr = data(:, 1);

indicator_modis = (forecast_modis ~= 0);
indicator_misr = (forecast_misr ~= 0);
indicator_dummy_whole = (forecast_misr == forecast_misr);

len=days*row*column;
Q1 = alpha1 * spdiags(indicator_misr, 0, len, len) + alpha2 * spdiags(indicator_modis, 0, len, len) + alpha3 * spdiags(indicator_dummy_whole, 0, len, len);


Q=2*(Q1+Q2);


%Save the row, column and the values of those elements who have non zero
%values which can be used later for sparse matrix creation

u = Q\b;


%Saving the u values and the missingValuesIndex
save('u.mat','u');



