row=9;
column=13;
days = 365;

beta1 = 1.8357;
beta2 = 1.6097;
%alpha1 = 10;
%alpha2 = 1;
%alpha1 = 3.2511;
%alpha2 = 2.8858;
alpha1 = 1.9795;
alpha2 = 1.7886;
alpha3 = 0.001;

%Default AOD value for the dummy predictor
default_predicted = 0;


%Load the data
[data missingValuesIndex] = loaddata();

%Initialize the b matrix
% alpha*indicatorfunction*AODvalue
b = zeros(days*column*row,1);
for i=1:days*row*column
    b(i,1) = data(i,1) * data(i,2) * alpha2 + data(i,3) * data(i,4) * alpha1 + default_predicted * alpha3;
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


for i=1:N
    x = linearMatrix(i,1);
    y = linearMatrix(i,2);
    for j=1:N
        x1 = linearMatrix(j,1);
        y1 = linearMatrix(j,2);
        
        if(i == j)
            %diagonal element
            sparseIndexX(index,1) = i;
            sparseIndexY(index,1) = j;
            
            %Add the beta1 and partial beta2 values
            % For the four corner elements number of neighbours = 2
            if i == 1 || i == N || (x == 1 && y == column) || (x == row && y == 1)
                %sparseIndexValue(index,1) = 2*beta1 + alpha1 + alpha2 +  beta2;
                sparseIndexValue(index,1) = 2*beta1 + beta2;
            elseif x >=2 && y >= 2 && x < row && y < column
                % For the elements in the middle number of neighbours = 4
                %sparseIndexValue(index,1) = 4*beta1 + alpha1 + alpha2 + beta2;
                sparseIndexValue(index,1) = 4*beta1 + beta2;
            else
                % For the elements at the edge but not the corner
                % number of neighbours = 3
                %sparseIndexValue(index,1) = 3*beta1 + alpha1 + alpha2 + beta2;
                sparseIndexValue(index,1) = 3*beta1 + beta2;
            end
            index = index + 1;
        else
            %if the elements are adjacent 
            if (abs(x-x1) == 1 && abs(y-y1) == 0) || (abs(x-x1) == 0 && abs(y-y1) == 1)
                sparseIndexX(index,1) = i;
                sparseIndexY(index,1) = j;
                sparseIndexValue(index,1) = -1*beta1;
                index = index + 1;
            end
        end
    end
end


[m n] = size(sparseIndexX);

sparseIndexX_WholeYear = zeros(m*365,1);
sparseIndexY_WholeYear = zeros(m*365,1);
sparseIndexValue_WholeYear = zeros(m*365,1);

sparseIndexX_WholeYear(1:m,1)=sparseIndexX(1:m,1);
sparseIndexY_WholeYear(1:m,1)=sparseIndexY(1:m,1);
sparseIndexValue_WholeYear(1:m,1)=sparseIndexValue(1:m,1);

for day=1:days-1
    sparseIndexX_WholeYear(day*m+1:day*m+m,1)=sparseIndexX(1:m,1)+day*N;
    sparseIndexY_WholeYear(day*m+1:day*m+m,1)=sparseIndexY(1:m,1)+day*N;
    sparseIndexValue_WholeYear(day*m+1:day*m+m,1)=sparseIndexValue(1:m,1);
end


%Populate the remaining beta2 values

%Already have m*365 entries because of beta1. Add new entries because of beta2
index = m*365+1;
j = N+1;
for i=1:N*364;
    sparseIndexX_WholeYear(index,1) = i;
    sparseIndexY_WholeYear(index,1) = j;
    sparseIndexValue_WholeYear(index,1) = -beta2;
    j=j+1;
    index = index + 1;
end

i=1;
for j=N+1:N*365
    sparseIndexX_WholeYear(index,1) = j;
    sparseIndexY_WholeYear(index,1) = i;
    sparseIndexValue_WholeYear(index,1) = -beta2;
    i=i+1;
    index = index + 1;
end

%Populate the remaining beta2 values in the diagonal elements
[m n] = size(sparseIndexX_WholeYear);
for i=1:m
    if (sparseIndexX_WholeYear(i,1) == sparseIndexY_WholeYear(i,1)) && sparseIndexX_WholeYear(i,1) > N && sparseIndexX_WholeYear(i,1) <= N*364
        sparseIndexValue_WholeYear(i,1) = sparseIndexValue_WholeYear(i,1) + beta2;
    end
end

%Populating the alpha values

%Find the diagonal element indices first
[m n] = size(sparseIndexX_WholeYear);
index = 1;
diagonalMatrixIndex = zeros(row*column*days,1);
for i=1:m
    % Got the diagonal element
    if sparseIndexX_WholeYear(i,1) ==  sparseIndexY_WholeYear(i,1)
        diagonalMatrixIndex(index,1) = i;
        index = index + 1;
    end
end
%Now populate the alpha values
[m n] = size(diagonalMatrixIndex);
for i=1:m
    sparseIndexValue_WholeYear(diagonalMatrixIndex(i,1),1) = sparseIndexValue_WholeYear(diagonalMatrixIndex(i,1),1) + data(i,2) * alpha2 + data(i,4) * alpha1 + alpha3;
end


%TODO increase the efficiency of this routine
%Possible method save the xindex and the corresponding values to be
%changed. Then modify the Q matrix directly
%[m n] = size(missingValuesIndex);
%for i=1:m
%    missingValue = missingValuesIndex(i,1);
%    
    %Find the index of the data points which are connected with this data
    %point with missing values
%    a= find(sparseIndexY_WholeYear == missingValue);
%    [m1 n1] = size(a);
%    for j=1:m1
%        yIndex = a(j,1);
%        
%        %Find the corresponsing x index
%        xIndex = sparseIndexX_WholeYear(j,1);
%        value = sparseIndexValue_WholeYear(j,1);
%        
        %Since we will be removing the points with missing values,we need
        %to remove the beta smoothing as well.xIndex gives the diagonal
        %elements which have the data point with missing attribute as a
        %smoother. Just change the value of the diagonal element
%        diagonalElementIndex = find(sparseIndexX_WholeYear==xIndex&sparseIndexY_WholeYear==xIndex);
%        sparseIndexValue_WholeYear(diagonalElementIndex,1) = sparseIndexValue_WholeYear(diagonalElementIndex,1) - value;
%    end
%end

% u=sigma*b
% inv(sigma) = 2*Q
% sigma = 1/2*inv(Q)
% Qu = b/2
% u = Q\(b/2);

% Create the precision matrix
Q = sparse(sparseIndexX_WholeYear', sparseIndexY_WholeYear', sparseIndexValue_WholeYear', N*365, N*365);

clear day days i index j linearMatrix m n row ...
    sparseIndexValue_WholeYear sparseIndexX_WholeYear...
    sparseIndexY_WholeYear sparseIndexX sparseIndexY sparseIndexValue x y x1 y1 diagonalMatrixIndex

%Remove the values from the precision matrix having both modis and misr
%data missing. Remove the spatial and temporal smoothing caused by these
%values as well

%Q(missingValuesIndex,:)=[];
%Q(:,missingValuesIndex)=[];


%b(missingValuesIndex,:) = [];

u  = Q\b;


%Saving the u values and the missingValuesIndex
save('u.mat','u');
save('missingValuesIndex.mat','missingValuesIndex');


