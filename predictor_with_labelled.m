row=9;
column=13;
days = 365;


%%%%%%%%%%%% Initialize the model parameters %%%%%%%%%%%%%%%%%%%%%%%%
beta1 = 1.4544;
beta2 = 1.2633;
alpha1 = 1.7905;
alpha2 = 1.6462;
alpha3 = 0.001;


%alpha1 = 10;
%alpha2 = 1;
%alpha1 = 3.2511;
%alpha2 = 2.8858;

%Default AOD value for the dummy predictor
default_predicted = 0;


%Load the data
[data missingValuesIndex collocatedData] = loaddata_with_labelled();

collocatedData = collocatedData.collocatedData;

%Initialize the b matrix
% alpha*indicatorfunction*AODvalue
b = zeros(days*column*row,1);
for i=1:days*row*column
    b(i,1) = data(i,1) * data(i,2) * alpha1 + data(i,3) * data(i,4) * alpha2 + default_predicted * alpha3;
end


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


%We calculate the precision matrix for a single day only
%Then replicate the same matrix through out the year

%Save the row, column and the values of those elements who have non zero
%values which can be used later for sparse matrix creation
sparseIndexX = zeros(N,1);
sparseIndexY = zeros(N,1);
sparseIndexValue = zeros(N,1);

index = 1;
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
            
            
            % For the four corner elements number of neighbours = 2
            if i == 1 || i == N || (x == 1 && y == column) || (x == row && y == 1)
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
                sparseIndexX(index,1) = i;
                sparseIndexY(index,1) = j;
                sparseIndexValue(index,1) = -1;
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


QSpatial = sparse(sparseIndexX_WholeYear', sparseIndexY_WholeYear', sparseIndexValue_WholeYear', N*365, N*365);

%%%%%%%%%%%%%%%% End of calculating QSpatial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%% Start calculating QTemporal
%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sparseIndexX = zeros(N,1);
sparseIndexY = zeros(N,1);
sparseIndexValue = zeros(N,1);
index = 1;

for i=1:N*days
    %diagonal element
    sparseIndexX(index,1) = i;
    sparseIndexY(index,1) = i;
    sparseIndexValue(index,1) = 1;
    index = index + 1;
end


[m n] = size(sparseIndexX);

sparseIndexX_WholeYear = zeros(m,1);
sparseIndexY_WholeYear = zeros(m,1);
sparseIndexValue_WholeYear = zeros(m,1);

sparseIndexX_WholeYear(1:m,1)=sparseIndexX(1:m,1);
sparseIndexY_WholeYear(1:m,1)=sparseIndexY(1:m,1);
sparseIndexValue_WholeYear(1:m,1)=sparseIndexValue(1:m,1);


%Populate the remaining beta2 values
index = m+1;
j = N+1;
for i=1:N*364;
    sparseIndexX_WholeYear(index,1) = i;
    sparseIndexY_WholeYear(index,1) = j;
    sparseIndexValue_WholeYear(index,1) = -1;
    j=j+1;
    index = index + 1;
end

i=1;
for j=N+1:N*365
    sparseIndexX_WholeYear(index,1) = j;
    sparseIndexY_WholeYear(index,1) = i;
    sparseIndexValue_WholeYear(index,1) = -1;
    i=i+1;
    index = index + 1;
end

%Populate the remaining beta2 values in the diagonal elements
[m n] = size(sparseIndexX_WholeYear);
for i=1:m
    if (sparseIndexX_WholeYear(i,1) == sparseIndexY_WholeYear(i,1)) && sparseIndexX_WholeYear(i,1) > N && sparseIndexX_WholeYear(i,1) <= N*364
        sparseIndexValue_WholeYear(i,1) = sparseIndexValue_WholeYear(i,1) + 1;
    end
end

QTemporal = sparse(sparseIndexX_WholeYear', sparseIndexY_WholeYear', sparseIndexValue_WholeYear', N*365, N*365);

%%%%%%%%%%%%%%%% End calculating QTemporal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%% Calculate Q with alpha values %%%%%%%%%%%%%%%%%%%%%
sparseIndexX_Q1 = zeros(row,1);
sparseIndexY_Q1 = zeros(row,1);
sparseIndexValue_Q1 = zeros(row,1);

%For Q1 matrix
%Populate the diagonal values alpha1*delta1(xi) + alpha2*delta2(xi) +
%alpha3
%Indicator variable for this default predictor is always 1
index = 1;

for i=1:days*row*column
    sparseIndexX_Q1(index,1) = i;
    sparseIndexY_Q1(index,1) = i;
    sparseIndexValue_Q1(index,1) = alpha1 * data(i,2) + alpha2 * data(i,4) + alpha3;
    index = index+1;
end
Q1 = sparse(sparseIndexX_Q1', sparseIndexY_Q1', sparseIndexValue_Q1', days*row*column, days*row*column);

%%%%%%%%%%%%%%%%% End of calculating Q1 with alpha values%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%Get the index for the collocated and uncollocated points to sort the Q
%%%%%%%%matrix%%%%%%%%%%%%
[row column] = size(collocatedData);
collocatedIndex = zeros(row,1);
index = 1;
sizeX = 13;
sizeY = 9;
for i=1:row
    data = collocatedData(index,:);
    collocatedIndex(index,1) = (data(1,8) - 1)*sizeX*sizeY + (data(1,9)-1)*sizeX + data(1,10);
    index = index + 1;
end

unCollocatedIndex = linearMatrix(:,4);
unCollocatedIndex(collocatedIndex,:)=[];

%%%%%%%%%End of getting the index for the collocated points
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate QLL_Spatial and QUU_Spatial and QUL_Spatial
QLL_Spatial = beta1*QSpatial(collocatedIndex,:);
QLL_Spatial = QLL_Spatial(:,collocatedIndex);
QUU_Spatial = beta1*QSpatial(unCollocatedIndex,:);
QUU_Spatial = QUU_Spatial(:,unCollocatedIndex);
QUL_Spatial = beta1*QSpatial(unCollocatedIndex,:);
QUL_Spatial = QUL_Spatial(:,collocatedIndex);

%Calculate QLL_Temporal and QUU_Temporal and QUL Temporal
QLL_Temporal = beta2*QTemporal(collocatedIndex,:);
QLL_Temporal = QLL_Temporal(:,collocatedIndex);
QUU_Temporal = beta2*QTemporal(unCollocatedIndex,:);
QUU_Temporal = QUU_Temporal(:,unCollocatedIndex);
QUL_Temporal = beta2*QTemporal(unCollocatedIndex,:);
QUL_Temporal = QUL_Temporal(:,collocatedIndex);


%Calculate Q1LL and Q1UU and Q1UL
Q1LL = Q1(collocatedIndex,:);
Q1LL = Q1LL(:,collocatedIndex);
Q1UU = Q1(unCollocatedIndex,:);
Q1UU = Q1UU(:,unCollocatedIndex);
Q1UL = Q1(unCollocatedIndex,:);
Q1UL = Q1UL(:,collocatedIndex);

QUU = Q1UU + QUU_Temporal + QUU_Spatial;
QLL = Q1LL + QLL_Temporal + QLL_Spatial;
QUL = Q1UL + QUL_Temporal + QUL_Spatial;


bLL = b(collocatedIndex,:);
bUU = b(unCollocatedIndex,:);

uLL = QLL\bLL;
uUU = QUU\bUU;
u = uUU + QUU\(QUL*(collocatedData(:,5) - bLL));

%Saving the u values
save('u.mat','u');



