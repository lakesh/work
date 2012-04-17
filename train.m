% Need to calculate derivate with respect to alpha1, alpha2, beta1 and
% beta2


%years = [2004 2005 2006 2007];
years = [2005 2006 2007];
dataPath =  '/work/';
filePrefix = 'collocated_MISR_MODIS_AERONET_';

%randomly initialize the alpha and beta values
alpha1 = 10;
alpha2 = 1;
alpha3 = 0.01; % for dummy predictor
beta1 = 0.5;
beta2 = 0.5;

%Default AOD value for the dummy predictor
default_predicted = 0;

row = 9;
column = 13;
numberOfDays = 365;
days = 365;

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

for j=1:N*days
    %diagonal element
    sparseIndexX(index,1) = j;
    sparseIndexY(index,1) = j;
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
k = N+1;
for j=1:N*364;
    sparseIndexX_WholeYear(index,1) = j;
    sparseIndexY_WholeYear(index,1) = k;
    sparseIndexValue_WholeYear(index,1) = -1;
    k=k+1;
    index = index + 1;
end

j=1;
for k=N+1:N*365
    sparseIndexX_WholeYear(index,1) = k;
    sparseIndexY_WholeYear(index,1) = j;
    sparseIndexValue_WholeYear(index,1) = -1;
    j=j+1;
    index = index + 1;
end

%Populate the remaining beta2 values in the diagonal elements
[m n] = size(sparseIndexX_WholeYear);
for j=1:m
    if (sparseIndexX_WholeYear(j,1) == sparseIndexY_WholeYear(j,1)) && sparseIndexX_WholeYear(j,1) > N && sparseIndexX_WholeYear(j,1) <= N*364
        sparseIndexValue_WholeYear(j,1) = sparseIndexValue_WholeYear(j,1) + 1;
    end
end

QTemporal = sparse(sparseIndexX_WholeYear', sparseIndexY_WholeYear', sparseIndexValue_WholeYear', N*365, N*365);

%%%%%%%%%%%%%%%% End calculating QTemporal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear sparseIndexX sparseIndexY sparseIndexValue sparseIndexValue_WholeYear sparseIndexX_WholeYear sparseIndexY_WholeYear




%Learning rate
learningRate = 0.001;

iterate = true;
numberOfIterations = 0;

while iterate == true
    numberOfIterations = numberOfIterations + 1;
    
    %Get the number of years
    [yearrow yearcolumn] = size(years);
    
    delta_log_alpha1 = 0;
    delta_log_alpha2 = 0;
    delta_log_beta1 = 0;
    delta_log_beta2 = 0;
    
    for i=1:yearcolumn
        year = years(i);
        load([dataPath filePrefix num2str(year) '.mat']);
        
        [row column] = size(collocatedData);
        
        %calculate the b matrix 
        b=zeros(row,1);
        for j=1:row
            b(j,1) = alpha1 * collocatedData(j,1) * collocatedData(j,2) + alpha2 * collocatedData(j,3) * collocatedData(j,4) + alpha3 * default_predicted;
        end
        
        %calculate the X matrix
        X = zeros(row,2);
        for j=1:row
            %Misr(xi) * IMisr(xi)
            X(j,1) = collocatedData(j,1) * collocatedData(j,2);
            %Modis(xi) * IModis(xi)
            X(j,2) = collocatedData(j,3) * collocatedData(j,4);
        end

        %calculate the E1 and E2 matrix
        E1 = zeros(row,row);
        E2 = zeros(row,row);

        for j=1:row
            E1(j,j) = collocatedData(j,2);
            E2(j,j) = collocatedData(j,4);
        end

        %calculate Sum(IModis(xi)*(yi^2-2*yi*Modis(xi))
        sumMODIS = 0;
        for j=1:row
            sumMODIS = sumMODIS + collocatedData(j,4) * (collocatedData(j,5) * collocatedData(j,5) - 2* collocatedData(j,5) * collocatedData(j,3));
        end


        %calculate Sum(Imisr(xi)*(yi^2-2*yi*Misr(xi))
        sumMISR = 0;
        for j=1:row
            sumMISR = sumMISR + collocatedData(j,2) * (collocatedData(j,5) * collocatedData(j,5) - 2* collocatedData(j,5) * collocatedData(j,1));
        end



        %Calculate the Q1 and Q2 matrix

        sparseIndexX_Q1 = zeros(row,1);
        sparseIndexY_Q1 = zeros(row,1);
        sparseIndexValue_Q1 = zeros(row,1);

        %For Q1 matrix
        %Populate the diagonal values alpha1*delta1(xi) + alpha2*delta2(xi) +
        %alpha3
        %Indicator variable for this default predictor is always 1
        index = 1;

        for j=1:row
            sparseIndexX_Q1(index,1) = j;
            sparseIndexY_Q1(index,1) = j;
            sparseIndexValue_Q1(index,1) = alpha1 * collocatedData(j,2) + alpha2 * collocatedData(j,4) + alpha3;
            index = index+1;
        end
        Q1 = sparse(sparseIndexX_Q1', sparseIndexY_Q1', sparseIndexValue_Q1', row, row);
        
        

        



        [row column] = size(collocatedData);

        %%%%Sort the QTemporal and QSpatial matrix based on collocated data
        %%%%%%%%%%%%%%%%

        %Get the index for the collocated points to sort the Q matrix
        collocatedIndex = zeros(row,1);
        index = 1;
        sizeX = 13;
        sizeY = 9;
        for j=1:row
            data = collocatedData(j,:);
            collocatedIndex(index,1) = (data(1,8) - 1)*sizeX*sizeY + (data(1,9)-1)*sizeX + data(1,10);
            index = index + 1;
        end

        %Calculate QLL_Spatial
        QLL_Spatial = QSpatial(collocatedIndex,:);
        QLL_Spatial = QLL_Spatial(:,collocatedIndex);

        %Calculate QLL_Temporal
        QLL_Temporal = QTemporal(collocatedIndex,:);
        QLL_Temporal = QLL_Temporal(:,collocatedIndex);

        %%%%%%%% End of sorting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Q2 = QLL_Spatial;
        Q3 = QLL_Temporal;

        %The Q matrix
        Q=Q1+beta1*Q2+beta2*Q3;


        %calculate the derivative
        delta_log_alpha1 = delta_log_alpha1 + (-sumMISR + 1/2 * trace(Q\E1) - X(:,1)'*(Q\b) + b'* (Q\(E1*(Q\b))) - b'*(Q\X(:,1)));
        delta_log_alpha2 = delta_log_alpha2 + (-sumMODIS + 1/2 * trace(Q\E2) - X(:,2)'*(Q\b) + b'* (Q\E2*(Q\b)) - b'*(Q\X(:,2)));
        delta_log_beta1 = delta_log_beta1 + (-beta1*collocatedData(:,5)'*Q2*collocatedData(:,5) + 1/2*trace(Q\Q2) - b'* (Q\(Q2*(Q\b))));
        delta_log_beta2 = delta_log_beta2 + (-beta2*collocatedData(:,5)'*Q3*collocatedData(:,5) + 1/2*trace(Q\Q3) - b'* (Q\(Q3*(Q\b))));

    end
    
    delta_log_alpha1 = alpha1*(delta_log_alpha1 - alpha1);
    delta_log_alpha2 = alpha2*(delta_log_alpha2 - alpha2);
    delta_log_beta1 = beta1*(delta_log_beta1 - beta1);
    delta_log_beta2 = beta2*(delta_log_beta2 - beta2);
    
    %Update the parameters
    log_alpha1 = log(alpha1) + learningRate * delta_log_alpha1;
    log_alpha2 = log(alpha2) + learningRate * delta_log_alpha2;
    log_beta1 = log(beta1) + learningRate * delta_log_beta1;
    log_beta2 = log(beta2) + learningRate * delta_log_beta2;

    alpha1_new = exp(log_alpha1);
    alpha2_new = exp(log_alpha2);
    beta1_new = exp(log_beta1);
    beta2_new = exp(log_beta2);

    difference = abs(alpha1_new-alpha1) + abs(alpha2_new-alpha2) + abs(beta1_new-beta1) + abs(beta2_new-beta2);

    alpha1 = alpha1_new;
    alpha2 = alpha2_new;
    beta1 = beta1_new;
    beta2 = beta2_new;

    if difference <= 0.001
        iterate = false;
    end
    %disp(sumMODIS);
end