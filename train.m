% Need to calculate derivate with respect to alpha1, alpha2, beta1 and
% beta2

%load the MISR,MODIS, AERONET collocated data
load('/work/collocated_MISR_MODIS_AERONET.mat');

[row column] = size(collocatedData);


%randomly initialize the alpha and beta values
alpha1 = 10;
alpha2 = 1;
alpha3 = 0.01; % for dummy predictor
beta1 = 0.5;
beta2 = 0.5;

%Default AOD value for the dummy predictor
default_predicted = 0;

%calculate the X matrix
X = zeros(row,2);
for i=1:row
    %Misr(xi) * IMisr(xi)
    X(i,1) = collocatedData(i,1) * collocatedData(i,2);
    %Modis(xi) * IModis(xi)
    X(i,2) = collocatedData(i,3) * collocatedData(i,4);
end

%calculate the E1 and E2 matrix
E1 = zeros(row,row);
E2 = zeros(row,row);

for i=1:row
    E1(i,i) = collocatedData(i,2);
    E2(i,i) = collocatedData(i,4);
end

%calculate Sum(IModis(xi)*(yi^2-2*yi*Modis(xi))
sumMODIS = 0;
for i=1:row
    sumMODIS = sumMODIS + collocatedData(i,4) * (collocatedData(i,5) * collocatedData(i,5) - 2* collocatedData(i,5) * collocatedData(i,3));
end


%calculate Sum(Imisr(xi)*(yi^2-2*yi*Misr(xi))
sumMISR = 0;
for i=1:row
    sumMISR = sumMISR + collocatedData(i,2) * (collocatedData(i,5) * collocatedData(i,5) - 2* collocatedData(i,5) * collocatedData(i,1));
end


%traceLog = zeros(10000000);
%traceIndex = 1;

iterate = true;
numberOfIterations = 0;
while iterate == true
    numberOfIterations = numberOfIterations + 1;
    %calculate the b matrix 
    b=zeros(row,1);
    for i=1:row
        b(i,1) = alpha1 * collocatedData(i,1) * collocatedData(i,2) + alpha2 * collocatedData(i,3) * collocatedData(i,4) + alpha3 * default_predicted;
    end


    %Calculate the Q1 and Q2 matrix

    %TODO calculate a third matrix Q3 for temporal correlation

    sparseIndexX_Q1 = zeros(row,1);
    sparseIndexY_Q1 = zeros(row,1);
    sparseIndexValue_Q1 = zeros(row,1);

    sparseIndexX_Q2 = zeros(1,1);
    sparseIndexY_Q2 = zeros(1,1);
    sparseIndexValue_Q2 = zeros(1,1);

    index = 1;

    %For spatial correlation (Q2)
    %Populate the -beta1 values
    for i=1:row
        for j=1:row
            if i ~= j && collocatedData(i,8) == collocatedData(j,8) && (abs(collocatedData(i,9)-collocatedData(j,9)) == 1 && abs(collocatedData(i,10)-collocatedData(j,10)) == 0) && (abs(collocatedData(i,10)-collocatedData(j,10)) == 1 && abs(collocatedData(i,9)-collocatedData(j,9)) == 0)
                sparseIndexX_Q2(index,1) = i;
                sparseIndexY_Q2(index,1) = j;
                sparseIndexValue_Q2(index,1) = -1;
                index = index+1;
            end
        end
    end
    Q2 = [];

    %check if there is any correlation
    if all(sparseIndexX_Q2) == 1
        %If yes initialize the Q2 matrix with the correlated data
        Q2= sparse(sparseIndexX_Q2', sparseIndexY_Q2', sparseIndexValue_Q2', row, row);
        %In Q2 matrix the absolute sum of the non diagonal elements == the value of the
        %diagonal element so initialize the diagonal element with that value
        for i=1:row
            position = find(sparseIndexX_Q2(:,1) == i);
            [row1 column1] = size(position);
            sparseIndexX_Q2(index,1) = i;
            sparseIndexY_Q2(index,1) = i;
            sparseIndexValue_Q2(index,1) = row1;
            index = index+1;
        end
    else
        Q2 = sparse(row, row);
    end



    %For Q1 matrix
    %Populate the diagonal values alpha1*delta1(xi) + alpha2*delta2(xi) +
    %alpha3
    %Indicator variable for this default predictor is always 1
    index = 1;

    for i=1:row
        sparseIndexX_Q1(index,1) = i;
        sparseIndexY_Q1(index,1) = i;
        sparseIndexValue_Q1(index,1) = alpha1 * collocatedData(i,2) + alpha2 * collocatedData(i,4) + alpha3;
        index = index+1;
    end
    Q1 = sparse(sparseIndexX_Q1', sparseIndexY_Q1', sparseIndexValue_Q1', row, row);

    %The Q matrix
    Q=Q1+Q2;




    %Learning rate
    n = 0.01;

    %Q_full = full(Q);
    %calculate the derivative
    delta_log_alpha1 = alpha1* (-sumMISR + 1/2 * trace(Q\E1) - X(:,1)'*(Q\b) + b'* (Q\(E1*(Q\b))) - b'*(Q\X(:,1)) - (alpha1));
    delta_log_alpha2 = alpha2* (-sumMODIS + 1/2 * trace(Q\E2) - X(:,2)'*(Q\b) + b'* (Q\E2*(Q\b)) - b'*(Q\X(:,2)) - (alpha2));
    delta_log_beta2 = beta2 * (beta2*collocatedData(:,5)'*Q2*collocatedData(:,5) + 1/2*trace(Q\Q2) - b'* (Q\(Q2*(Q\b))));

    %Update the parameters
    log_alpha1 = log(alpha1) + n * delta_log_alpha1;
    log_alpha2 = log(alpha2) + n * delta_log_alpha2;
    log_beta2 = log(beta2) + n * delta_log_beta2;

    alpha1_new = exp(log_alpha1);
    alpha2_new = exp(log_alpha2);
    beta2_new = exp(log_beta2);
    
    difference = abs(alpha1_new-alpha1) + abs(alpha2_new-alpha2) + abs(beta2_new-beta2);
    
    alpha1 = alpha1_new;
    alpha2 = alpha2_new;
    beta2 = beta2_new;
    
    if difference <= 0.0001
        iterate = false;
    end
    %traceLog(traceIndex,1) = difference;
    %traceIndex = traceIndex + 1;
    disp(difference);
end