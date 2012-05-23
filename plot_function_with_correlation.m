function_array=zeros(50,50,50);

alpha1=1:3:150;
alpha2=1:3:150;
beta1=1:3:150;
[X Y Z] = meshgrid(alpha1,alpha2,beta1);
        
maxValue = -10000;
maxValueAlpha1 = 0;
maxValueAlpha2 = 0;
maxValueBeta1 = 0;

%Initialize the alpha parameters
alpha1=1;
alpha2=2;
alpha3 = 0.001;
beta1 = 1;


iteration = 1;
iter_alpha1=1;
iter_alpha2=1;
iter_beta1=1;

%Total number of years
total_years = 1;


% Grid size and number of days
row = 9;
column = 13;
N=row*column;
numberOfDays = 365;
days=365;



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

%index = 1;
%linearMatrix = zeros(500,4);
%for i=1:500
%    linearMatrix(index,4) = index;
%    index = index+1;
%end

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

clear sparseIndexValue sparseIndexValue_WholeYear sparseIndexX sparseIndexX_WholeYear sparseIndexY sparseIndexY_WholeYear


%Mark the labelled and unlabelled points
labelled_indexes = cell(total_years,1);
unlabelled_indexes = cell(total_years,1);

for i=1:total_years
    data = colloc_data;
    labelled_indexes{i}.labelled_index = find(data(:,3) ~= 0);
    unlabelled_indexes{i}.unlabelled_index = linearMatrix(:,4);
    unlabelled_indexes{i}.unlabelled_index(labelled_indexes{i}.labelled_index,:)=[];
end

for alpha1=1:3:150
    tic;
    for alpha2=1:3:150
        for beta1=1:3:150
            truth_aeronet = colloc_data(:,3);
            forecast_misr = colloc_data(:,1);
            forecast_modis = colloc_data(:,2);

            indicator_modis = (forecast_modis ~= 0);
            indicator_misr = (forecast_misr ~= 0);
            indicator_dummy = (forecast_misr == forecast_misr);

            [row column] = size(colloc_data);
            len = row;

            b = 2 * (alpha1 * spdiags(indicator_misr, 0, len, len) * forecast_misr + alpha2 * spdiags(indicator_modis, 0, len, len)* forecast_modis );

            Q1 = alpha1*spdiags(indicator_misr, 0, len, len) + alpha2*spdiags(indicator_modis, 0, len, len) + alpha3 *spdiags(indicator_dummy, 0, len, len);
            Q2 = beta1*QSpatial;

            sigma_inv = 2*(Q1+Q2);

            Q_LL = sigma_inv(labelled_indexes{i}.labelled_index,:);
            Q_LL = Q_LL(:,labelled_indexes{i}.labelled_index);
            Q_UU = sigma_inv(unlabelled_indexes{i}.unlabelled_index,:);
            Q_UU = Q_UU(:,unlabelled_indexes{i}.unlabelled_index);
            Q_LU = sigma_inv(labelled_indexes{i}.labelled_index,:);
            Q_LU = Q_LU(:, unlabelled_indexes{i}.unlabelled_index);
            Q_UL = sigma_inv(unlabelled_indexes{i}.unlabelled_index,:);
            Q_UL = Q_UL(:, labelled_indexes{i}.labelled_index);

            sigma_star_inv = Q_LL - (Q_LU / Q_UU) * Q_UL;

            u = sigma_star_inv\(b(labelled_indexes{1}.labelled_index,:));

            F = -0.5*(truth_aeronet-u)'*sigma_star_inv*(truth_aeronet-u) + 0.5*(2 * sum(log(diag(chol(sigma_star_inv)))));

            function_array(iter_alpha1,iter_alpha2, iter_beta1) = F;
            iter_beta1 = iter_beta1+1;

            if F > maxValue 
                maxValue = F;
                maxValueAlpha1 = iter_alpha1;
                maxValueAlpha2 = iter_alpha2;
                maxValueBeta1 = iter_beta1;
            end
        end
        iter_alpha2=iter_alpha2+1;
        iter_beta1=1;
    end
    iter_alpha2=1;
    iter_beta1=1;
    iter_alpha1 = iter_alpha1+1;
    toc;
end

