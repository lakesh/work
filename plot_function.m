function_array=zeros(150,150);

alpha1=1:150;
alpha2=1:150;
[X Y] = meshgrid(alpha1,alpha2);
        
maxValue = 0;
maxValueAlpha1 = 0;
maxValueAlpha2 = 0;

%Initialize the alpha parameters
alpha1=1;
alpha2=2;
alpha3 = 0.001;


iteration = 1;
iter_alpha1=1;
iter_alpha2=1;

%Total number of years
total_years = 1;


% Grid size and number of days
row = 9;
column = 13;
N=row*column;
numberOfDays = 365;

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

index = 1;
linearMatrix = zeros(500,4);
for i=1:500
    linearMatrix(index,4) = index;
    index = index+1;
end
%Mark the labelled and unlabelled points
labelled_indexes = cell(total_years,1);
unlabelled_indexes = cell(total_years,1);

for i=1:total_years
    data = colloc_data;
    labelled_indexes{i}.labelled_index = find(data(:,3) ~= 0);
    unlabelled_indexes{i}.unlabelled_index = linearMatrix(:,4);
    unlabelled_indexes{i}.unlabelled_index(labelled_indexes{i}.labelled_index,:)=[];
end

for alpha1=1:150
    tic;
    for alpha2=1:150
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

        sigma_inv = 2*Q1;
        
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

        F = -0.5*(truth_aeronet(labelled_indexes{1}.labelled_index,:)-u)'*sigma_star_inv*(truth_aeronet(labelled_indexes{1}.labelled_index,:)-u) + 0.5*(2 * sum(log(diag(chol(sigma_star_inv)))));
        
        function_array(iter_alpha1,iter_alpha2) = F;
        iter_alpha2 = iter_alpha2+1;
        
        if F > maxValue 
            maxValue = F;
            maxValueAlpha1 = iter_alpha1;
            maxValueAlpha2 = iter_alpha2;
        end
    end
    iter_alpha2=1;
    iter_alpha1 = iter_alpha1+1;
    toc;
end

