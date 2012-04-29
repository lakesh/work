function [alpha1, alpha2] = train_CRF_NIPS_v2_lakesh(colloc_data, init_params, learning_rate)

% this version calculates predictions for both association and temporal
% interaction potential using the structure location_details found by
% calling:
% location_details = get_all_separate_inputs(quality_matches_modis_time_loc_qa_mod470_aero440, quality_matches_misr_time_loc_qab_misr446b_misr446_aero440);

% the same as version without NIPS, only that the derivative is better
% very simple, since here there are no neighbourhoods

% function [CRF_pred, params, mae] = CRF_train_v1(data_truth, data_pred,
% data_last_week, init_params)
%
% � Trains CRF, calculates parameters and predicts times. Returns error.
% � Practically the same as the first version, only here we don't use
% � linear regression predictions (provided by Coric), but use history and
% � other info as suggested by professor.

% definitions

alpha1 = init_params;
alpha2 = init_params;
beta1 = init_params;

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
%%%%%%%%%%%%%%%% End of calculating QSpatial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%years=[2004 2005 2006 2007];
years=[2003 2004 2005 2007];
%years=[2003 2004 2005 2006];

total_years = 4;

gradient_alpha1_array=zeros(total_years,1);
gradient_alpha2_array=zeros(total_years,1);
gradient_beta1_array=zeros(total_years,1);

colloc_datas = cell(total_years,1);
colloc_datas{1} = load('/work/collocated_MISR_MODIS_AERONET_2003.mat');
colloc_datas{2} = load('/work/collocated_MISR_MODIS_AERONET_2004.mat');
colloc_datas{3} = load('/work/collocated_MISR_MODIS_AERONET_2005.mat');
colloc_datas{4} = load('/work/collocated_MISR_MODIS_AERONET_2007.mat');

colloc_indexes = cell(total_years,1);
uncolloc_indexes = cell(total_years,1);

len=zeros(4,1);
for i=1:total_years
   len(i) =  size(colloc_datas{i}.collocatedData, 1);
end

for i=1:total_years
    colloc_indexes{i}.colloc_index = zeros(len(i),1);
    index = 1;
    sizeX = 13;
    sizeY = 9;
    for j=1:len(i)
        data = colloc_datas{i}.collocatedData(j,:);
        colloc_indexes{i}.colloc_index(index,1) = (data(1,8) - 1)*sizeX*sizeY + (data(1,9)-1)*sizeX + data(1,10);
        index = index + 1;
    end
    uncolloc_indexes{i}.uncolloc_index = linearMatrix(:,4);
    uncolloc_indexes{i}.uncolloc_index(colloc_indexes{i}.colloc_index,:)=[];
end

QLL_Spatial = cell(total_years,1);
QLU_Spatial = cell(total_years,1);
QUL_Spatial = cell(total_years,1);
QUU_Spatial = cell(total_years,1);

for i=1:total_years
    QLL_Spatial{i} = QSpatial(colloc_indexes{i}.colloc_index,:);
    QLL_Spatial{i} = QLL_Spatial{i}(:,colloc_indexes{i}.colloc_index);
    QUU_Spatial{i} = QSpatial(uncolloc_indexes{i}.uncolloc_index,:);
    QUU_Spatial{i} = QUU_Spatial{i}(:,uncolloc_indexes{i}.uncolloc_index);
    QUL_Spatial{i} = QSpatial(uncolloc_indexes{i}.uncolloc_index,:);
    QUL_Spatial{i} = QUL_Spatial{i}(:,colloc_indexes{i}.colloc_index);
    QLU_Spatial{i} = QSpatial(colloc_indexes{i}.colloc_index,:);
    QLU_Spatial{i} = QLU_Spatial{i}(:,uncolloc_indexes{i}.uncolloc_index);
end
    
Q2=cell(total_years,1);
for i=1:total_years
     Q2{i} = (QLL_Spatial{i}-QLU_Spatial{i}*(QUU_Spatial{i}\QUL_Spatial{i}));
     Q2{i} = (QLL_Spatial{i}-QLU_Spatial{i}*(QUU_Spatial{i}\QUL_Spatial{i}));
     Q2{i} = (QLL_Spatial{i}-QLU_Spatial{i}*(QUU_Spatial{i}\QUL_Spatial{i}));
     Q2{i} = (QLL_Spatial{i}-QLU_Spatial{i}*(QUU_Spatial{i}\QUL_Spatial{i}));
end

% train CRF
while true
    tic;
        
    for i=1:total_years
        year = years(i);
        %load(['/work/collocated_MISR_MODIS_AERONET_' num2str(year) '.mat']);
        colloc_data = colloc_datas{i}.collocatedData;
        len = size(colloc_data, 1);
        truth_aeronet = colloc_data(:, 5);
        forecast_modis = colloc_data(:, 3);
        forecast_misr = colloc_data(:, 1);

        indicator_modis = (forecast_modis ~= 0);
        indicator_misr = (forecast_misr ~= 0);
        indicator_dummy = (forecast_misr == forecast_misr);

        Q1 = alpha1 * spdiags(indicator_misr, 0, len, len) + alpha2 * spdiags(indicator_modis, 0, len, len) + alpha3 * spdiags(indicator_dummy, 0, len, len);
      
        %Q2 = beta1*(QLL_Spatial{i}-QLU_Spatial{i}*(QUU_Spatial{i}\QUL_Spatial{i}));
        Q2_i = beta1*Q2{i};

        sigma_inv = 2 * (Q1 + Q2_i);
        
        b = 2 * (alpha1 * spdiags(indicator_misr, 0, len, len) * forecast_misr + alpha2 * spdiags(indicator_modis, 0, len, len) * forecast_modis);
   
        loc_prediction = sigma_inv\b;


        % association parameters alpha 
        gradient_alpha1 = -((truth_aeronet - loc_prediction)' * spdiags(indicator_misr, 0, len, len) * (truth_aeronet - loc_prediction)) + ...
            2 * (forecast_misr' - loc_prediction'*spdiags(indicator_misr, 0, len, len)) * (truth_aeronet - loc_prediction) + ...
            trace(sigma_inv\spdiags(indicator_misr, 0, len, len));

        gradient_alpha2 = -((truth_aeronet - loc_prediction)' * spdiags(indicator_misr, 0, len, len) * (truth_aeronet - loc_prediction)) + ...
            2 * (forecast_modis' - loc_prediction'*spdiags(indicator_modis, 0, len, len)) * (truth_aeronet - loc_prediction) + ...
            trace(sigma_inv\spdiags(indicator_modis, 0, len, len)); 

        % interaction parameter beta
        gradient_beta1 = ...
         - (truth_aeronet + loc_prediction)' * Q2_i/beta1 * (truth_aeronet - loc_prediction) + ...
         trace(sigma_inv\(Q2_i/beta1));

        gradient_alpha1_array(i) = gradient_alpha1;
        gradient_alpha2_array(i) = gradient_alpha2;
        gradient_beta1_array(i) = gradient_beta1;
        
    end
    
    % apply gradient information
    alpha1_new = exp(log(alpha1) + learning_rate * alpha1 * (sum(gradient_alpha1_array)));
    alpha2_new = exp(log(alpha2) + learning_rate * alpha2 * (sum(gradient_alpha2_array)));
    beta1_new = exp(log(beta1) + learning_rate * beta1 * (sum(gradient_beta1_array)));

    delta_alpha1 = abs(alpha1_new - alpha1);
    delta_alpha2 = abs(alpha2_new - alpha2);
    delta_beta1 = abs(beta1_new - beta1);
    
    alpha1 = alpha1_new;
    alpha2 = alpha2_new;
    beta1 = beta1_new;
    
    disp(delta_alpha1);
    
    %if (delta_alpha1 < 0.00000001 && delta_alpha2 < 0.00000001 && delta_beta1 < 0.00000001)
    if (delta_alpha1 < 0.00001 && delta_alpha2 < 0.00001 && delta_beta1 < 0.00001)
    %if (delta_alpha1 < 0.0000000001 && delta_alpha2 < 0.00000000001 && delta_beta1 < 0.00000000001)
        break;
    end;
end;