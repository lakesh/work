function [alpha1s, alpha2s, betas] = train_CRF_NIPS_v2(location_details, init_params, learning_rate)

% this version calculates predictions for both association and temporal
% interaction potential using the structure location_details found by
% calling:
% location_details = get_all_separate_inputs(quality_matches_modis_time_loc_qa_mod470_aero440, quality_matches_misr_time_loc_qab_misr446b_misr446_aero440);

% the same as version without NIPS, only that the derivative is better
% very simple, since here there are no neighbourhoods

% function [CRF_pred, params, mae] = CRF_train_v1(data_truth, data_pred,
% data_last_week, init_params)
%
%   Trains CRF, calculates parameters and predicts times. Returns error.
%   Practically the same as the first version, only here we don't use
%   linear regression predictions (provided by Coric), but use history and
%   other info as suggested by professor.

% definitions
num_iters = 20;

num_locations = length(location_details);

alpha1 = init_params;
alpha2 = init_params;
beta = init_params;

gradient_alpha1 = zeros(1, num_locations);
gradient_alpha2 = zeros(1, num_locations);
gradient_beta   = zeros(1, num_locations);

alpha1s = zeros(1, num_iters);
alpha2s = zeros(1, num_iters);
betas = zeros(1, num_iters);

num_points_modis = 0;
num_points_misr = 0;
num_points_aeronet = 0;
for loc = 1 : num_locations
    num_points_modis = location_details(loc).count_measurements_modis;
    num_points_misr = location_details(loc).count_measurements_misr;
    num_points_aeronet = location_details(loc).count_measurements_aeronet;
end;
modis_misr_ratio = num_points_modis / num_points_misr;

% train CRF

while true
    tic;
        collocatedData = colloc_datas{1}.collocated_misr_modis_aeronet;    
        [row column] = size(collocatedData);
        len = row;
        truth_aeronet = collocatedData(:,3);
        forecast_modis = collocatedData(:,2);
        forecast_misr = collocatedData(:,1);
        
        % indicators, if 0 then there is no prediction for that instance
        indicator_modis = (forecast_modis ~= 0);
        indicator_misr = (forecast_misr ~= 0);
        indicator_aeronet = (truth_aeronet ~= 0); 
        indicator_dummy = (truth_aeronet == truth_aeronet);
        
        alpha3 = 0.001;
        forecast_dummy = 0;
        Q1 = alpha1 * spdiags(indicator_modis, 0, len, len) + alpha2 * spdiags(indicator_misr, 0, len, len) +  alpha3 * spdiags(indicator_dummy, 0, len, len);
        
        
        %for i=1:row
        %    Q1(i,i) = Q1(i,i) + 1*0.001;
        %end
        %Q2 = beta * aaaa_create_sparse_tridiagonal_matrix_for_beta(len);

        sigma_inv = 2 * (Q1);
        b = 2 * (alpha1 * forecast_modis + alpha2 * forecast_misr);
%         sparse_tridiag_sigma = aaaa_inverse_tridiagonal_get_only_diagonals(sigma_inv);
        
        %loc_prediction = computeMi(sigma_inv, b);
        loc_prediction = sigma_inv\b;

        % association parameters alpha        
        gradient_alpha1 = -0.5 * sum(((truth_aeronet - loc_prediction) .* indicator_modis) .^ 2) + ...
            (2 * forecast_modis - (loc_prediction .* indicator_modis))' * (truth_aeronet - loc_prediction) + ...
            0.5 * computeTrace(sigma_inv, spdiags(indicator_modis, 0, len, len));
%             0.5 * diag(sparse_tridiag_sigma)' * indicator_modis;

        gradient_alpha2 = -0.5 * sum(((truth_aeronet - loc_prediction) .* indicator_misr) .^ 2) + ...
            (2 * forecast_misr - (loc_prediction .* indicator_misr))' * (truth_aeronet - loc_prediction) + ...
            0.5 * computeTrace(sigma_inv, spdiags(indicator_misr, 0, len, len));
%             0.5 * diag(sparse_tridiag_sigma)' * indicator_misr;

       
    
    alpha1_new = exp(log(alpha1) + learning_rate * alpha1 * (gradient_alpha1));
    alpha2_new = exp(log(alpha2) +  learning_rate * alpha2 * (gradient_alpha2));
    
    difference = abs(alpha1-alpha1_new) + abs(alpha2-alpha2_new);
    if difference <= 0.0001
        break;
    end
    
    disp(difference);
    
    % apply gradient information    
    alpha1 = exp(log(alpha1) + learning_rate * alpha1 * (gradient_alpha1));
    alpha2 = exp(log(alpha2) +  learning_rate * alpha2 * (gradient_alpha2));
    
   
    
    %alpha1s(iteration) = alpha1
    %alpha2s(iteration) = alpha2
    %iteration = iteration+1;
    
end;