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
%   Trains CRF, calculates parameters and predicts times. Returns error.
%   Practically the same as the first version, only here we don't use
%   linear regression predictions (provided by Coric), but use history and
%   other info as suggested by professor.

% definitions
num_iters = 10000;


alpha1_array=[];
alpha2_array=[];

alpha1 = init_params;
alpha2 = init_params;
% beta = init_params;

colloc_data = colloc_datas{1}.collocated_misr_modis_aeronet;
alpha3 = 0.001;
len = size(colloc_data, 1);
iter = 0;
% train CRF
while true
    tic;
    truth_aeronet = colloc_data(:, 3);
    forecast_modis = colloc_data(:, 2);
    forecast_misr = colloc_data(:, 1);
    
    indicator_modis = (forecast_modis ~= 0);
    indicator_misr = (forecast_misr ~= 0);
    indicator_dummy = (forecast_misr == forecast_misr);
 
    Q1 = alpha1 * spdiags(indicator_misr, 0, len, len) + alpha2 * spdiags(indicator_modis, 0, len, len) + alpha3 * spdiags(indicator_dummy, 0, len, len);
    %Q2 = beta * aaaa_create_sparse_tridiagonal_matrix_for_beta(len);
    Q2=0;
    
    sigma_inv = 2 * (Q1 + Q2);
    b = 2 * (alpha1 * spdiags(indicator_misr, 0, len, len) * forecast_misr + alpha2 * spdiags(indicator_modis, 0, len, len) * forecast_modis);
%         sparse_tridiag_sigma = aaaa_inverse_tridiagonal_get_only_diagonals(sigma_inv);

    loc_prediction = sigma_inv\b;
    

    % association parameters alpha        
    gradient_alpha1 = -((truth_aeronet - loc_prediction)' * spdiags(indicator_misr, 0, len, len) * (truth_aeronet - loc_prediction)) + ...
        2 * (forecast_misr' - loc_prediction'*spdiags(indicator_misr, 0, len, len)) * (truth_aeronet - loc_prediction) + ...
        trace(sigma_inv\spdiags(indicator_misr, 0, len, len));

    gradient_alpha2 = -((truth_aeronet - loc_prediction)' * spdiags(indicator_modis, 0, len, len) * (truth_aeronet - loc_prediction)) + ...
        2 * (forecast_modis' - loc_prediction'*spdiags(indicator_modis, 0, len, len)) * (truth_aeronet - loc_prediction) + ...
        trace(sigma_inv\spdiags(indicator_modis, 0, len, len));

    % interaction parameter beta
    %gradient_beta = ...
    % - (truth_aeronet + loc_prediction)' * aaaa_create_sparse_tridiagonal_matrix_for_beta(len) * (truth_aeronet - loc_prediction) + ...
    %    trace(sigma_inv\ aaaa_create_sparse_tridiagonal_matrix_for_beta(len));
 
 
    % apply gradient information    
    alpha1_new = exp(log(alpha1) + learning_rate * alpha1 * (gradient_alpha1));
    alpha2_new = exp(log(alpha2) + learning_rate * alpha2 * (gradient_alpha2));
    %beta_new = exp(log(beta) + learning_rate * beta * (gradient_beta-0.01*beta));
    
    delta_alpha1 = abs(alpha1_new-alpha1);
    delta_alpha2 = abs(alpha2_new-alpha2);
    %delta_beta = abs(beta_new - beta);
    
    
    alpha1 = alpha1_new;
    alpha2 = alpha2_new;
    %beta = beta_new;
    
    
    [delta_alpha1 delta_alpha2 alpha1 alpha2]
    %pause
    %alpha1_array(iteration) = alpha1;
    %alpha2_array(iteration) = alpha2;
    
    %if (delta_alpha1 < 0.00000001 && delta_alpha2 < 0.00000001)
    %    break;
    %end;
    %if (delta_alpha1 < 0.0000000001 && delta_alpha2 < 0.00000000001 && delta_beta < 0.00000000001)
    %    break;
    %end;
    if (delta_alpha1 < 0.000000001 && delta_alpha2 < 0.000000001) %&& delta_beta < 0.000001)
       break; 
    end
    
    iter = iter + 1;
    
%     [truth_aeronet'; forecast_modis'; forecast_misr'; loc_prediction']
 	
%     pause;
%     beta = exp(log(beta) + learning_rate * beta * gradient_beta);
end;