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
% beta = init_params;

alpha3 = 0.001;


years=[2004 2005 2006 2007];
total_years = 4;

gradient_alpha1_array=zeros(total_years,1);
gradient_alpha2_array=zeros(total_years,1);
% train CRF
while true
    tic;
    
    for i=1:total_years
        year = years(i);
        load(['/work/collocated_MISR_MODIS_AERONET_' num2str(year) '.mat']);
        colloc_data = collocatedData;
        len = size(colloc_data, 1);
        truth_aeronet = colloc_data(:, 5);
        forecast_modis = colloc_data(:, 3);
        forecast_misr = colloc_data(:, 1);

        indicator_modis = (forecast_modis ~= 0);
        indicator_misr = (forecast_misr ~= 0);
        indicator_dummy = (forecast_misr == forecast_misr);

        Q1 = alpha1 * spdiags(indicator_misr, 0, len, len) + alpha2 * spdiags(indicator_modis, 0, len, len) + alpha3 * spdiags(indicator_dummy, 0, len, len);
        Q2 = 0;%beta * aaaa_create_sparse_tridiagonal_matrix_for_beta(len);

        sigma_inv = 2 * (Q1 + Q2);
        b = 2 * (alpha1 * spdiags(indicator_misr, 0, len, len) * forecast_misr + alpha2 * spdiags(indicator_modis, 0, len, len) * forecast_modis);
    % � � � � sparse_tridiag_sigma = aaaa_inverse_tridiagonal_get_only_diagonals(sigma_inv);

        loc_prediction = sigma_inv\b;


        % association parameters alpha � � � �
        gradient_alpha1 = -((truth_aeronet - loc_prediction)' * spdiags(indicator_misr, 0, len, len) * (truth_aeronet - loc_prediction)) + ...
            2 * (forecast_misr' - loc_prediction'*spdiags(indicator_misr, 0, len, len)) * (truth_aeronet - loc_prediction) + ...
            trace(sigma_inv\spdiags(indicator_misr, 0, len, len));

        gradient_alpha2 = -((truth_aeronet - loc_prediction)' * spdiags(indicator_misr, 0, len, len) * (truth_aeronet - loc_prediction)) + ...
            2 * (forecast_modis' - loc_prediction'*spdiags(indicator_modis, 0, len, len)) * (truth_aeronet - loc_prediction) + ...
            trace(spdiags(sigma_inv\indicator_modis, 0, len, len));

        % interaction parameter beta
    % � � gradient_beta = ...
    % � � � � - ((truth_aeronet + loc_prediction) .* indicator_aeronet)' * aaaa_create_sparse_tridiagonal_matrix_for_beta(len) * ((truth_aeronet - loc_prediction) .* indicator_aeronet) + ...
    % � � � � 0.5 * computeTrace(sigma_inv, 2 * aaaa_create_sparse_tridiagonal_matrix_for_beta(len));

        gradient_alpha1_array(i) = gradient_alpha1;
        gradient_alpha2_array(i) = gradient_alpha2;
        
    end
    
    % apply gradient information � �
    alpha1_new = exp(log(alpha1) + learning_rate * alpha1 * (sum(gradient_alpha1_array)-0.01*alpha1));
    alpha2_new = exp(log(alpha2) + learning_rate * alpha2 * (sum(gradient_alpha2_array)-0.01*alpha2));

    delta_alpha1 = abs(alpha1_new - alpha1);
    delta_alpha2 = abs(alpha2_new - alpha2);
    alpha1 = alpha1_new;
    alpha2 = alpha2_new;

    disp(delta_alpha1);

    if (delta_alpha1 < 0.0000000001 && delta_alpha2 < 0.00000000001)
        break;
    end;
end;