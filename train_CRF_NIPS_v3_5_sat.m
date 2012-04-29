function [alpha1s, alpha2s, alpha3s, alpha4s, alpha5s, betas] = train_CRF_NIPS_v3_5_sat(location_details, init_params, learning_rate)
% logdet(A) = 2 * sum(log(diag(chol(A)))); 
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
num_iters = 2000;
threshold = 0.0001;
regularization = 1;

num_locations = length(location_details);

alpha1 = init_params(1);
alpha2 = init_params(2);
alpha3 = init_params(3);
alpha4 = init_params(4);
alpha5 = init_params(5);
beta = init_params(6);

gradient_alpha1 = zeros(1, num_locations);
gradient_alpha2 = zeros(1, num_locations);
gradient_alpha3 = zeros(1, num_locations);
gradient_alpha4 = zeros(1, num_locations);
gradient_alpha5 = zeros(1, num_locations);
gradient_beta   = zeros(1, num_locations);
alpha1s = zeros(1, num_iters);
alpha2s = zeros(1, num_iters);
alpha3s = zeros(1, num_iters);
alpha4s = zeros(1, num_iters);
alpha5s = zeros(1, num_iters);
betas = zeros(1, num_iters);

% num_points_modis = 0;
% num_points_misr = 0;
% num_points_aeronet = 0;
% for loc = 1 : num_locations
%     num_points_modis = location_details(loc).count_measurements_modis;
%     num_points_misr = location_details(loc).count_measurements_misr;
%     num_points_aeronet = location_details(loc).count_measurements_aeronet;
% end;
% modis_misr_ratio = 1;%num_points_modis / num_points_misr;

% train CRF
load good_perm;
for iteration = 1 : num_iters  
    tic;
    count_NA_USA = 0;
    count_NA_USA_valid = 0;
    
    objective_func_first_part = 0;
    objective_func_second_part = 0;
    blah = [];
    for loc = 1 : num_locations
        % take only NE USA locations
        geo_loc = location_details(1, loc).location;
        if ((geo_loc(1) < -70) && (geo_loc(1) > -79) && (geo_loc(2) > 38.5) && (geo_loc(2) < 43))
            count_NA_USA = count_NA_USA + 1;
        else
            continue;
        end;
        
        num_points_modis = location_details(loc).count_measurements_modis;
        num_points_misr = location_details(loc).count_measurements_misr;
        num_points_omi = location_details(loc).count_measurements_omi;
        num_points_seawifs = location_details(loc).count_measurements_seawifs;
        num_points_caliop = location_details(loc).count_measurements_caliop;
        num_points_aeronet = location_details(loc).count_measurements_aeronet;
        
        if (num_points_aeronet == 0)
            % if no true readings here then move on
            continue;
        else
            if (sum([num_points_modis num_points_misr num_points_omi num_points_seawifs num_points_caliop] ~= 0) < 2)
                % if only one satellite here then move on
                continue;
            end;
            count_NA_USA_valid = count_NA_USA_valid + 1;
        end;
%         blah = [blah; loc];
        if (isempty(find(temp_perm(1 : 50) == count_NA_USA_valid, 1)))
            continue;
        end;
%         
%         num_points_modis1(count_NA_USA_valid) = location_details(loc).count_measurements_modis;
%         num_points_misr1(count_NA_USA_valid) = location_details(loc).count_measurements_misr;
%         num_points_omi1(count_NA_USA_valid) = location_details(loc).count_measurements_omi;
%         num_points_seawifs1(count_NA_USA_valid) = location_details(loc).count_measurements_seawifs;
%         num_points_caliop1(count_NA_USA_valid) = location_details(loc).count_measurements_caliop;
%         num_points_aeronet1(count_NA_USA_valid) = location_details(loc).count_measurements_aeronet;
        
        len = location_details(loc).total_days;
        truth_aeronet = location_details(loc).measurements_aeronet;
        forecast_modis = location_details(loc).measurements_modis;
        forecast_misr = location_details(loc).measurements_misr;
        forecast_omi = location_details(loc).measurements_omi;
        forecast_seawifs = location_details(loc).measurements_seawifs;
        forecast_caliop = location_details(loc).measurements_caliop;
        
        % indicators, if 0 then there is no prediction for that instance
        indicator_modis = (forecast_modis ~= 0);
        indicator_misr = (forecast_misr ~= 0);
        indicator_omi = (forecast_omi ~= 0);
        indicator_seawifs = (forecast_seawifs ~= 0);
        indicator_caliop = (forecast_caliop ~= 0);
        indicator_labeled = (truth_aeronet ~= 0);
        indicator_unlabeled = (truth_aeronet == 0);
        
        Q1 = alpha1 * spdiags(indicator_modis, 0, len, len) + alpha2 * spdiags(indicator_misr, 0, len, len) + alpha3 * spdiags(indicator_omi, 0, len, len) + ...
            alpha4 * spdiags(indicator_seawifs, 0, len, len) + alpha5 * spdiags(indicator_caliop, 0, len, len);
        Q2 = beta * aaaa_create_sparse_tridiagonal_matrix_for_beta(len);

        sigma_inv = 2 * (Q1 + Q2);
        b = 2 * (alpha1 * forecast_modis + alpha2 * forecast_misr + alpha3 * forecast_omi + alpha4 * forecast_seawifs + alpha5 * forecast_caliop);
        
        Q_LL = sigma_inv(indicator_labeled, indicator_labeled);
        Q_UU = sigma_inv(indicator_unlabeled, indicator_unlabeled);
        Q_LU = sigma_inv(indicator_labeled, indicator_unlabeled);
        Q_UL = sigma_inv(indicator_unlabeled, indicator_labeled);
        
%         sparse_tridiag_sigma = aaaa_inverse_tridiagonal_get_only_diagonals(sigma_inv);
        
        loc_prediction = computeMi(sigma_inv, b);
        mju_star = loc_prediction(indicator_labeled);
        b_star = b(indicator_labeled);
        sigma_star_inv = Q_LL - (Q_LU / Q_UU) * Q_UL;
        
        objective_func_first_part(count_NA_USA_valid) = (truth_aeronet(indicator_labeled) - mju_star)' * sigma_star_inv * (truth_aeronet(indicator_labeled) - mju_star);
%         objective_func_second_part = objective_func_second_part * det(inv(sigma_inv));
%         objective_func_second_part(count_NA_USA_valid) = log(det(sigma_star_inv));    % this faster
        objective_func_second_part(count_NA_USA_valid) = 2 * sum(log(diag(chol(sigma_star_inv)))); 
        
        % association parameters alpha
        alpha1_sigma_inv = 2 * spdiags(indicator_modis, 0, len, len);
        alpha1_sigma_inv_derivative = alpha1_sigma_inv(indicator_labeled, indicator_labeled) - (alpha1_sigma_inv(indicator_labeled, indicator_unlabeled) / Q_UU) * Q_UL + ...
            (Q_LU / Q_UU) * (alpha1_sigma_inv(indicator_unlabeled, indicator_unlabeled) / Q_UU) * Q_UL - (Q_LU / Q_UU) * alpha1_sigma_inv(indicator_unlabeled, indicator_labeled);        
        gradient_alpha1(loc) = -0.5 * (truth_aeronet(indicator_labeled) - mju_star)' * alpha1_sigma_inv_derivative * (truth_aeronet(indicator_labeled) - mju_star) + ...
            (2 * forecast_modis(indicator_labeled)' - mju_star' * alpha1_sigma_inv_derivative) * (truth_aeronet(indicator_labeled) - mju_star) + ...
            0.5 * trace(sigma_star_inv \ alpha1_sigma_inv_derivative);
%             0.5 * diag(sparse_tridiag_sigma)' * indicator_modis;

        alpha2_sigma_inv = 2 * spdiags(indicator_misr, 0, len, len);
        alpha2_sigma_inv_derivative = alpha2_sigma_inv(indicator_labeled, indicator_labeled) - (alpha2_sigma_inv(indicator_labeled, indicator_unlabeled) / Q_UU) * Q_UL + ...
            (Q_LU / Q_UU) * (alpha2_sigma_inv(indicator_unlabeled, indicator_unlabeled) / Q_UU) * Q_UL - (Q_LU / Q_UU) * alpha2_sigma_inv(indicator_unlabeled, indicator_labeled);        
        gradient_alpha2(loc) = -0.5 * (truth_aeronet(indicator_labeled) - mju_star)' * alpha2_sigma_inv_derivative * (truth_aeronet(indicator_labeled) - mju_star) + ...
            (2 * forecast_misr(indicator_labeled)' - mju_star' * alpha2_sigma_inv_derivative) * (truth_aeronet(indicator_labeled) - mju_star) + ...
            0.5 * trace(sigma_star_inv \ alpha2_sigma_inv_derivative);

        alpha3_sigma_inv = 2 * spdiags(indicator_omi, 0, len, len);
        alpha3_sigma_inv_derivative = alpha3_sigma_inv(indicator_labeled, indicator_labeled) - (alpha3_sigma_inv(indicator_labeled, indicator_unlabeled) / Q_UU) * Q_UL + ...
            (Q_LU / Q_UU) * (alpha3_sigma_inv(indicator_unlabeled, indicator_unlabeled) / Q_UU) * Q_UL - (Q_LU / Q_UU) * alpha3_sigma_inv(indicator_unlabeled, indicator_labeled);        
        gradient_alpha3(loc) = -0.5 * (truth_aeronet(indicator_labeled) - mju_star)' * alpha3_sigma_inv_derivative * (truth_aeronet(indicator_labeled) - mju_star) + ...
            (2 * forecast_omi(indicator_labeled)' - mju_star' * alpha3_sigma_inv_derivative) * (truth_aeronet(indicator_labeled) - mju_star) + ...
            0.5 * trace(sigma_star_inv \ alpha3_sigma_inv_derivative);

        alpha4_sigma_inv = 2 * spdiags(indicator_seawifs, 0, len, len);
        alpha4_sigma_inv_derivative = alpha4_sigma_inv(indicator_labeled, indicator_labeled) - (alpha4_sigma_inv(indicator_labeled, indicator_unlabeled) / Q_UU) * Q_UL + ...
            (Q_LU / Q_UU) * (alpha4_sigma_inv(indicator_unlabeled, indicator_unlabeled) / Q_UU) * Q_UL - (Q_LU / Q_UU) * alpha4_sigma_inv(indicator_unlabeled, indicator_labeled);        
        gradient_alpha4(loc) = -0.5 * (truth_aeronet(indicator_labeled) - mju_star)' * alpha4_sigma_inv_derivative * (truth_aeronet(indicator_labeled) - mju_star) + ...
            (2 * forecast_seawifs(indicator_labeled)' - mju_star' * alpha4_sigma_inv_derivative) * (truth_aeronet(indicator_labeled) - mju_star) + ...
            0.5 * trace(sigma_star_inv \ alpha4_sigma_inv_derivative);

        alpha5_sigma_inv = 2 * spdiags(indicator_caliop, 0, len, len);
        alpha5_sigma_inv_derivative = alpha5_sigma_inv(indicator_labeled, indicator_labeled) - (alpha5_sigma_inv(indicator_labeled, indicator_unlabeled) / Q_UU) * Q_UL + ...
            (Q_LU / Q_UU) * (alpha5_sigma_inv(indicator_unlabeled, indicator_unlabeled) / Q_UU) * Q_UL - (Q_LU / Q_UU) * alpha5_sigma_inv(indicator_unlabeled, indicator_labeled);        
        gradient_alpha5(loc) = -0.5 * (truth_aeronet(indicator_labeled) - mju_star)' * alpha5_sigma_inv_derivative * (truth_aeronet(indicator_labeled) - mju_star) + ...
            (2 * forecast_caliop(indicator_labeled)' - mju_star' * alpha5_sigma_inv_derivative) * (truth_aeronet(indicator_labeled) - mju_star) + ...
            0.5 * trace(sigma_star_inv \ alpha5_sigma_inv_derivative);

        beta_sigma_inv = 2 * aaaa_create_sparse_tridiagonal_matrix_for_beta(len);
        beta_sigma_inv_derivative = beta_sigma_inv(indicator_labeled, indicator_labeled) - (beta_sigma_inv(indicator_labeled, indicator_unlabeled) / Q_UU) * Q_UL + ...
            (Q_LU / Q_UU) * (beta_sigma_inv(indicator_unlabeled, indicator_unlabeled) / Q_UU) * Q_UL - (Q_LU / Q_UU) * beta_sigma_inv(indicator_unlabeled, indicator_labeled); 
        gradient_beta(loc) = -0.5 * (truth_aeronet(indicator_labeled) + mju_star)' * beta_sigma_inv_derivative * (truth_aeronet(indicator_labeled) - mju_star) + ...
            0.5 * trace(sigma_star_inv \ beta_sigma_inv_derivative);
    end;
    
    objective_func = - 0.5 * sum(objective_func_first_part) + 0.5 * sum(objective_func_second_part);
    
    % apply gradient information    
    alpha1_new = exp(log(alpha1) + learning_rate * alpha1 * (sum(gradient_alpha1) - regularization * alpha1));
    alpha2_new = exp(log(alpha2) + learning_rate * alpha2 * (sum(gradient_alpha2) - regularization * alpha2));
    alpha3_new = exp(log(alpha3) + learning_rate * alpha3 * (sum(gradient_alpha3) - regularization * alpha3));
    alpha4_new = exp(log(alpha4) + learning_rate * alpha4 * (sum(gradient_alpha4) - regularization * alpha4));
    alpha5_new = exp(log(alpha5) + learning_rate * alpha5 * (sum(gradient_alpha5) - regularization * alpha5));
    beta_new = exp(log(beta) + learning_rate * beta * (sum(gradient_beta) - regularization * beta));
    
    delta_alpha1 = abs(alpha1_new - alpha1);
    delta_alpha2 = abs(alpha2_new - alpha2);
    delta_alpha3 = abs(alpha3_new - alpha3);
    delta_alpha4 = abs(alpha4_new - alpha4);
    delta_alpha5 = abs(alpha5_new - alpha5);
    delta_beta = abs(beta_new - beta);
    alpha1 = alpha1_new;
    alpha2 = alpha2_new;
    alpha3 = alpha3_new;
    alpha4 = alpha4_new;
    alpha5 = alpha5_new;
    beta = beta_new;
    
    objective_func
    [alpha1, alpha2, alpha3, alpha4, alpha5, beta]
    [delta_alpha1, delta_alpha2, delta_alpha3, delta_alpha4, delta_alpha5, delta_beta]
    
    if ((delta_alpha1 < threshold) && (delta_alpha2 < threshold) && (delta_alpha3 < threshold) && (delta_alpha4 < threshold) && (delta_alpha5 < threshold) && (delta_beta < threshold))
        break;
    end;
    
    alpha1s(iteration) = alpha1;
    alpha2s(iteration) = alpha2;
    alpha3s(iteration) = alpha3;
    alpha4s(iteration) = alpha4;
    alpha5s(iteration) = alpha5;
    betas(iteration) = beta;
    toc;
end;