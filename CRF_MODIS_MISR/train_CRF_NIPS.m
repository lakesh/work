function [alpha1, alpha2, prediction] = train_CRF_NIPS(truth, forecast1, forecast2, init_params, learning_rate)

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
num_iters = 1000;
len = length(forecast1);

alpha1 = init_params;
alpha2 = init_params;

prediction = computeMju(alpha1, alpha2, forecast1, forecast2);

% train CRF
for iterations = 1 : num_iters  
    Q1 = alpha1 + alpha2;
    Q2 = zeros(1);      
    
    sigma_inv = 2 * (Q1 + Q2);
    sigma = inv(sigma_inv);

    % association parameters alpha        
    gradient = - sum((truth - prediction) .^ 2) + 2 * sum((forecast1 - prediction) .* (truth - prediction)) + len * sigma;
    log_alpha1 = log(alpha1) + learning_rate * alpha1 * gradient;
    alpha1 = exp(log_alpha1);

    gradient = - sum((truth - prediction) .^ 2) + 2 * sum((forecast2 - prediction) .* (truth - prediction)) + len * sigma;
    log_alpha2 = log(alpha2) + learning_rate * alpha2 * gradient;
    alpha2 = exp(log_alpha2);
    
    prediction = computeMju(alpha1, alpha2, forecast1, forecast2);
end;

function mju = computeMju(alpha1, alpha2, forecast1, forecast2)
mju = (alpha1 * forecast1 + alpha2 * forecast2) / (alpha1 + alpha2);
return;