while true
    tic;
        [row column] = size(collocatedData);
        len = row;
        truth_aeronet = collocatedData(:,5);
        forecast_modis = collocatedData(:,3);
        forecast_misr = collocatedData(:,1);
        
        % indicators, if 0 then there is no prediction for that instance
        indicator_modis = (forecast_modis ~= 0);
        indicator_misr = (forecast_misr ~= 0);
        indicator_aeronet = (truth_aeronet ~= 0); 
        indicator_dummy = (truth_aeronet == truth_aeronet);
        
        alpha3 = 0.001;
        forecast_dummy = 0;
        Q1 = alpha1 * spdiags(indicator_modis, 0, len, len) + alpha2 * spdiags(indicator_misr, 0, len, len);
        
        

        sigma_inv = 2 * (Q1);
        b = 2 * (alpha1 * forecast_modis + alpha2 * forecast_misr);
%         sparse_tridiag_sigma = aaaa_inverse_tridiagonal_get_only_diagonals(sigma_inv);
        
        %loc_prediction = computeMi(sigma_inv, b);
        loc_prediction = (alpha1*forecast_modis+alpha2*forecast_misr)/(alpha1 + alpha2);

        % association parameters alpha        
        %gradient_alpha1 = -0.5 * (truth_aeronet-loc_prediction)'*2*spdiags(indicator_modis, 0, len, len)* (truth_aeronet-loc_prediction)+ ...
        %    (2 * forecast_modis' - loc_prediction'*2*spdiags(indicator_modis, 0, len, len)) * (truth_aeronet - loc_prediction) + ...
        %    0.5 * trace(sigma_inv\(2*(spdiags(indicator_modis, 0, len, len))));
%       %      0.5 * diag(sparse_tridiag_sigma)' * indicator_modis;

        %gradient_alpha2 = -0.5 * (truth_aeronet-loc_prediction)'*2*spdiags(indicator_misr, 0, len, len)* (truth_aeronet-loc_prediction)+ ...
        %    (2 * forecast_misr' - loc_prediction'*2*spdiags(indicator_misr, 0, len, len)) * (truth_aeronet - loc_prediction) + ...
        %    0.5 * trace(sigma_inv\(2*(spdiags(indicator_misr, 0, len, len))));
%             0.5 * diag(sparse_tridiag_sigma)' * indicator_misr;

        gradient_alpha1 = -sum((truth_aeronet - loc_prediction) .^ 2) + ...
        2 * (forecast_modis - loc_prediction)' * (truth_aeronet - loc_prediction) + ...
        0.5 * trace(sigma_inv\(2 * speye(len)));

        gradient_alpha2 = -sum((truth_aeronet - loc_prediction) .^ 2) + ...
            2 * (forecast_misr - loc_prediction)' * (truth_aeronet - loc_prediction) + ...
            0.5 * trace(sigma_inv\(2 * speye(len)));
    
        alpha1_new = exp(log(alpha1) + learning_rate * alpha1 * (gradient_alpha1));
        alpha2_new = exp(log(alpha2) +  learning_rate * alpha2 * (gradient_alpha2));

        difference = abs(alpha1-alpha1_new) + abs(alpha2-alpha2_new);
        if difference <= 0.0001
            break;
        end
        
        pause
        
        disp(difference);

        % apply gradient information    
        alpha1 = exp(log(alpha1) + learning_rate * alpha1 * (gradient_alpha1));
        alpha2 = exp(log(alpha2) +  learning_rate * alpha2 * (gradient_alpha2));
    toc;
end