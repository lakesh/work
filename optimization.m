% train CRF
function [f df] = optimization(x)

    alpha1=x(1);
    alpha2=x(2);
    beta1=x(3);
    total_years=1;
    
    gradient_alpha1_array=zeros(total_years,1);
    gradient_alpha2_array=zeros(total_years,1);
    f_array=zeros(total_years,1);
    gradient_beta1_array=zeros(total_years,1);
    
    global colloc_datas labelled_indexes unlabelled_indexes len alpha3 QSpatial;
    for i=1:total_years
    
        colloc_data = colloc_datas{i}.collocated_misr_modis_aeronet;
        labelled = find(colloc_data(:,3) ~= 0);

        truth_aeronet = colloc_data(labelled, 3);
        forecast_modis = colloc_data(labelled, 2);
        forecast_misr = colloc_data(labelled, 1);


        indicator_modis = (forecast_modis ~= 0);
        indicator_misr = (forecast_misr ~= 0);

        indicator_modis_whole = (colloc_data(:,2) ~= 0);
        indicator_misr_whole = (colloc_data(:,1) ~= 0);
        indicator_dummy_whole = (colloc_data(:,2) == colloc_data(:,2));

        Q1 = alpha1 * spdiags(indicator_misr_whole, 0, len, len) + alpha2 * spdiags(indicator_modis_whole, 0, len, len) + alpha3 * spdiags(indicator_dummy_whole, 0, len, len);
        Q2 = beta1*QSpatial;

       % sigma_inv = 2 * (Q1);
        sigma_inv = 2 * (Q1+Q2);

        Q_LL = sigma_inv(labelled_indexes{i}.labelled_index,:);
        Q_LL = Q_LL(:,labelled_indexes{i}.labelled_index);
        Q_UU = sigma_inv(unlabelled_indexes{i}.unlabelled_index, :);
        Q_UU = Q_UU(:,unlabelled_indexes{i}.unlabelled_index);
        Q_LU = sigma_inv(labelled_indexes{i}.labelled_index, :);
        Q_LU = Q_LU(:, unlabelled_indexes{i}.unlabelled_index);
        Q_UL = sigma_inv(unlabelled_indexes{i}.unlabelled_index,:);
        Q_UL = Q_UL(:, labelled_indexes{i}.labelled_index);

        sigma_star_inv = Q_LL - (Q_LU / Q_UU) * Q_UL;

        labelled_size = size(forecast_misr,1);
        b = 2 * (alpha1 * spdiags(indicator_misr, 0, labelled_size, labelled_size) * forecast_misr + alpha2 * spdiags(indicator_modis, 0, labelled_size, labelled_size) * forecast_modis);

        loc_prediction = sigma_star_inv\b;


        % association parameters alpha 
        alpha1_sigma_inv = 2 * spdiags(indicator_misr_whole, 0, len, len);
        alpha1_sigma_inv_LL = alpha1_sigma_inv(labelled_indexes{i}.labelled_index,:);
        alpha1_sigma_inv_LL = alpha1_sigma_inv_LL(:,labelled_indexes{i}.labelled_index);
        alpha1_sigma_inv_LU = alpha1_sigma_inv(labelled_indexes{i}.labelled_index,:);
        alpha1_sigma_inv_LU = alpha1_sigma_inv_LU(:,unlabelled_indexes{i}.unlabelled_index);
        alpha1_sigma_inv_UL = alpha1_sigma_inv(unlabelled_indexes{i}.unlabelled_index,:);
        alpha1_sigma_inv_UL = alpha1_sigma_inv_UL(:,labelled_indexes{i}.labelled_index);
        alpha1_sigma_inv_UU = alpha1_sigma_inv(unlabelled_indexes{i}.unlabelled_index,:);
        alpha1_sigma_inv_UU = alpha1_sigma_inv_UU(:,unlabelled_indexes{i}.unlabelled_index);

        alpha2_sigma_inv = 2 * spdiags(indicator_modis_whole, 0, len, len);
        alpha2_sigma_inv_LL = alpha2_sigma_inv(labelled_indexes{i}.labelled_index,:);
        alpha2_sigma_inv_LL = alpha2_sigma_inv_LL(:,labelled_indexes{i}.labelled_index);
        alpha2_sigma_inv_LU = alpha2_sigma_inv(labelled_indexes{i}.labelled_index,:);
        alpha2_sigma_inv_LU = alpha2_sigma_inv_LU(:,unlabelled_indexes{i}.unlabelled_index);
        alpha2_sigma_inv_UL = alpha2_sigma_inv(unlabelled_indexes{i}.unlabelled_index,:);
        alpha2_sigma_inv_UL = alpha2_sigma_inv_UL(:,labelled_indexes{i}.labelled_index);
        alpha2_sigma_inv_UU = alpha2_sigma_inv(unlabelled_indexes{i}.unlabelled_index,:);
        alpha2_sigma_inv_UU = alpha2_sigma_inv_UU(:,unlabelled_indexes{i}.unlabelled_index);

        beta1_sigma_inv = 2*QSpatial;
        beta1_sigma_inv_LL = beta1_sigma_inv(labelled_indexes{i}.labelled_index,:);
        beta1_sigma_inv_LL = beta1_sigma_inv_LL(:,labelled_indexes{i}.labelled_index);
        beta1_sigma_inv_LU = beta1_sigma_inv(labelled_indexes{i}.labelled_index,:);
        beta1_sigma_inv_LU = beta1_sigma_inv_LU(:,unlabelled_indexes{i}.unlabelled_index);
        beta1_sigma_inv_UL = beta1_sigma_inv(unlabelled_indexes{i}.unlabelled_index,:);
        beta1_sigma_inv_UL = beta1_sigma_inv_UL(:,labelled_indexes{i}.labelled_index);
        beta1_sigma_inv_UU = beta1_sigma_inv(unlabelled_indexes{i}.unlabelled_index,:);
        beta1_sigma_inv_UU = beta1_sigma_inv_UU(:,unlabelled_indexes{i}.unlabelled_index);


        %alpha1_sigma_inv_derivative = alpha1_sigma_inv_LL - (alpha1_sigma_inv_LU / Q_UU) * Q_UL + ...
        %    (Q_LU / Q_UU) * (alpha1_sigma_inv_UU / Q_UU) * Q_UL - (Q_LU / Q_UU) * alpha1_sigma_inv_UL;        
        alpha1_sigma_inv_derivative = alpha1_sigma_inv_LL - (Q_LU/Q_UU)*alpha1_sigma_inv_UL - ((alpha1_sigma_inv_LU  - (Q_LU/Q_UU) * alpha1_sigma_inv_UU) / Q_UU) * Q_UL;

        gradient_alpha1 = -0.5 * ((truth_aeronet - loc_prediction)' * alpha1_sigma_inv_derivative * (truth_aeronet - loc_prediction)) + ...
            (2 * forecast_misr' - loc_prediction'*alpha1_sigma_inv_derivative) * (truth_aeronet - loc_prediction) + ...
            0.5 * trace(sigma_star_inv \ alpha1_sigma_inv_derivative);
        clear alpha1_sigma_inv_LL alpha1_sigma_inv_LU alpha1_sigma_inv_UL alpha1_sigma_inv_UU

        %[gradient_alpha1]

        %alpha2_sigma_inv_derivative = alpha2_sigma_inv_LL - (alpha2_sigma_inv_LU / Q_UU) * Q_UL + ...
        %    (Q_LU / Q_UU) * (alpha2_sigma_inv_UU / Q_UU) * Q_UL - (Q_LU / Q_UU) * alpha2_sigma_inv_UL;        
        alpha2_sigma_inv_derivative = alpha2_sigma_inv_LL - (Q_LU/Q_UU)*alpha2_sigma_inv_UL - ((alpha2_sigma_inv_LU  - (Q_LU/Q_UU) * alpha2_sigma_inv_UU) / Q_UU) * Q_UL;

        gradient_alpha2 = -0.5 * ((truth_aeronet - loc_prediction)' * alpha2_sigma_inv_derivative * (truth_aeronet - loc_prediction)) + ...
            (2 * forecast_modis' - loc_prediction'*alpha2_sigma_inv_derivative) * (truth_aeronet - loc_prediction) + ...
            0.5 * trace(sigma_star_inv \ alpha2_sigma_inv_derivative);
        %[gradient_alpha2]
        clear alpha2_sigma_inv_LL alpha2_sigma_inv_LU alpha2_sigma_inv_UL alpha2_sigma_inv_UU

        % interaction parameter beta
        %beta1_sigma_inv_derivative = beta1_sigma_inv_LL - (beta1_sigma_inv_LU / Q_UU) * Q_UL + ...
        %   (Q_LU / Q_UU) * (beta1_sigma_inv_UU / Q_UU) * Q_UL - (Q_LU / Q_UU) * beta1_sigma_inv_UL; 
        beta1_sigma_inv_derivative = beta1_sigma_inv_LL - (Q_LU/Q_UU)*beta1_sigma_inv_UL - ((beta1_sigma_inv_LU  - (Q_LU/Q_UU) * beta1_sigma_inv_UU) / Q_UU) * Q_UL;

        gradient_beta1 = ...
         - 0.5 * (truth_aeronet + loc_prediction)' * beta1_sigma_inv_derivative * (truth_aeronet - loc_prediction) + ...
         0.5 * trace(sigma_star_inv \ beta1_sigma_inv_derivative);
        %[gradient_beta1]
        
        clear beta1_sigma_inv_LL beta1_sigma_inv_LU beta1_sigma_inv_UL beta1_sigma_inv_UU
        
        gradient_alpha1_array(i) = gradient_alpha1;
        gradient_alpha2_array(i) = gradient_alpha2;
        gradient_beta1_array(i) = gradient_beta1;     
        f_array(i) =  -(-0.5*(truth_aeronet-loc_prediction)'*sigma_star_inv*(truth_aeronet-loc_prediction) + 0.5*(2 * sum(log(diag(chol(sigma_star_inv))))));
        
    end
    
    f = sum(f_array);%-(-0.5*(truth_aeronet-loc_prediction)'*sigma_star_inv*(truth_aeronet-loc_prediction) + 0.5*(2 * sum(log(diag(chol(sigma_star_inv))))));
    gradient_alpha1 = sum(gradient_alpha1_array);
    gradient_alpha2 = sum(gradient_alpha2_array);
    gradient_beta1 = sum(gradient_beta1_array);
    
    df = [gradient_alpha1 gradient_alpha2 gradient_beta1];
    [alpha1 alpha2 beta1]
    %[gradient_alpha1_array(1) gradient_alpha1_array(2) gradient_alpha1_array(3) gradient_alpha1_array(4)]
    
    

end

