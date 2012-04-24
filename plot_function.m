function_array=zeros(150,150);

alpha1=1:300;
alpha2=1:300;
[X Y] = meshgrid(alpha1,alpha2);
        
iteration = 1;

iter_alpha1=1;
iter_alpha2=1;

for alpha1=1:300
    for alpha2=1:300
        truth_aeronet = colloc_data(:,5);
        forecast_misr = colloc_data(:,1);
        forecast_modis = colloc_data(:,3);

        indicator_modis = (forecast_modis ~= 0);
        indicator_misr = (forecast_misr ~= 0);

        [row column] = size(colloc_data);
        len = row;

        b = 2 * (alpha1 * spdiags(indicator_misr, 0, len, len) * forecast_misr + alpha2 * spdiags(indicator_modis, 0, len, len)* forecast_modis);

        Q1 = alpha1*spdiags(indicator_misr, 0, len, len) + alpha2*spdiags(indicator_modis, 0, len, len);

        sigma_inv = 2*Q1;

        u = sigma_inv\b;

        F = -0.5*(truth_aeronet-u)'*sigma_inv*(truth_aeronet-u) - 0.5*log(det(inv(sigma_inv)));
        
        function_array(iter_alpha1,iter_alpha2) = F;
        iter_alpha2 = iter_alpha2+1;
    end
    iter_alpha2=1;
    iter_alpha1 = iter_alpha1+1;
end

