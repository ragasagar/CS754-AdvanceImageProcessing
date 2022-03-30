function questionc(s50, s51,s52)
    
    %Random angles between o to pie 
    r_angles_set1 = unifrnd(0, 180, 1, 18);
    r_angles_set2 = unifrnd(0, 180, 1, 18);
    r_angles_set3 = unifrnd(0,180,1,18);
   
    %Measurment matrix for radon transform
    m_slice_50 = radon(s50, r_angles_set1);
    m_slice_51 = radon(s51, r_angles_set2);
    
    %Vectorized y vector
    y_50 = m_slice_50(:);
    y_51 = m_slice_51(:);
    
    %coupled measurement
    y_coupled = [y_50; y_51];
    m = size(y_coupled, 1);
    n = size(s50(:), 1) + size(s51(:), 1);
    data_size = size(s50,1);
    projection_size = size(m_slice_50, 1);
    
    A = fm_coupled_2slides(@idct2, projection_size, data_size, r_angles_set1, r_angles_set2);
    At = fm_coupled_2slidest(@dct2, projection_size, data_size, r_angles_set1, r_angles_set2);
    lambda = 0.01;
    rel_tol = 0.01;
    
    [beta, status] = l1_ls(A, At, m, n, y_coupled, lambda,rel_tol);
    beta1 = beta(1:n/2);
    delta_beta1 = beta(1+n/2:end);
    beta2 = beta1 + delta_beta1;
    
    coupled_r_2_s50 = idct2(reshape(beta1, data_size, data_size));
    coupled_r_2_s51 = idct2(reshape(beta2, data_size, data_size));
    figure();
    montage([coupled_r_2_s50,coupled_r_2_s51]);
    colormap('gray');
    title("Coupled CS-based Reconstruction - Slice 50 and slice51");
    
    
    %3slides
    m_slice_52 = radon(s52, r_angles_set3);
    y_52 = m_slice_52(:);
   
    % coupling all three slices
    y_coupled = [y_50; y_51; y_52];
    m = size(y_coupled, 1);
    n = size(s50(:), 1) + size(s51(:), 1) + size(s52(:), 1);
    
    %Initializing forwardmodel, 2 objects for fm and fm transpose
    A = fm_coupled_3slides(@idct2, projection_size, data_size, r_angles_set1, r_angles_set2, r_angles_set3);
    At = fm_coupled_3slidest(@dct2, projection_size, data_size, r_angles_set1, r_angles_set2, r_angles_set3);
    
    %calculating value using library
    [beta, status] = l1_ls(A, At, m, n, y_coupled, lambda, rel_tol);
    beta1 = beta(1:n/3);
    delta_beta1 = beta(1+n/3:2*n/3);
    delta_beta2 = beta(1+2*n/3:end);
    beta2 = beta1 + delta_beta1;
    beta3 = beta1 + delta_beta1 + delta_beta2;
    
    coupled_r_3_s50 = idct2(reshape(beta1, data_size, data_size));
    coupled_r_3_s51 = idct2(reshape(beta2, data_size, data_size));
    coupled_r_3_s52 = idct2(reshape(beta3, data_size, data_size));
    
    figure();
    montage([coupled_r_3_s50,coupled_r_3_s51,coupled_r_3_s52]);
    colormap('gray');
    title("Coupled CS-based Reconstruction - Slice 50, slice51, slice52");
end