function questionb(s50, s51)
    s50_size = size(s50,1);
    lambda = 0.01;
    rel_tol = 0.01;
    % creating  parallel beam tomographic projections at
    % at any 18 randomly angle
    rng(4)
    r_angles=unifrnd(0,180,1,18);
    m_slice_50 = radon(s50, r_angles);
    m_slice_51 = radon(s51, r_angles);
    addpath('./l1_ls_matlab');
   
    x50 = s50(:) ;%vectorized form of signal
    n = size(x50,1);
    y_50 = m_slice_50(:);
    m = size(y_50, 1);
    s50_size = size(s50,1);
    projection_size = size(m_slice_50,1);
    A = fm(@idct2, projection_size, s50_size, r_angles);
    At = fmt(@dct2, projection_size, s50_size, r_angles);
    [x_pred, status] = l1_ls(A, At, m, n, y_50,lambda,rel_tol);
    reconst_ics_s50 = idct2(reshape(x_pred, s50_size, s50_size));
    figure();
    imshow(reconst_ics_s50);
    colormap('gray');
    title("Independent CS-based Reconstruction - Slice 50");
    %s51 processing
    x51 = s51(:) ;%vectorized form of signal
    n = size(x51,1);
    y_51 = m_slice_51(:);
    m = size(y_51, 1);
    projection_size = size(m_slice_51,1);
    A = fm(@idct2, projection_size, s50_size, r_angles);
    At = fmt(@dct2, projection_size, s50_size, r_angles);
    [x_pred_51, status] = l1_ls(A, At, m, n, y_51,lambda,rel_tol);
    reconst_ics_s51 = idct2(reshape(x_pred_51, s50_size, s50_size));
    figure();
    imshow(reconst_ics_s51);
    colormap('gray');
    title("Independent CS-based Reconstruction - Slice 51");
end