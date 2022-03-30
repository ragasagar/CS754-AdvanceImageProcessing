%%
clear;
%% Reading images and perfroming ista
p=8;
x = double(imread("barbara256.png"));

rng(4);
dct_2d = dctmtx(p);
dct = kron(dct_2d', dct_2d');
phi = randn(32,64);
x_est = zeros(size(x));
count = zeros(size(x));
A = phi*dct;
eigen = max(eig(A'*A));
alpha = eigen + 2;
[row, col]= size(x);
for i = 1:row-p+1
    for j=1:col-p+1
        y_patch = x(i:i+p-1, j:j+p-1);
        y_patch  = phi*y_patch(:);% vectorizing the y_patch array
        lambda = 1;
        theta_est=zeros(p^2, 1);
        for iter = 1: 200
            theta_est = soft(theta_est + (1/alpha)*A'*(y_patch-A*theta_est), lambda/(2*alpha));
        end
        x_est_patch = dct*theta_est;
        x_est_patch = reshape(x_est_patch, p, p);
        x_est(i:i+p-1, j:j+p-1) = x_est(i:i+p-1, j:j+p-1) + x_est_patch;
        count(i:i+p-1, j:j+p-1) = count(i:i+p-1, j:j+p-1) + ones(p,p);
    end
end
x_est = x_est./count;
rmse = norm(x-x_est, 'fro')/norm(x, 'fro')

%%
figure;
subplot(1,2,1), imagesc (single (x)); 
title('Original Image')
colormap ('gray');

subplot(1,2,2), imagesc (single (x_est)); 
title('Reconstructed Image')
colormap ('gray');