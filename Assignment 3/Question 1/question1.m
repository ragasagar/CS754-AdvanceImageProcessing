%%
clear;
%% Reading images and perfroming ista
p=8;
x = imread("barbara256.png");
x = double(x);

rng(0);
n = 3*randn(size(x));
y = x+n;

dct_2d = dctmtx(p);
dct = kron(dct_2d', dct_2d');
phi = eye(p^2);
x_est = zeros(size(x));
count = zeros(size(x));
A = phi*dct;
eigen = max(eig(A'*A));
alpha = eigen + 1;
[row, col]= size(x);
for i = 1:row-p+1
    for j=1:col-p+1
        y_patch = y(i:i+p-1, j:j+p-1);
        y_patch = y_patch(:);% vectorizing the y_patch array
        lambda = 1;
        theta_est=zeros(p^2, 1);
        for iter = 1: 100
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
error = norm(x- y, 'fro')/norm(x, 'fro')

%% Displaying the images
figure;
subplot(1,3,1), imagesc (single (x)); 
title('Original Image')
colormap ('gray');

subplot(1,3,2), imagesc (single (y)); 
title('Noisy Image')
colormap ('gray');
axis tight;

subplot(1,3,3), imagesc (single (x_est)); 
title('Reconstructed Image')
colormap ('gray');
