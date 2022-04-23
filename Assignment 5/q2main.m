%% Constant values
close all;
clear;
clc;
m=[40 50 64 80 100 120];
n=128;
alpha_o=[0,3];
alpha = [0 1 2 3 4];
%% Main code
rng(0);
% Random orthogonal matrix
        % U = RandOrthMat(n);
 [U S V]= svd(rand(n));
 rmse = calculateError(m, n, alpha_o, U);
 figure;
 for i = 1:size(alpha_o,2)
     plot(m, log(rmse(:,i))+0.01, '-*'); 
     hold on
 end
xlabel('measurement size');
ylabel("log(RMSE)");
hold off;
legend('alpha 0', 'alpha 3');


%% comparing with increasing alpha values
 rmse = calculateError(m, n, alpha, U);
 figure;
 for i = 1:size(alpha,2)
     plot(m, log(rmse(:,i))+0.01, '-*'); 
     hold on
 end
xlabel('measurement size');
ylabel("log(RMSE)");
hold off;
legend('alpha 0', 'alpha 1', 'alpha 2', 'alpha 3', 'alpha 4');

%% Function for calculation
function rmse = calculateError(m, n,alpha,U)
    rmse = zeros(size(m,2), size(alpha,2));
    for i = 1: size(alpha, 2)
        A = diag((1:n).^-alpha(i));
        rng(1);
        sum_x = U*A*U';
        %processing each m with 10 signals
         for j = 1:size(m,2)
             for k = 1:10
                 x = sum_x*rand(n,1); %% Generate 10 signals from N (0, Σx).
                 phi  = sqrt(1/m(j)) * randn(m(j), n);
                 m_x = phi*x;
                 sigma = 0.01 * mean(abs(m_x));
                 y = m_x + sigma*rand(m(j),1);
                 reconst_x = (inv((phi'*phi)/(2*sigma^2) +  sum_x^-1/2))*phi'*y/(2*sigma^2);
                 rmse(j,i) = rmse(j,i) + norm(reconst_x - x)/norm(x);
             end
             rmse(j,i) = rmse(j,i)/10
         end
    end
end



