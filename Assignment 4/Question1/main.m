%% Intial values
clear;
m = 200;
n= 500;
r_m = .9*m;
v_m = .1*m;
close all;
addpath('./l1_ls_matlab');
%% Setting x value
x = zeros(n,1);
random_index = randperm(n,18);
x(random_index) = randi([0,1000],18,1);

A = (rand(200,500) < 0.5);
A = (2/sqrt(m))*A - (1/sqrt(m));


y = A*x;
sigma = 0.05 * sum(abs(A*x))/m;
y_noise = y + sigma*randn(size(y));

m_indices = 1:m;
r_indices  = randperm(m, r_m);
v_indices = setdiff(m_indices,  r_indices);
R = A(r_indices,:);
Rt = R';
V = A(v_indices,:);

y_r =y_noise(r_indices);
y_v = y_noise(v_indices);

lambdas = [0.0001 0.0005  0.001 0.005  0.01 0.05 0.1 0.5 1 2 5 10 15 20 30 50 100];
errors = zeros(size(lambdas,1));
rmses = zeros(size(lambdas,1));
for i = 1:size(lambdas,2)
    [x_pred, status] = l1_ls(R, Rt, m, n, y_r,lambdas(i), 0.01);
   
    %validation error check
    error = sum((y_v - V*x_pred)' * (y_v - V*x_pred))/size(V,1);
    errors(i) = error;
    
    %rmse error check
    rmses(i) = norm(x_pred-x)/norm(x);
end
figure;
plot(log(lambdas), errors );
set(gca, 'XTick',log(lambdas)); 
xtickangle(90)


title('VE vs logarithm lambdas');
ylabel('VE');
xlabel('log(lambda)');

figure;
plot(log(lambdas), rmses );
set(gca, 'XTick',log(lambdas)); 
xtickangle(90)


title('RMSE vs logarithm lambdas');
xlabel('log(lambda)');
ylabel('RMSE');

%% Coincident case:
R = A;
V = A;
y_r = y_noise;
y_v = y_noise;
errors_co = zeros(11,1);
rmses_co = zeros(11, 1);
for i = 1:size(lambdas, 2)
    [x_pred, status] = l1_ls(R,y_r,lambdas(i), 0.01);
    
    %validation error check
    y_pred = V*x_pred;
    error = sum((y_v - V*x_pred)' * (y_v - V*x_pred))/size(V,1);
    errors_co(i) = error;
    
    %rmse error check
   rmses_co(i) = norm(x_pred-x)/norm(x);
end
figure;
plot(log(lambdas), errors_co );
set(gca, 'XTick',log(lambdas));
xtickangle(90)


title('VE vs logarithm lambdas coincident set');
ylabel('VE');
xlabel('log(lambda)');

figure;
plot(log(lambdas), rmses_co );
set(gca, 'XTick',log(lambdas)) ;
xtickangle(90)


title('RMSE vs logarithm lambdas coincident set');
xlabel('log(lambda)');
ylabel('RMSE');

