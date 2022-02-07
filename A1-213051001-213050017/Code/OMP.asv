%% OMP Implementation
function [theta_star, T] = OMP(A, y, noise,p)
    r = y;
    theta = zeros(size(A,2), 1);
    T = [];
    i = 0;
    normA = normc(A);
    while(norm(r)^2>noise && i<p)
       
       % Argmax problem
       [~ , j] = max(abs(r' *normA));
       
       T = union(T, j);
       i = i + 1;
       
       theta(T)  = pinv(A(:,T))*y;

       r = y - A(:, T)* theta(T);
    end
    
    theta_star = theta;
end