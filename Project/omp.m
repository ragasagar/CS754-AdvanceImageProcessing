%% OMP Implementation with niave implementation.
function [theta_star, T] = omp(A, y, noise,p)
    r = y;
    theta = zeros(size(A,2), 1);
    T = [];
    i = 0;
    normA = normc(A);
    while(norm(r)^2>noise && i<p)
       
       % Argmax problem
       [~ , j] = max(abs(r' *normA));
   
       T = [T j];
       i = i + 1;
       
       theta(T)  = pinv(normA(:,T))*y;

       r = y - normA(:, T)* theta(T);
    end
    
    theta_star = theta;
end