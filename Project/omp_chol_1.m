%% OMP implementation using CHOLesky decomposition
function [theta_star, T] = omp_chol_1(A, u, noise, max_iter)
        r  = u;
        theta_star = zeros(size(A,2), 1);
        T = [];
        i = 1;
        normA = normc(A);
        V = [];
        y=zeros(size(A,1), 1);
        T = [];
        while(norm(r)^2> noise)
            % Projection phase
             if i == 1
                 p = (normA' * r);
             end
             
             %Selection step
             [~, j] = max(abs(p));
             T = [T, j];
             
             %Decompostion
             if i == 1
                 V = (1);
             else 
                 s = normA(:,T(1: end-1))' * normA(:, j);
                 z = V\s;
                 V = [V zeros(i-1, 1); z' sqrt(1 - z'*z)];
             end
             
             y(T) = V \ (normA(:, T)'*u);
             theta_star(T) = V' \ y(T);
                 
             r = u - normA(:, T)* theta_star(T);
                 
             p = (r' * normA);
             i = i+1;
        end 
end