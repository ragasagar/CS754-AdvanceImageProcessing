%% OMP implementation using CHOLesky decomposition
function [theta_star, T] = omp_chol_2(A, u, noise, max_iter)
        r  = u;
        theta_star = zeros(size(A,2),1);
        T = [];
        i = 1;
        normA = normc(A);
        V = [];
        y = zeros(size(A,2),1);
        
        r = norm(r);
        while (r^2 > noise)
            % Projection phase
             if i == 1
                 p0 = u' * normA ;
                 p=p0;
             else
                %chol 2 modification
                G = normA' * normA;
                p = p0 - (G(:,T) * theta_star(T))';
             end
             
             %Selection step
             [~, j] = max(abs(p));
             T = [T j];
             
             %Decompostion
              if i == 1
                 V = (1);
             else
                 s = normA(:,T(1: i-1))' * normA(:, j);
                 z = V \ s;
                 V = [V zeros(i-1, 1); z' sqrt(1 - z'*z)];
             end
             
             y(T) = V \ (normA(:, T)'*u);
             
             theta_star(T) = V' \ y(T);
             r = sqrt(u'*u - y(T)'*y(T));
             
             i = i+1;
        end
end