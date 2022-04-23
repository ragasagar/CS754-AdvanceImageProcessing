%% OMP implementation using QR_1 decomposition
function [theta_star, T] = omp_qr_1(A, u, noise, max_iter)
        r = u;
        theta_star = zeros(size(A,2),1);
        y = zeros(size(A,2),1);
        T = [];
        i = 1;
        normA = normc(A);
        Q = [];
        R=[];
        r = norm(r);
        while(r^2 > noise && i <= size(A,1))
            % Projection phase
             if i == 1
                 p = u' * normA;
             end
             
             %Selection step
             [~, j] = max(abs(p));
             T = [T, j];
              
             %Decompostion
             if i == 1
                 R = norm(normA(:,j));
                 Q = normA(:,j)/R;
             else
                 w = Q' * normA(:, j);
                 beta = sqrt(normA(:,j)'*normA(:,j) - w'*w);
                 Q = [Q (normA(:,j) - Q*w)/beta];
                 R = [R w; zeros(i-1, 1)' beta];
                 
                 r = sqrt(r'*r - (abs(u'*Q(:, i)))^2);
             end
             i = i+1;
             % Updating the projection vector
             p = p - ((Q(:, i-1)' * u ) * normA' * Q(:, i-1))';
             y(T) = R' \ (normA(:, T)' * u);
            theta_star(T) = R \  y(T);
        end
        
end