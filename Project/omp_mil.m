%% OMP with matrix inversion lemma method
function [theta_star, T] = omp_mil(A, u, noise, max_iter)
        r  = u;
        theta_star = zeros(size(A,2),1);
        T = [];
        i = 1;
        normA = normc(A);
        psi = [];
        delta = [];
        r = norm(r);
        while(r^2 > noise && i < size(A, 1))
            
            % Projection step
            if i == 1
                p0 = u' * normA;
                p = p0;
            else
                if i==2
                    h = [];
                end
                G = normA' * normA;
                p = p + ( delta * G(:, T) * [h' -1]' )' ;
            end
            
            % Selection of the index part
             [val_k, j] = max(abs(p));
             T = [T, j];
             
             % Updation steps
             if i == 1
                 lambda = 1/(normA(:,j)'*normA(:,j));
                 delta = lambda * val_k;
                 psi = normA(:,T) * (normA(:, T)' * normA(:, T))^-1;
                 theta_star(T) = psi' * u;
             else
                 h = psi' * normA( :,  j);
                 temp = normA(:, T( 1 : i-1))*h ;
                 lambda = 1/ ((normA(:,j)'*normA(:,j)) - temp' * temp );
                 v = normA(:,j) - normA(:, T( 1 : i-1))*h;
                 psi = [ (psi - (lambda * v)*h'),  lambda*v];

                 %Calculating the delta value:
                 delta = lambda * ( p0(j) - h' * p0(T(1: i-1) )' ); 
                 r  = sqrt( r^2 + (abs(delta)^2)/lambda) - (2* real(delta * p(j)) );
                 theta_star(T) = [theta_star(T(1: i-1)); 0 ] - delta*[h; -1];
             end
             i = i+1;
        end
        
end