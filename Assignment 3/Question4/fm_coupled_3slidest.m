classdef fm_coupled_3slidest < forward_model    
    methods
        function obj = fm_coupled_3slidest(transform, m,n, angles, angles_2, angles_3)
            obj@forward_model(transform, m, n, angles,...
                angles_2,angles_3);
        end
        
        function product = mtimes(At, y)
            size_y = size(y);
            %Converting vectorized to matrix format
            %Iradon transform of the y
            %Transforming the vectors
            y1 = reshape(y(1:size_y/3), At.m, size(At.angles, 2));
            x1 = iradon(y1, At.angles, 'linear','Ram-Lak', 1, At.n);
            beta1 = At.transform(x1);
            
            y2 = reshape( y(1+size_y/3:2*size_y/3), At.m, size(At.angles_2, 2));
            d_x1 = iradon(y2, At.angles_2, 'linear','Ram-Lak', 1, At.n);
            d_beta1 = At.transform(d_x1);
            
            y3 = reshape(y(1+2*size_y/3:end), At.m, size(At.angles_3, 2));
            d_x2 = iradon(y3, At.angles_3, 'linear','Ram-Lak', 1, At.n);
            d_beta2 = At.transform(d_x2);
            
            product = [beta1(:) + d_beta1(:) + d_beta2(:);d_beta1(:);d_beta2(:)];
        end
    end
end
        