classdef fm_coupled_2slidest < forward_model    
    methods
        function obj = fm_coupled_2slidest(transform,m,n, angles, angles_2)
            obj@forward_model(transform, m, n, angles,angles_2);
        end
        
        function product = mtimes(A, y)
            size_y = size(y);
            
            y1 = reshape(y(1:size_y/2), A.m, size(A.angles,2));
            y2 = reshape(y(0.5*size_y+1:end), A.m, size(A.angles_2,2));
            
            beta1 = iradon(y1, A.angles, 'linear','Ram-Lak', 1, A.n);
            delta_beta1 = iradon(y2, A.angles_2, 'linear', 'Ram-Lak', 1, A.n);
            
            x1 = A.transform(beta1);
            delta_x1 = A.transform(delta_beta1);
            
            product = [x1(:) + delta_x1(:); delta_x1(:)];
        end
    end
end