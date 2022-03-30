classdef fm_coupled_2slides < forward_model
    methods
        function obj = fm_coupled_2slides(transform,m,n, angles, angles_2)
            obj@forward_model(transform, m, n, angles,angles_2);
        end
        
        function product = mtimes(A,x)
            size_x = size(x);
            
            x1 = reshape(x(1:size_x/2), A.n, A.n);
            delta_x1 = reshape(x(0.5*size_x+1:end), A.n, A.n);
            
            beta1 = A.transform(x1);
            delta_beta1 = A.transform(delta_x1);
            
            R1U_beta1 = radon(beta1, A.angles);
            R2U_beta1 = radon(beta1, A.angles_2);
            
            R2U_delta_beta1 = radon(delta_beta1, A.angles_2);
            product = [R1U_beta1(:); R2U_beta1(:) + R2U_delta_beta1(:)];
        end
    end
end