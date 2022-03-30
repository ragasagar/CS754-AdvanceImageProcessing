classdef fm_coupled_3slides < forward_model
    methods
        function obj = fm_coupled_3slides(transform, m, n, angles, angles_2, angles_3)
            obj@forward_model(transform, m, n, angles,...
                angles_2,angles_3);

        end
        
        function product = mtimes(A, x)
            size_x = size(x);            
            beta1 = reshape( x(1:size_x/3), A.n, A.n);
            d_beta1 = reshape(x(size_x/3+1:2*size_x/3), A.n, A.n);
            d_beta2 = reshape(x(2*size_x/3+1:end), A.n, A.n);
            
            beta1 = A.transform(beta1);
            beta1_deta1 = A.transform(d_beta1);
            beta2_delta2 = A.transform(d_beta2);
            
            R1x1 = radon(beta1, A.angles);
            R2x1 = radon(beta1, A.angles_2);
            R3x1 = radon(beta1, A.angles_3);
            R2delta_x1 = radon(beta1_deta1, A.angles_2);
            R3delta_x2 = radon(beta2_delta2, A.angles_3);
            
            product = [R1x1(:); R2x1(:) + R2delta_x1(:); R3x1(:)+ R3delta_x2(:)];
        end
    end
end
            
        