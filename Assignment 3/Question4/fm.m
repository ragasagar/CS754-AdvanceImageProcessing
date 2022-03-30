classdef fm < forward_model
    
    methods
        function obj = fm(transform, m,n, angles)
            obj@forward_model(transform, m, n, angles);
        end
        
        function product = mtimes(A, x)
            beta = A.transform(reshape(x, A.n, A.n));
            product = radon(beta, A.angles);
            product = product(:);
        end
    end
end
            
            