classdef fmt < forward_model
    methods
        function obj = fmt(transform, m,n, angles)
            obj@forward_model(transform, m, n, angles);
        end
        
        function product = mtimes(At, y)
            y = reshape(y, At.m, size(At.angles,2));
            
            result  =iradon(y, At.angles, 'linear','Ram-Lak', 1, At.n);
            x = At.transform(result);
            product = x(:);
        end
    end
end
            
            