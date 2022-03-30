classdef forward_model
    properties
        transform
        n
        m
        angles
        angles_2
        angles_3
    end
    
    methods
        function obj = forward_model(transform, m,n, angles, ...
                angles_2,angles_3)
            if nargin < 5
                angles_2 = zeros();
                angles_3 = zeros();
            elseif nargin < 6
                 angles_3 = zeros();
            end
            obj.transform = transform;
            obj.m = m;
            obj.n = n;
            obj.angles = angles;
            obj.angles_2 = angles_2;
            obj.angles_3 = angles_3;
        end

    end
    methods (Abstract)
        product = mtimes(A, x)
    end
end
            
            