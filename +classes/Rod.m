classdef Rod < classes.RigidBody
    
    properties
        L;          % length of the rod
        dAB;        % minimal distance between adjacent spheres
        dAB_eff;    % effective distance between adjacent spheres
    end
    
    methods
        
        function obj = Rod(L, a, lambda, rho, position, setup)
            arguments
                L (1,1) double = 30;
                a (1,1) double = 1;
                lambda (1,1) double = 2.5
                rho (1,1) double = 1030
                position.xG (3,1) double = [0;0;0]
                position.rot (3,3) double = eye(3)
                setup.setup (1,1) double = 0
            end
            
            dAB = lambda*a;
            n = floor(L/dAB) + 1; % number of spheres along the rod
            dAB_eff = L/(n - 1);
            z = linspace(0, L, n) - L/2; % center at z = 0
            if floor(n) == 1
                z = 0;
            end
            N = length(z);
            
            % placing spheres
            xx = zeros(3,N);
            for k = 1:N
                xx(:, k) = [0; 0; z(k)];
            end

            % readjusting position
            current_xG = classes.ConstrainedAggregate.getCenterOfMass(xx, a, rho);
            X_ID = mat2cell(xx - current_xG, 3, ones(1, N)); 
            for k = 1:N
                xx(: ,k) = position.rot*xx(:, k) + position.xG - current_xG;
            end
            
            % object initialization
            obj@classes.RigidBody(xx, a, rho, setup.setup);
            obj.L = L;
            obj.dAB = dAB;
            obj.dAB_eff = dAB_eff;
            obj.X_ID = X_ID;
            obj.Rot_init = position.rot;
        end

    end
        
end