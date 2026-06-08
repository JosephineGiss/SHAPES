classdef ChiralSwimmer_Rigid < classes.RigidBody
    
    properties
        helix;      % store properties of the helix
        rad;        % radius of the head
    end

    methods 
        
        function obj = ChiralSwimmer_Rigid(L, n, rad, a, lambda, rho, shape, position, setup)
            
            arguments
                L (1,1) double = 10
                n (1,1) double {mustBeInteger} = 3 
                rad (1,1) double = 1 % radius of the helix and the big sphere
                a (1,1) double = 0.1 % radius of the small spheres
                lambda (1,1) double = 2.5 % spacing factor
                rho (1,1) double = 1030
                shape.handedness (1,1) double = 1 % (1 == right-handed, -1 == left-handed)
                position.xG (3,1) double = [0;0;0]
                position.rot (3,3) double = eye(3)
                setup.setup (1,1) double = 0
            end
            helix = classes.CircularHelix(L, n, rad, a, lambda, rho, "handedness", shape.handedness);
            
            % placing spheres
            xx(:, 1) = [0; 0; -(rad + L/2)];
            xx(:, 2:helix.N + 1) = cell2mat(helix.x);

            % readjusting position
            current_xG = classes.ConstrainedAggregate.getCenterOfMass(xx, [rad, a*ones(1, helix.N)], rho);
            % X_ID = xx - current_xG;
            X_ID = mat2cell(xx - current_xG, 3, ones(1,helix.N + 1)); 
            for k = 1:helix.N + 1
                xx(:, k) = position.rot*(xx(:, k) - current_xG) + position.xG;
            end
            helix.xG = position.rot * (helix.xG - current_xG) + position.xG;
            
            % instantiating object
            obj@classes.RigidBody(xx, [rad, a*ones(1, helix.N)], rho, setup.setup);
            obj.helix = helix;
            obj.rad = rad;
            obj.X_ID = X_ID;
            
        end

    end

end