classdef CircularArc < classes.RigidBody

    properties
        rad;        % radius of the circular arc
        thetaf;     % angle of the circular arc
        theta;      % angles used to place spheres
        dAB;        % minimal distance between adjacent spheres
        dAB_eff;    % effective distance between adjacent spheres
    end
    
    methods
        
        function obj = CircularArc(rad, a, thetaf, lambda, rho, position, setup)
            
            arguments
                rad (1,1) double = 1
                a (1,1) double = 0.05
                thetaf (1,1) double = pi
                lambda (1,1) double = 2.3
                rho (1,1) double = 1030
                position.xG (3,1) double = [0;0;0]
                position.rot (3,3) double = eye(3)
                setup.setup (1,1) double = 0
            end
            
            dAB = lambda*a; % minimal distance between adjacent spheres
            theta_min = dAB/rad; % minimal angle between the center of two spheres
            N = floor(thetaf/theta_min + 1); % maximal number of spheres
            theta = linspace(0,thetaf,N);

            % placing spheres
            xx = zeros(3, N);
            for k = 1:N
                xx(:, k) = [rad*cos(theta(k)); rad*sin(theta(k)); 0];
            end
            dAB_eff = norm(xx(:, 2) - xx(:, 1));
            
            % readjusting position
            current_xG = classes.ConstrainedAggregate.getCenterOfMass(xx, a, rho);
            X_ID = mat2cell(xx - current_xG, 3, ones(1, N));
            for k = 1:N
                xx(:, k) = position.rot*(xx(:, k) - current_xG) + position.xG;
            end
            
            % instantiating object
            obj@classes.RigidBody(xx, a, rho, setup.setup);
            obj.rad = rad;
            obj.thetaf = thetaf;
            obj.theta = num2cell(theta);
            obj.dAB = dAB;
            obj.dAB_eff = dAB_eff;
            obj.X_ID = X_ID;
            obj.Rot_init = position.rot;
            
        end
        
    end

end