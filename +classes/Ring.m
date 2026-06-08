classdef Ring < classes.RigidBody

    properties
        rad; % radius of the ring
        dAB; % minimal distance between adjacent spheres
    end
    
    methods
        
        function obj = Ring(rad, a, lambda, rho, position, setup)
            
            arguments 
                    rad (1,1) double = 1
                    a (1,1) double = 0.1
                    lambda (1,1) double = 2.5
                    rho (1,1) double = 1030
                    position.xG (3,1) double = [0;0;0]
                    position.rot (3,3) double = eye(3)
                    setup.setup (1,1) double = 0
            end
            
            dAB = lambda*a;

            m = floor(2*pi/(acos(1 - dAB^2/(2*rad^2))) + 1);
            % symmetry
            if mod(m,2) == 0 
                m = m - 1;
            end
            
            thetas = linspace(0, 2*pi, m);
            thetas(end) = [];
            
            N = length(thetas);
            
            xx = zeros(3,N);
            for k = 1:N
                xx(:,k) = [rad*sin(thetas(k)); 0; rad*cos(thetas(k))];
            end

            % readjusting position
            current_xG = classes.ConstrainedAggregate.getCenterOfMass(xx, a, rho);
            for k = 1:N
                xx(:, k) = position.rot*xx(:, k) + position.xG - current_xG;
            end

            % object initialization
            obj@classes.RigidBody(xx, a, rho, setup.setup);

            obj.rad = rad;
            obj.dAB = dAB;
            
        end

    end

end