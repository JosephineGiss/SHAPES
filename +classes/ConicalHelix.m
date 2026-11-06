classdef ConicalHelix < classes.RigidBody

    properties
        handedness;
        L;          % length of the spiral along the z-axis
        n;          % number of full turns
        alpha;      % radius of the helix
        beta;       % pitch/(2*pi)
        dAB;        % minimal distance between adjacent spheres
        dAB_eff;    % effective distance between adjacent spheres
    end
    
    methods
        
        function obj = ConicalHelix(L, n, alpha, a, lambda, rho, shape, position, setup)
            
            arguments
                L (1,1) double = 15
                n (1,1) double = 3 
                alpha (1,1) double = 0.5 % arg from helix parameterization
                a (1,1) double = 1 
                lambda (1,1) double = 2.5 % spacing factor
                rho (1,1) double = 1030
                shape.handedness (1,1) double = 1 % (1 == right-handed, -1 == left-handed)
                shape.T0 (1,1) double = 0
                position.xG (3,1) double = [0;0;0]
                position.rot (3,3) double = eye(3)
                setup.setup (1,1) double = 0
            end
            
            dAB = lambda*a;
            beta = L/(2*pi*n); 
            T = 2*pi*n;
            theta = linspace(shape.T0, T, 10000);
            v = sqrt(alpha^2*(1 + theta.^2) + beta^2);
            s = cumtrapz(theta, v); % curvilinear abscissa: s(i) = int_theta(1)^theta(end) v*dt; 
            Lc = s(end) - s(1); % curvilinear length
            Ns = floor(Lc/dAB + 1); % maximum number of spheres along the helix
            s_uniform = linspace(0, Lc, Ns); 
            t = interp1(s, theta, s_uniform); % t values to get uniform distribution of spheres along the curve

            % placing spheres
            xx = zeros(3, Ns - 1);
            for k = 1:Ns - 1 % to get center of mass aligned with the z-axis in case of full turns
                xx(:, k) = [alpha*t(k)*cos(t(k)), shape.handedness*alpha*t(k)*sin(t(k)), beta*t(k)];
            end
            dAB_eff = norm(xx(:, end) - xx(:, end - 1));
            
            % readjusting position
            current_xG = classes.ConstrainedAggregate.getCenterOfMass(xx, a, rho);
            X_ID = mat2cell(xx - current_xG, 3, ones(1, Ns - 1)); 
            for k = 1:Ns - 1
                xx(:, k) = position.rot*(xx(:, k) - current_xG) + position.xG;
            end
            
            % instantiating object
            obj@classes.RigidBody(xx, a, rho, setup.setup);
            obj.handedness = shape.handedness;
            obj.L = L;
            obj.alpha = alpha;
            obj.n = n;
            obj.beta = beta;
            obj.dAB = dAB;
            obj.dAB_eff = dAB_eff;
            obj.X_ID = X_ID;
            obj.Rot_init = position.rot;

        end

            
        
    end
    

end