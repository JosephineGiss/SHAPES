classdef ChiralDipole < classes.RigidBody

    properties
        L;
        rod;
        right_handed_helix;
        left_handed_helix;
    end
    
    methods
        
        function obj = ChiralDipole(L_rod, L_helix, n, alpha, a, lambda, rho, position, setup)
            
            arguments 
                L_rod (1,:) double = 20
                L_helix (1,1) double = 10
                n (1,1) double = 3
                alpha (1,1) = 1
                a (1,1) double = 1
                lambda (1,1) double = 2.5 % spacing factor
                rho (1,1) double = 1030
                position.xG (3,1) double = [0;0;0]
                position.rot (3,3) double = eye(3)
                setup.setup (1,1) double = 0
            end
            
            rod = classes.Rod(L_rod, a, lambda);

            xGR = rod.x{end} + [0; 0; rod.dAB_eff*rod.a{end}];
            xGL = rod.x{1} - [0; 0; rod.dAB_eff*rod.a{1}];

            rot_helix1 = aux.rot_u([0 1 0], pi);
            rot_helix2 = aux.rot_u([0 0 1], pi);
           
            right_handed_helix = classes.ConicalHelix(L_helix, n, alpha, a, lambda, 0, handedness = 1, xG = xGR);
            left_handed_helix = classes.ConicalHelix(L_helix, n, alpha, a, lambda, 0, handedness = -1, xG = xGL, rot = rot_helix2*rot_helix1);

            xx_combined = [cell2mat(rod.x), cell2mat(right_handed_helix.x), cell2mat(left_handed_helix.x)];
            [xmin, xmax] = bounds(xx_combined(3,:));
            L = xmax - xmin;

            N = size(xx_combined, 2);
            xx = zeros(3, N);
            
            % readjusting position
            for k = 1:N
                xx(:,k) = position.rot*xx_combined(:,k) + position.xG;
            end
            
            % instantiating object
            obj@classes.RigidBody(xx, a, rho, setup.setup);
            obj.L = L;
            obj.rod = rod;
            obj.right_handed_helix = right_handed_helix;
            obj.left_handed_helix = left_handed_helix;

        end

    end

end