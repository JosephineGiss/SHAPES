classdef ScrewPropeller < classes.RigidBody
    % Screw Propeller from Happel & Brenner 1983
    % Low Reynolds number hydodynamics
    % Figure 5-4.1
    
    properties
        disk1  % first disk
        disk2  % second disque
        h      % semi-length of the connecting rod
    end

    methods

        function obj = ScrewPropeller(rad, h, theta, a, lambda, rho, position, setup)

            arguments 
                rad (1,:) double = 12
                h (1,1) double = 30
                theta (1,1) double {mustBeInRange(theta, 0, 1.5708)} = pi/3 % pi is not literal and cannot be used in arguments block
                a (1,1) double = 1;
                lambda (1,1) double = 2.5 % spacing factor
                rho (1,1) double = 1030
                position.xG (3,1) double = [0;0;0]
                position.rot (3,3) double = eye(3)
                setup.setup (1,1) double = 0
            end
            
            if isscalar(rad)
                rad = [rad rad];
            end
            
            Rot_disk1 = aux.rot_u([1 0 0], theta);   % rotation of the first disk 
            Rot_disk2 = aux.rot_u([1 0 0], -theta);  % rotation of the second disk 

            disk1 = classes.Disk(rad(1), a, lambda, rho, xG = [-h/2; 0; 0], rot = Rot_disk1, setup = setup.setup);
            disk2 = classes.Disk(rad(2), a, lambda, rho, xG = [h/2; 0; 0], rot = Rot_disk2, setup = setup.setup);
            
            xx_combined = [cell2mat(disk1.x), cell2mat(disk2.x)];
            
            N = size(xx_combined, 2);
            xx = zeros(3, N);
            
            % placing spheres
            Rot = aux.rot_u([0 1 0], pi/2); % so we have the formalism from H&B
            for k = 1:N
                xx(:,k) = Rot*xx_combined(:,k);
            end
            
            % readjusting position
            current_xG = classes.ConstrainedAggregate.getCenterOfMass(xx, a, rho);
            X_ID = mat2cell(xx - current_xG, 3, ones(1, N)); 
            for k = 1:N
                xx(:, k) = position.rot*(xx(:, k) - current_xG) + position.xG;
            end
            
            % instantiating object
            obj@classes.RigidBody(xx, a, rho, setup.setup);
            obj.disk1 = disk1;
            obj.disk2 = disk2;
            obj.X_ID = X_ID;
            obj.Rot_init = position.rot;

        end

    end


end