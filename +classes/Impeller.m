classdef Impeller < classes.RigidBody 
    % Impeller from Happel & Brenner 1983
    % Low Reynolds number hydodynamics
    % Figure 5-4.4
    
    properties
        disk1  % first disk
        disk2  % second disque
        h      % semi-length of the connecting rod
    end

    methods

        function obj = Impeller(rad, h, theta, a, lambda, rho, position, setup)

            arguments 
                rad (1,:) double = 12 
                h (1,1) double = 30
                theta (1,1) double {mustBeInRange(theta, 0, 1.5708)} = pi/4 % pi is not literal and cannot be used in arguments block
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

            Rot_disk1 = aux.rot_u([0 1 0], theta);   % rotation of the first disk 
            Rot_disk2 = aux.rot_u([0 1 0], -theta);  % rotation of the second disk 

            disk1 = classes.Disk(rad(1), a, lambda, rho, xG = [-h/2; 0; 0], rot = Rot_disk1, setup = setup.setup);
            disk2 = classes.Disk(rad(2), a, lambda, rho, xG = [h/2; 0; 0], rot = Rot_disk2, setup = setup.setup);
            
            xx_combined = [cell2mat(disk1.x), cell2mat(disk2.x)];
            
            N = size(xx_combined, 2);
            xx = zeros(3, N);
            
            % placing spheres
            Rot1 = aux.rot_u([0 0 1], pi/2);
            Rot2 = aux.rot_u([0 1 0], -pi/2);
            Rot = Rot2*Rot1; % so we have the formalism from H&B
            for k = 1:N
                xx(:,k) = Rot*xx_combined(:,k);
            end
            
            % readjusting position
            for k = 1:N
                xx(:,k) = position.rot*xx(:,k) + position.xG;
            end

            % instantiating object
            obj@classes.RigidBody(xx, a, rho, setup.setup);
            obj.disk1 = disk1;
            obj.disk2 = disk2;
            obj.h = h;

        end

    end


end