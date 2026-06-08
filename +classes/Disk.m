classdef Disk < classes.RigidBody

    properties (Access = public)
        rad;        % radius of the disk
        dAB;        % minimal distance between adjacent spheres
        N_circles;  % number of circles
        dr_eff;     % real spacement between circles
        dAB_eff;    % real spacement between spheres for each circle
    end

    methods 
        function obj = Disk(rad, a, lambda,  rho, position, setup)

            arguments 
                rad (1,1) double 
                a (1,1) double = 1;
                lambda (1,1) double = 2.5 % spacing factor
                rho (1,1) double = 1030
                position.xG (3,1) double = [0;0;0]
                position.rot (3,3) double = eye(3)
                setup.setup (1,1) double = 0
            end
            
            dAB = lambda*a; 
            i = 0; % count the number of spheres
            j = 0; % count the number of circles
            n = floor(rad/dAB + 1); % max number of circles
            M = cell(1,n); % spheres per circle
            xx = zeros(3, 5000); 
            dAB_eff = zeros(1, n);
            dr_eff = (rad)/(n - 1);
            for r = linspace(0, rad, n) % radius of circles
                j = j + 1;
                if r == 0 % one sphere at the center
                    i = i + 1;
                    xx(:,i) = [0;0;0];
                    M{j} = 1; 
                else
                    m = floor(2*pi/(acos(1 - dAB^2/(2*r^2))) + 1); % number of spheres par circle
                    % symmetry
                    if mod(m,2) == 0 % to have an even number of spheres per circle, we need m to be odd for now
                        m = m - 1;
                    end
                    thetas = linspace(0, 2*pi, m);
                    thetas(end) = []; % we remove the element of thetas which is 2pi because we already have an angle 0
                    for theta = thetas
                        i = i + 1;
                        xx(:,i) = [r*cos(theta); r*sin(theta); 0];
                    end
                    M{j} = m - 1; % Because we removed one element
                    dAB_eff(j) = norm([r*cos(thetas(1)); r*sin(thetas(1)); 0] - [r*cos(thetas(2)); r*sin(thetas(2)); 0]);
                end
            end
            

            N = i; % total number of small spheres used to make the disk
            xx = xx(:, 1:N);
            for k = 1:N
                xx(:, k) = position.rot*xx(:, k) + position.xG;
            end
            
            % instantiating object
            obj@classes.RigidBody(xx,a,rho,setup.setup);
            obj.N_circles = j; 
            obj.rad = rad;
            obj.dAB = dAB;
            obj.dr_eff = dr_eff;
            obj.dAB_eff = dAB_eff;
            
        end
    end

end