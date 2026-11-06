classdef UBentDisk < classes.Disk
    
    properties
        Rc; % radius of curvature (R = inf for flat disk, smaller R = more curved)
    end
    
    methods
        function obj = UBentDisk(rad, Rc, a, lambda, rho, position, setup)
            
            arguments 
                rad (1,1) double = 12
                Rc (1,:) double = [6.3 5.8 6.9]
                a (1,1) double = 0.2367/2
                lambda (1,1) double = 5
                rho (1,1) double = 1030
                position.xG (3,1) double = [0;0;0]
                position.rot (3,3) double = eye(3)
                setup.setup (1,1) double = 0
            end
            
            if any(Rc <= 0)
                error('Radius of curvature R must be positive ou infinite');
            end
            
            % Calling constructor => flat disk initialization
            % we do not compute tensors for the disk
            guard = classes.SphereAssembly.suspendCalc(); %#ok<NASGU>
            obj@classes.Disk(rad, a, lambda, rho, setup = setup.setup);
            
            obj.Rc = Rc;
            if ~isscalar(Rc)
                points = linspace(0, obj.rad, length(Rc));
                poly_Rc = getRadiusOfCurvatureFit(obj, points, Rc);
                obj.Rc = poly_Rc;
            end

            % apply bending if Rc is not infinite
            if ~isinf(Rc)
                if isscalar(Rc)
                    xx = obj.applyBending(cell2mat(obj.x), Rc);
                else
                    xx = obj.applyBending(cell2mat(obj.x), Rc, poly_Rc);
                end
            end
            
            % xx0 = xx;
            % for k = 1:obj.N
            %     xx(:,k) = position.rot*xx(:,k) + position.xG - sum([obj.m{:}].*xx0./obj.M, 2);
            % end
            % readjusting position
            current_xG = classes.ConstrainedAggregate.getCenterOfMass(xx, a, rho);
            X_ID = mat2cell(xx - current_xG, 3, ones(1, obj.N)); 
            for k = 1:obj.N
                xx(:, k) = position.rot*(xx(:, k) - current_xG) + position.xG;
            end

            
            %%% Updating space properties

            if obj.setup == 1 
                % checking all the spheres are on the same side of the wall
                wall_conditions = [nnz(xx(1,:)>0) == obj.N, nnz(xx(1,:)<0) == obj.N];
                if not(or(wall_conditions(1), wall_conditions(2)))
                    % either all x components of the spheres position are positive
                    % or all x components of the spheres position are negative
                    % if not -> error
                    error('All the spheres must be on the same side of the wall.')
                end
                if any(abs(xx(1,:))<a')
                    error('Some spheres are overlapping the wall.')
                end
            end

            obj.x = num2cell(xx, 1);
            DurlofskyTensors(obj);
            if obj.setup == 0
                Calculation_gdM(obj);
            elseif obj.setup == 1
                WallTensors(obj)
                Calculation_gdM_Wall(obj)
            end
            
            obj.xG = position.xG;
            obj.X_ID = X_ID;
            obj.Rot_init = position.rot;
        end
        
        function poly_Rc = getRadiusOfCurvatureFit(obj, varargin) 
            % polynomial approximation of the radius of curvature along the
            % radius of the disk using discrete points
            
            if nargin == 1
                % to model Miara's larger disk (2024)
                y = obj.rad;
                points = [0 y/2 y];
                RC = [6.3 5.8 6.9]; 
            else
                points = varargin{1};   % ideally opt for 3 points
                RC = varargin{2};       % radius of curvature values
            end

            poly_Rc = polyfit(points, RC, 2); % quadratic fit

        end

        function xx_bent = applyBending(~, xx, Rc, poly_Rc)
            % to get the U-shaped disk
            
            N = size(xx, 2);
            xx_bent = zeros(3, N);
            
            for i = 1:N
                pos = xx(:,i);
                
                x = pos(1); 
                y = pos(2);
                % z = pos(3);

                if nargin == 4
                    Rc = polyval(poly_Rc, abs(y));
                end
                
                % theta = atan(y/R);
                theta = y/Rc;
                
                % New coordinates
                x_new = x;
                y_new = Rc*sin(theta);
                z_new = Rc*(1 - cos(theta));
                
                pos_bent = [x_new; y_new; z_new];
                
                xx_bent(:,i) = pos_bent;
            end
        end
    end
end