classdef RigidBody < classes.ConstrainedAggregate
    
    properties
        % position of the spheres relative to the center of mass in the
        % reference frame (when Rot_init = eye(3))
        % see definition in the constructor
        X_ID;  
    end

    methods

        function obj = RigidBody(xx, a, rho, setup)
            % xx : positions of spheres, 3*N matrix (ith column for the ith sphere)
            % a : radius of the spheres, either a scalar (monodisperse) or N components vector (polydisperse)
            % rho : mass density of the spheres, either a scalar (homogenous) or N components vector (heterogeneous)
            % setup : 0 or 1

            if nargin == 0
                xx = [[0;0;0], [5;0;0]];         
                a = 1;
                rho = 1030;
                setup = 0;
            elseif nargin == 2
                rho = 1030;
                setup = 0;
            elseif nargin == 3
                setup = 0;
            end
            
            obj@classes.ConstrainedAggregate(xx, a, rho, setup);
            obj.X_ID = mat2cell(cell2mat(obj.x) - obj.xG, 3, ones(1,obj.N)); 
        end
        
        function setPosition(obj, xG, rot, ~)
            
            obj.xG = xG; % updated center of mass
            for i = 1:obj.N
                obj.x{i} = rot*(obj.X_ID{i}) + obj.xG; % updated positions
            end
            
            % we skip the following when considering a system of aggregates
            if ~obj.calledFromIH
                % updated grand resistance matrix
                DurlofskyTensors(obj); 
                if obj.setup == 0 
                    Calculation_gdM(obj); 
                elseif obj.setup == 1
                    WallTensors(obj) 
                    Calculation_gdM_Wall(obj)
                end
            end
            obj.calledFromIH = false;

        end

        function [t_out, y_out] = Dynamics(obj, FText, flow, CI, options) % 1 rigid particle

            arguments % everything given in the lab frame

                obj (1,1) classes.RigidBody
                % Resistance tensors precalculated in the reference frame (Rot_init = eye(3))
                % Sphere positions are evaluated thanks to X_ID and Rot(t)

                FText (:,1) double = zeros(6,1) % buoyancy force in the lab frame (i.e e_i basis)
                flow.mu (1,1) double = 0.01
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)

                % CI.angles (3,1) double = [0; pi/4; 0] % [phi; theta; psi]
                % theta must not be zero (or need to add a new
                % parameterization)
                % angles MUST be specified if compute_angles = "yes"
                % I've included Euler angles into compute_angles

                CI.q0 (4,1) double = aux.matrice_to_quaternion(obj.Rot_init);
                % DO NOT forget to change the Rot_init property if the
                % initial configuration != reference configuration
                % or assign the corresponding q0

                CI.xG0 (3,1) double = obj.xG
                options.tspan double = 0:0.04:100
                options.compute_angles (4,1) double = [0;0;0;0]
                % 1st component is 0 ("no") or 1 ("yes")
                % the three remaining components are Euler angles : [phi; theta; psi]
                % in the initial configuration
                % if compute_angles is set to "yes", q0 is automatically
                % calculated from Euler angles and so does not require to be
                % set (will be rewritten otherwise)
            end

            if options.compute_angles(1) == 1
                % q0 and Euler angles must coincide if options.compute_angles == yes
                CI.q0 = aux.matrice_to_quaternion(aux.rot_euler(options.compute_angles(2),options.compute_angles(3),options.compute_angles(4)));
                CI_ = [CI.q0; CI.xG0; options.compute_angles(2:4)];
            else
                CI_ = [CI.q0; CI.xG0];
            end
            opts = odeset('RelTol',10^(-6),'AbsTol',10^(-6),'OutputFcn', @aux.checkingProgression);

            [t_out, y_out] = ode45(@(t,y) Stokes_Forces(obj, t, y, FText,...
                flow.mu, flow.U0, flow.Omega_inf, flow.E_inf),...
                options.tspan, CI_, opts);
        end
        
        function [uG, omegaG] = getVelocities(obj, FText, xG, quat, ~, flow) % lab frame + non inertial %%% check if calling the method without times works
            
            arguments % everything given in the lab frame
                obj (1, 1) classes.RigidBody % N independent spheres
                FText (6, 1) double  % buoyancy force in the lab frame (i.e e_i basis) 6x1 column vector
                xG (3, 1) double % center of mass of the solid
                quat (4, 1) double % quaternions
                ~ % times : unused for the RigidBody class
                flow.mu (1, 1) double = 1
                flow.U0 (3, 1) double = zeros(3,1)
                flow.Omega_inf (3, 1) double = zeros(3,1)
                flow.E_inf (9, 1) double = zeros(9,1)
            end

            % ambient flow (lab frame)
            u_inf = @(x) flow.U0 + cross(flow.Omega_inf, x)...
                + aux.recTransfo(flow.E_inf)*x; 
            
            Rot = aux.quaternion_to_matrice(quat); % updated rotation state
            
            if obj.setup == 0
                % Here, we use mobile frame for convenience.
                % Resistance tensors are constant over time.
                % Object properties are not updated.

                % You need to place the center of mass of your particle at 
                % the origin and use its default orientation then compute 
                % its resistance tensors.

                % if the center of mass is not at the origin or if the
                % aggregate is not in its default orientation
                if or(any(abs(obj.xG) > 10^(-14)), any(any(obj.Rot_init ~= eye(3))))
                    setPosition(obj, zeros(3,1), eye(3))
                    getResistanceTensors(obj)
                % if resistance tensors are not computed yet
                elseif isempty(obj.xres)
                    getResistanceTensors(obj)
                % if res tensors are not computed from the center of mass
                elseif any(obj.xres ~= obj.xG)
                    getResistanceTensors(obj)
                end
                
                small_res_mat = [obj.A, obj.Bt; obj.B, obj.C];
                FT = [Rot'*FText(1:3); Rot'*FText(4:6)]; % force/torque in the mobile frame (i.e t_i basis)

                GHE = [aux.recTransfo(obj.Gt,"39")*(aux.recTransfo( Rot'*aux.recTransfo(flow.E_inf)*Rot) );...
                        aux.recTransfo(obj.Ht,"39")*(aux.recTransfo( Rot'*aux.recTransfo(flow.E_inf)*Rot) )];
                RHS = (1/flow.mu)*FT + GHE; 
                Sol = small_res_mat\RHS;

                % solid velocities in the mobile frame
                UG = Sol(1:3) + Rot'*u_inf(xG); 
                OmegaG = Sol(4:6) + Rot'*flow.Omega_inf;

                % solid velocities in the lab frame 
                uG = Rot*UG;
                omegaG = Rot*OmegaG;

            elseif obj.setup == 1
                
                setPosition(obj, xG, Rot)
                getResistanceTensors(obj, xG) 
    
                small_res_mat = [obj.A, obj.Bt; obj.B, obj.C];
    
                ghe = [aux.recTransfo(obj.Gt,"39")*flow.E_inf;...
                        aux.recTransfo(obj.Ht,"39")*flow.E_inf];
                RHS = (1/flow.mu)*FText(1:6) + ghe; 
                sol = small_res_mat\RHS; 
    
                % solid velocities in the lab frame
                uG = sol(1:3) + u_inf(xG); 
                omegaG = sol(4:6) + flow.Omega_inf;

            end

        end
        
        function varargout = getFTS(obj, uG , omegaG, flow)

            arguments
                obj (1,1) classes.RigidBody 
                uG (3,1) double
                omegaG (3,1) double
                flow.mu (1,1) double = 1
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)
            end
            
            % ambient flow (lab frame)
            u_inf = @(x) flow.U0 + cross(flow.Omega_inf, x)...
                + aux.recTransfo(flow.E_inf)*x; 
            
            gdU = zeros(3*obj.N,1);
            gdOmega = zeros(3*obj.N,1);
            gdE = zeros(9*obj.N,1);
            
            % spheres velocity
            for i=1:obj.N
                obj.u{i} = uG + cross(omegaG,obj.x{i} - obj.xG);
                obj.omega{i} = omegaG;
            end 
            
            % velocities vector
            for i = 1:obj.N
                gdU(3*i-2:3*i) = obj.u{i} - u_inf(obj.x{i});
                gdOmega(3*i-2:3*i) = obj.omega{i} - flow.Omega_inf;
                gdE(9*i-8:9*i) = zeros(9,1) - flow.E_inf;
            end
            
            res = ((1/flow.mu)*obj.gdM)\[gdU;gdOmega;gdE];
            
            for i = 1:obj.N
                obj.Force{i} = res(3*i - 2:3*i);
                obj.Torque{i} = res(3*obj.N + 3*i - 2:3*obj.N + 3*i);
                obj.Stresslet{i} = res(6*obj.N + 9*i - 8:6*obj.N + 9*i);
            end 
            if nargout > 0
                F = zeros(3,1);
                T = zeros(3,1);
                S = zeros(3,3);
                for i = 1:obj.N
                    F = F + obj.Force{i};
                    T = T + obj.Torque{i} + cross((obj.x{i} - obj.xG),obj.Force{i});
                    S = S + aux.recTransfo(obj.Stresslet{i}) + (1/2)*((obj.x{i} - obj.xG)*obj.Force{i}')...
                                + (1/2)*(obj.Force{i}*(obj.x{i} - obj.xG)') - (1/3)*(obj.x{i} - obj.xG)'*obj.Force{i}*eye(3);
                end
                varargout{1} = F;
                varargout{2} = T;
                varargout{3} = S;
            end
        end
        
        function y_dot = Stokes_Forces(obj, t, y, FText, mu, U0, Omega_inf, E_inf)
            % differential equation to solve the dynamics
            % y is a column vector
            % y(1:4) quaternion q
            % y(5:7) coordinates of the center of mass xG
            % if computed_angles == "yes", y(8:10) Euler angles
            % y(8) = phi, y(9) = theta, y(10) = psi
            
            [uG, omegaG] = getVelocities(obj, FText, y(5:7), y(1:4),...
                "mu", mu, "U0", U0, "Omega_inf", Omega_inf, "E_inf", E_inf); % lab frame
            
            Rot = aux.quaternion_to_matrice(y(1:4));

            % outputs
            dq_dt = aux.quaternion_dot(y(1:4), omegaG, "lab"); % new rotation matrix
            dr_dt = uG; % translational velocity in the lab frame
            OmegaG = Rot'*omegaG;

            % Euler angles option
            if length(y) == 10
                % matrix that characterizes the time evolution of Euler angles
                invEA_mat = [sin(y(10))/sin(y(9)), cos(y(10))/sin(y(9)), 0;...
                    cos(y(10)), -sin(y(10)), 0;...
                    -sin(y(10))*cos(y(9))/sin(y(9)), -cos(y(10))*cos(y(9))/sin(y(9)), 1];
                % output for Euler angles
                dangles_dt = invEA_mat*OmegaG;

                y_dot = [dq_dt; dr_dt; dangles_dt];
                disp("t = " + t)
            elseif length(y) == 7
                y_dot = [dq_dt; dr_dt];
                disp("t = " + t)
            end

        end

        function Animation(obj, t_out, y_out, color, video, settings)
            arguments
                obj (1,1) classes.RigidBody
                t_out double
                y_out double
                color (:, 3) double = [0.75, 0.7, 0.65]
                video.video (1,1) string {mustBeMember(video.video,["yes", "no"])} = "no"
                video.filename (1,1) string = ""
                video.framerate (1,1) {mustBeNumeric} = 25
                video.orient (1,2) double = [-37.5,30]
                video.layout (1,1) string {mustBeMember(video.layout,["fixed", "mobile"])} = "mobile"
                settings.preset string {mustBeMember(settings.preset, ["default", "sketch"])} = "default"
                settings.center_of_mass string {mustBeMember(settings.center_of_mass, ["on", "off"])} = "off"
            end

            if all(size(color) == [1,3])
                color = [0.75, 0.7, 0.65].*ones(obj.N, 1);
            end

            T = length(t_out); % number of time steps
            n_faces = 40;
            [X0,Y0,Z0] = sphere(n_faces);
            [end_i, end_j] = size(X0);
            
            if video.video == "yes"
                
                if video.filename == ""
                    error("Give a name to the video file. (filemane=)")
                end
                
                filename = video.filename;
                v = VideoWriter(filename);
                v.Quality = 95;
                v.FrameRate = video.framerate;
                open(v);
                
            end
            
            xx = cell2mat(obj.x);
            xmax = max(abs(xx(1,:) - obj.xG(1))) + 2*obj.a{1};
            ymax = max(abs(xx(2,:) - obj.xG(2))) + 2*obj.a{1};
            zmax = max(abs(xx(3,:) - obj.xG(3))) + 2*obj.a{1};
            scale = norm([xmax; ymax; zmax]);

            for h = 1:T
                xG_now = y_out(h,5:7)'; % current center of mass
                Rot = aux.quaternion_to_matrice(y_out(h,1:4)); % current rotation
                
                X_new = cell(1, obj.N);
                Y_new = cell(1, obj.N);
                Z_new = cell(1, obj.N);

                X_new_sub = zeros(size(X0));
                Y_new_sub = zeros(size(Y0));
                Z_new_sub = zeros(size(Z0));

                for k = 1:obj.N
                    X1 = obj.a{k}*X0; 
                    Y1 = obj.a{k}*Y0; 
                    Z1 = obj.a{k}*Z0; 
                    for i = 1:end_i
                        for j = 1:end_j
                            new_grid = Rot*([X1(i,j) Y1(i,j) Z1(i,j)]') + xG_now + Rot*obj.X_ID{k};
                            % ^ first rotate the sphere pattern (mesh) then
                            % place at the center of mass 
                            % then place at the right position
                            X_new_sub(i,j) = new_grid(1);
                            Y_new_sub(i,j) = new_grid(2);
                            Z_new_sub(i,j) = new_grid(3);
                        end
                    end
                    X_new{k} = X_new_sub;
                    Y_new{k} = Y_new_sub;
                    Z_new{k} = Z_new_sub;
                end
                
                %%% Plot frames
                % Delete the contents of the axis instead of creating new 
                % objects to lighten MATLAB memories
                cla; 

                for i = 1:obj.N
                    s = surf(X_new{i}, Y_new{i}, Z_new{i});
                    set(s, 'FaceColor', color(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.7)
                    hold on
                end
                
                view(video.orient(1), video.orient(2))
                if video.layout == "mobile"
                    axis_vec = [xG_now(1)-scale xG_now(1)+scale xG_now(2)-scale xG_now(2)+scale xG_now(3)-scale xG_now(3)+scale];
                elseif video.layout == "fixed"
                    axis_vec = [min(y_out(:, 5))-scale max(y_out(:, 5))+scale...
                        min(y_out(:, 6))-scale max(y_out(:, 6))+scale...
                        min(y_out(:, 7))-scale max(y_out(:, 7))+scale];
                end
                
                if obj.setup == 1
                    lim = axis_vec;
                    if lim(1) > 0
                        lim(1) = 0.1;
                    elseif lim(2) < 0
                        lim(2) = 0;
                    end
                    [wall_y, wall_z] = meshgrid(linspace(lim(3), lim(4), n_faces + 1), linspace(lim(5), lim(6), n_faces + 1)); % Generate x and y data
                    wall_x = zeros(size(wall_y, 1)); % Generate z data
                    surf(wall_x, wall_y, wall_z)
                end
                hold off

                axis equal
                axis(axis_vec)

                if settings.preset == "default"
                    camlight headlight
                    lighting gouraud
                    material dull
                    grid on
                    set(gca, 'FontSize', 12, 'LineWidth', 1.2)
                    xlabel('$x$', 'FontSize', 15, 'Interpreter', 'latex')
                    ylabel('$y$', 'FontSize', 15, 'Interpreter', 'latex')
                    zlabel('$z$', 'FontSize', 15, 'Interpreter', 'latex')
                elseif settings.preset == "sketch"
                    camlight headlight
                    lighting gouraud
                    material shiny
                    grid off
                    box off
                    axis off
                end
            
            set(gcf, 'Color', 'w')
            hold off

                if mod(h, 100) == 0
                    disp(['Frame ' num2str(h) '/' num2str(T)])
                end
                
                drawnow
                
                if video.video == "yes"
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end

            end

            if video.video == "yes"
                close(v);
            end

        end

    end

end