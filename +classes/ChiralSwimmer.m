classdef ChiralSwimmer < classes.ConstrainedAggregate

    properties 
        xB;         % position of the body sphere
        handedness;
        L;          % length of the helix along the z-axis
        rad;        % radius of the helix and the head
        n;          % number of full turns
        beta;          % pitch/(2*pi)

        theta0;     % phase angle                                              
        theta_dot;  % angular velocity of the helix
        theta_alpha; % ref angle of a sphere alpha on the helix 
        z_alpha;    % height of a sphere alpha along the helix
    end

    properties (Access = protected)
        % active forcing and velocities
        Fa; 
        Ta;
        theta_v; 
        theta_omega; 
    end
    
    methods

        function obj = ChiralSwimmer(theta_dot, L, n, rad, a, lambda, rho, shape, position, setup)

            arguments
                theta_dot (1,1) double 
                L (1,1) double = 10
                n (1,1) double {mustBeInteger} = 3 
                rad (1,1) double = 1 % radius of the helix and the big sphere
                a (1,1) double = 0.1 % radius of the small spheres
                lambda (1,1) double = 2.5 % spacing factor
                rho (1,2) double = [1030, 0]
                shape.handedness (1,1) double = 1 % (1 == right-handed, -1 == left-handed)
                position.xB (3,1) double = [0;0;0]
                position.rot (3,3) double = eye(3)
                position.theta0 (1,1) double = 0
                setup.setup (1,1) double = 0
            end
            rot_helix = aux.rot_u([0 0 1], position.theta0);
            
            dAB = lambda*a;
            beta = L/(2*pi*n); 
            c = sqrt(rad^2 + beta^2);
            Lc = 2*pi*n*c; % curvilinear length
            Ns = floor(Lc/dAB + 1); % maximum number of spheres along the helix
            theta = linspace(0, 2*pi*n, Ns);
            rho = [rho(1), rho(2)*ones(1, Ns - 1)];

            % placing spheres
            xx = zeros(3, Ns);
            xx(:, 1) = [0; 0; -(rad + L/2)]; % head position
            for k = 2:Ns % to get center of mass aligned with the z-axis in case of full turns
                xx(:, k) = [rad*cos(theta(k - 1)); shape.handedness*rad*sin(theta(k - 1)); beta*theta(k - 1)];
            end
            % here we use the center of mass of the body as the mass of the tail is zero
            current_xG_helix = classes.ConstrainedAggregate.getCenterOfMass(xx(:, 2:Ns), a, rho(1)); 
            for k = 2:Ns
                xx(:, k) = rot_helix*(xx(:, k) - current_xG_helix);
            end
            
            % readjusting position
            xB = xx(:, 1);
            z_alpha = num2cell(xx(3,:) - xB(3)); 
            for k = 1:Ns
                xx(:, k) = position.xB + position.rot*(xx(:, k) - xB);
            end
            
            % instantiating object
            obj@classes.ConstrainedAggregate(xx, [rad, a*ones(1, Ns - 1)], rho, setup.setup);

            obj.xB = xx(:, 1);
            obj.rad = rad;
            
            obj.handedness = shape.handedness;
            obj.L = L;
            obj.n = n;
            obj.beta = beta;
            obj.z_alpha = z_alpha;
            obj.theta0 = position.theta0;
            obj.theta_dot = theta_dot;
            obj.theta_alpha = num2cell([NaN, theta(1:end - 1)]);
            obj.Rot_init = position.rot;

            obj.theta_v = cell(1, obj.N);
            obj.theta_v{1} = zeros(3, 1);
            obj.theta_omega = cell(1, obj.N);
            obj.theta_omega{1} = zeros(3, 1);
            
        end
        
        function setPosition(obj, xB, rot, t)
            % modify obj properties in order to reposition the swimmer
            % according to its new body position, rotation and swimming 
            % law(t) in the lab frame
            
            % properties modified : xB, x, xG, t_i, theta, theta_dot, 
            % theta_v, theta_omega, gdM
            
            obj.xB = xB; % updated body position
            obj.x{1} = xB;
            t1 = cos(obj.theta0 + obj.theta_dot*t)*rot(:, 1) + sin(obj.theta0 + obj.theta_dot*t)*rot(:, 2);
            t2 = -sin(obj.theta0 + obj.theta_dot*t)*rot(:, 1) + cos(obj.theta0 + obj.theta_dot*t)*rot(:, 2);
            for i = 2:obj.N
                obj.x{i} = xB + obj.z_alpha{i}*rot(:, 3) +...
                    obj.rad*(cos(obj.theta_alpha{i})*t1 + obj.handedness*sin(obj.theta_alpha{i})*t2); % updated positions
            end
            obj.xG = sum([obj.m{:}].*cell2mat(obj.x)./obj.totalMass, 2);
            
            for i = 2:obj.N
                obj.theta_omega{i} = obj.theta_dot*rot(:, 3);
                obj.theta_v{i} = obj.rad*obj.theta_dot*(-obj.handedness*sin(obj.theta_alpha{i})*t1 + cos(obj.theta_alpha{i})*t2);
            end
            
            % we skip the following when considering a system of aggregates
            if ~obj.calledFromIH
                % updated grand resistance matrix
                DurlofskyTensors(obj); 
                if obj.setup == 0 
                    Calculation_gdM(obj); 
                elseif obj.setup == 1
                    WallTensors(obj); 
                    Calculation_gdM_Wall(obj);
                end
            end
            obj.calledFromIH = false;

        end
        function [uB, omegaB] = getVelocities(obj, FText, xB, qB, t, flow) % flagella mass neglected

            arguments % everything given in the lab frame
                obj (1,1) classes.ConstrainedAggregate % N interdependent spheres
                FText (6,1) double % external forces and torques acting on the swimmer
                xB (3,1) double % body position of the swimmer
                qB (4,1) double % quaternions
                t (1,1) double % time (to get the right angles evolution)
                flow.mu (1,1) double = 1 
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)
            end
            % ambient flow (lab frame)
            u_inf = @(x) flow.U0 + cross(flow.Omega_inf, x)...
                + aux.recTransfo(flow.E_inf)*x; 
            
            % Update the positions of the spheres 
            Rot = aux.quaternion_to_matrice(qB); % updated rotation state
            setPosition(obj, xB, Rot, t)
            
            % extra terms for forces and torques arising from swimming
            obj.Fa = zeros(3,1);
            obj.Ta = zeros(3,1);
            [~, ~, ~, ~, ~, R_FU, R_FOmega, R_TU, R_TOmega] = getResistanceTensors(obj, xB);
            small_res_mat = [obj.A, obj.Bt; obj.B, obj.C];
            
            for i = 1:obj.N
                for j = 1:obj.N
                    obj.Fa = obj.Fa + (R_FU{i,j}*obj.theta_v{j} + R_FOmega{i,j}*obj.theta_omega{j});
                    obj.Ta = obj.Ta + (aux.hat_(obj.x{i} - xB)*R_FU{i,j} + R_TU{i,j})*obj.theta_v{j} + (aux.hat_(obj.x{i} - xB)*R_FOmega{i,j} + R_TOmega{i,j})*obj.theta_omega{j};
                end
            end
            FTa = [obj.Fa; obj.Ta];
            FTtotal = FText - FTa;
            ghe = [aux.recTransfo(obj.Gt,"39")*flow.E_inf;...
                    aux.recTransfo(obj.Ht,"39")*flow.E_inf];
            RHS = (1/flow.mu)*FTtotal + ghe; % Right Hand Side
            sol = small_res_mat\RHS; 
            
            % solid velocities in the lab frame
            uB = sol(1:3) + u_inf(xB); 
            omegaB = sol(4:6) + flow.Omega_inf;
        end
        function getFTS(obj, uB, omegaB, flow) 

            arguments % everything given in the lab frame
                obj (1,1) classes.ChiralSwimmer 
                uB (3,1) double = zeros(3,1) % velocity of the body
                omegaB (3,1) double = zeros(3,1) % rotation of the body
                flow.mu (1,1) double = 1 
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)
            end

            % ambient flow (lab frame)
            u_inf = @(x) flow.U0 + cross(flow.Omega_inf, x) + aux.recTransfo(flow.E_inf)*x; 

            gdU = zeros(3*obj.N,1);
            gdOmega = zeros(3*obj.N,1);
            gdE = zeros(9*obj.N,1);

            % spheres velocity
            for i=1:obj.N
                obj.u{i} = uB - aux.hat_(obj.x{i} - obj.xB)*omegaB + obj.theta_v{i};
                obj.omega{i} = omegaB + obj.theta_omega{i};
            end 

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
        end

        function [t_out,y_out] = Dynamics(obj, FText, flow, CI, options)
            
            arguments % everything given in the lab frame
                
                obj (1,1) classes.ConstrainedAggregate 

                FText (6,1) double = zeros(6,1) % buoyancy force in the lab frame (i.e e_i basis)
                flow.mu (1,1) double = 1
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)

                CI.qB (4,1) double = aux.matrice_to_quaternion(obj.Rot_init);
                % DO NOT forget to change the Rot_init property if the
                % initial configuration != reference configuration
                % or assign the corresponding q0
                
                CI.xB (3,1) double = obj.xB
                
                options.tspan double = 0:0.04:100
            end
            
            CI_ = [CI.qB; CI.xB];
            opts = odeset('RelTol',10^(-6),'AbsTol',10^(-6));
            
            [t_out,y_out] = ode45(@(t,y) Stokes_Forces(obj, t, y, FText,...
                flow.mu, flow.U0, flow.Omega_inf, flow.E_inf),...
                options.tspan, CI_, opts);
            
        end
        function y_dot = Stokes_Forces(obj, t, y, FText, mu, U0, Omega_inf, E_inf)
            % differential equation to solve the dynamics
            % y is a column vector
            % y(1:4) quaternion qB
            % y(5:7) coordinates of the body xB
            
            [uB, omegaB] = getVelocities(obj, FText, y(5:7), y(1:4), t,...
                    "mu", mu, "U0", U0, "Omega_inf", Omega_inf, "E_inf", E_inf);

            % outputs
            dq_dt = aux.quaternion_dot(y(1:4), omegaB, "lab"); % new rotation matrix
            dr_dt = uB; % translational velocity in the lab frame
            
            y_dot = [dq_dt; dr_dt];

            disp("t = " + t)
        end
        
        function Animation(obj, t_out, y_out, color, video) 
            arguments
                obj (1,1) classes.ChiralSwimmer
                t_out double
                y_out double
                color double = [0.75, 0.7, 0.65]
                video.option (1,1) string = "no"
                video.filename (1,1) string = "video_cyl"
                video.profile (1,1) string = "MPEG-4"
                video.framerate (1,1) {mustBeNumeric} = 25
                video.orient (1,2) double = [-37.5,30]
                video.layout (1,1) string {mustBeMember(video.layout,["fixed", "mobile"])} = "mobile"
            end

            if all(size(color) == [1,3])
                color = color.*ones(obj.N, 1);
            end
            
            T = length(t_out); % number of time steps
            [X0,Y0,Z0] = sphere(10);
            [end_i, end_j] = size(X0);

            if video.option == "yes"

                if video.filename == ""
                    error("Give a name to the video file. (filename=)")
                end

                filename = video.filename;
                v = VideoWriter(filename, video.profile);
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
                xB_now = y_out(h,5:7)'; % current body position
                rot = aux.quaternion_to_matrice(y_out(h,1:4)); % current rotation
                t = t_out(h);
                
                X_new = cell(1,obj.N);
                Y_new = cell(1,obj.N);
                Z_new = cell(1,obj.N);
                
                X_new_sub = zeros(size(X0));
                Y_new_sub = zeros(size(Y0));
                Z_new_sub = zeros(size(Z0));
                
                xx = zeros(3, obj.N);
                xx(:, 1) = xB_now;
                t1 = cos(obj.theta0 + obj.theta_dot*t)*rot(:, 1) + sin(obj.theta0 + obj.theta_dot*t)*rot(:, 2);
                t2 = -sin(obj.theta0 + obj.theta_dot*t)*rot(:, 1) + cos(obj.theta0 + obj.theta_dot*t)*rot(:, 2);
                for i = 2:obj.N
                    xx(:, i) = xB_now + obj.z_alpha{i}*rot(:, 3) +...
                        obj.rad*(cos(obj.theta_alpha{i})*t1 + obj.handedness*sin(obj.theta_alpha{i})*t2); % updated positions
                end
                
                for k = 1:obj.N
                    X1 = obj.a{k}*X0; 
                    Y1 = obj.a{k}*Y0; 
                    Z1 = obj.a{k}*Z0; 
                    for i = 1:end_i
                        for j = 1:end_j
                            new_grid = rot*([X1(i,j) Y1(i,j) Z1(i,j)]') + xx(:,k);
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
                
                hold off
                ax = gca;
                ax.FontSize = 20;
                axis equal
                xlabel('$x$', Interpreter='latex', FontSize=27)
                ylabel('$y$', Interpreter='latex', FontSize=27)
                zlabel('$z$', Interpreter='latex', FontSize=27)
                
                view(video.orient(1), video.orient(2))
                if video.layout == "mobile"
                    axis_vec = [xB_now(1)-scale xB_now(1)+scale xB_now(2)-scale xB_now(2)+scale xB_now(3)-scale xB_now(3)+scale];
                elseif video.layout == "fixed"
                    axis_vec = [min(y_out(:, 5))-scale max(y_out(:, 5))+scale...
                        min(y_out(:, 6))-scale max(y_out(:, 6))+scale...
                        min(y_out(:, 7))-scale max(y_out(:, 7))+scale];
                end

                if obj.setup == 1
                    hold on
                    lim = axis_vec;
                    if lim(1) > 0
                        lim(1) = 0;
                    elseif lim(2) < 0
                        lim(2) = 0;
                    end
                    [wall_y, wall_z] = meshgrid(linspace(lim(3), lim(4), 21), linspace(lim(5), lim(6), 21)); % Generate y and z data
                    wall_x = zeros(size(wall_y, 1)); % Generate x data
                    surf(wall_x, wall_y, wall_z)
                    hold off
                end

                axis equal
                axis(axis_vec)
                
                set(gcf, 'Color', 'w')
                hold off

                if mod(h, 100) == 0
                    disp(['Frame ' num2str(h) '/' num2str(T)])
                end
                
                if video.option == "yes"
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end

            end

            if video.option == "yes"
                close(v);
            end

        end
        
    end

end