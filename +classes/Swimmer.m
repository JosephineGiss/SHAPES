classdef Swimmer < classes.ConstrainedAggregate

    properties 
        xB;         % position of the body sphere
        Nf;         % number de flagellum
        Nm;         % number of spheres per flagella
        l;          % length of a flagella section 
        alpha;      % contact angle of the flagella on the body
        theta;      % relative angle between each segment 
        theta_dot;  % angular velocity of each segment

        theta_laws; % swimming law for each segment
        % first row contains theta evolution
        % second row contains theta_dot evolution

        % mobile frame of each segment + first element is the rotation matrix of the swimmer
        t_i; 

        % 1st row of the ith column is the index of the first sphere that 
        % is part of the ith flagellum, 2nd row is for the last sphere
        fl_idx; 

        mirror; % flagella parameterization can be mirrored
    end

    properties (Access = protected)
        % active forcing and velocities
        Fa;
        Ta;
        theta_v; 
        theta_omega; 
    end
    
    methods (Static, Access = protected)
        
        function [x, t_i] = flagellaGeometry(xB, rot, a1, alpha, theta, lseg, fl_idx)
            Nf = size(fl_idx, 2);
            N = fl_idx(2, end);
            x = zeros(3, N);  
            x(:,1) = xB;
            t_i = cell(1, N);   
            t_i{1} = rot;
            for i = 1:Nf
                idx = fl_idx(1,i):fl_idx(2,i);
                ang = alpha(i) + cumsum(theta(idx - 1)); % absolute angle of each segment
                T = rot(:,1)*cos(ang) + rot(:,2)*sin(ang);   
                l0 = rot(:,1)*cos(alpha(i)) + rot(:,2)*sin(alpha(i));
                x(:,idx) = xB + a1*l0 + cumsum(T*lseg, 2);
                t_i(idx) = num2cell(T, 1);
            end
        end
        
        function theta = normalizeTheta(theta, Nm, mirror)
            Nf = numel(Nm);
            if isempty(theta), theta = 0; end
            assert (iscell(theta) || (isnumeric (theta) && isvector(theta)), ...
                "theta must be a scalar, a vector, or a cell.")
            if iscell(theta)  
                assert(numel(theta) == Nf, "theta cell must contain values for each flagellum.")
                for i = 1:Nf
                    if isempty(theta{i})
                        theta{i} = 0;
                    end
                    assert(isnumeric(theta{i}), "All theta values must be scalar or empty.")
                    if isscalar(theta{i})
                        theta{i} = theta{i}*ones(1, Nm(i));
                    end
                    assert(numel(theta{i}) == Nm(i), "theta{%d} must have %d values", i, Nm(i))
                    theta{i} = theta{i}(:).';
                end
                theta = [theta{:}];
            elseif isscalar(theta)
                theta = theta*ones(1, sum(Nm));
            elseif numel(theta) == sum(Nm)           
                theta = theta(:).';
            elseif all(Nm == numel(theta))        
                theta = cell2mat(arrayfun(@(i) mirror(i)*theta(:).', 1:Nf, 'UniformOutput', false));
            else
                error("Theta length is incompatible with Nm")
            end
        end
        
    end
    
    methods

        function obj = Swimmer(Nf, Nm, alpha, a, l, rho, position, setup)

            arguments
                Nf (1, 1) double = 2
                Nm (1, :) double = 10
                alpha (1, :) double = [3*pi/7, 4*pi/7]
                a (1, :) double = [0.0035, 0.00035];
                l (1, :) double = 0.0010 
                rho (1, :) double = [1030, 0]
                position.xB (3, 1) double = [0; 0; 0]
                position.rot (3, 3) double = eye(3)
                position.mirror (1,:) double = ones(1, Nf)
                position.theta = []
                setup.setup (1, 1) double = 0
            end
            
            if isscalar(Nm) % same number of spheres per flagella
                Nm = Nm*ones(1, Nf);
            end
            assert(numel(Nm) == Nf, "You are reporting more spheres per flagellum than there are flagella.")
            assert(numel(position.mirror) == Nf)
            position.theta = classes.Swimmer.normalizeTheta(position.theta, Nm, position.mirror);
            
            N = sum(Nm) + 1;
            fl_idx = [cumsum(Nm) - Nm + 2; cumsum(Nm) + 1];
            
            if length(rho) == 2
                rho = [rho(1), rho(2)*ones(1, N - 1)];
            end

            if isscalar(position.theta)
                position.theta = position.theta*ones(1, N - 1);
            end

            if isscalar(a) % all the spheres have the same radius
                a = a*ones(1, N);
            elseif length(a) == 2 % body radius then flagellum section radius
                b = a;
                a = zeros(1,N);
                a(1) = max(b); % body radius
                a(2:N) = min(b); % flagellum thickness
            elseif length(a) == Nf + 1 % each flagella has a different section radius
                b = a;
                a = zeros(1, N);
                a(1) = b(1);
                for i = 1:Nf
                    idx = [0, cumsum(Nm)];
                    a(idx(i) + 2:idx(i + 1) + 1) = b(i + 1);
                end
            elseif length(a) ~= N
                error("Input radiuses (a variable) don't fit the number of spheres or flagellum")
            end
            
            if length(alpha) ~= Nf
                error("Please enter all the contact angles (alpha variables)")
            end

            [xx, t_i] = classes.Swimmer.flagellaGeometry(position.xB, position.rot, a(1), alpha, position.theta, l, fl_idx);
            
            obj@classes.ConstrainedAggregate(xx, a, rho, setup.setup);
            
            obj.xB = position.xB;
            obj.Nf = Nf;
            obj.Nm = Nm;
            obj.l = l;
            obj.alpha = alpha;
            obj.theta = position.theta;
            obj.theta_dot = zeros(1, N - 1);
            obj.fl_idx = fl_idx;
            obj.mirror = position.mirror;

            obj.theta_v = cell(1, N);
            obj.theta_v{1} = zeros(3, 1);
            obj.theta_omega = cell(1, N);
            obj.theta_omega{1} = zeros(3, 1);
            
            obj.Rot_init = position.rot;
            obj.t_i = t_i;
            
        end
        
        function setPosition(obj, xB, rot, t)
            % modify obj properties in order to reposition the swimmer
            % according to its new x0, rotation and swimming law(t) in the
            % lab frame
            
            % properties modified : xB, x, xG, t_i, theta, theta_dot, 
            % theta_v, theta_omega, gdM
            
            obj.xB = xB;
            obj.x{1} = xB;
            obj.t_i{1} = rot;
            
            for k = 1:obj.N - 1
                obj.theta(k) = obj.theta_laws{1, k + 1}(t);
                obj.theta_dot(k) = obj.theta_laws{2, k + 1}(t);
            end
            
            theta_om = zeros(1, obj.N);
            
            [x, t_i] = classes.Swimmer.flagellaGeometry(xB, rot, obj.a{1}, obj.alpha, obj.theta, obj.l, obj.fl_idx); 
            obj.x = num2cell(x,1); 
            obj.t_i = t_i;

            for i = 1:obj.Nf
                
                idx = obj.fl_idx(1, i):obj.fl_idx(2, i);
                theta_om(idx) = cumsum(obj.theta_dot(idx - 1));
                obj.theta_omega(idx) = num2cell(rot(:,3)*theta_om(idx), 1);
                T = cell2mat(obj.t_i(idx));
                obj.theta_v(idx) = num2cell(cumsum(obj.l.*cross(rot(:,3).*theta_om(idx), T), 2), 1);

            end
            obj.xG = sum([obj.m{:}].*cell2mat(obj.x)./obj.totalMass, 2);

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
                obj (1,1) classes.Swimmer 
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
            F_theta = zeros(3,1);
            T_theta = zeros(3,1);
            [~, ~, ~, ~, ~, R_FU, R_FOmega, R_TU, R_TOmega] = getResistanceTensors(obj, xB);
            small_res_mat = [obj.A, obj.Bt; obj.B, obj.C];
            
            for i = 1:obj.N
                for j = 1:obj.N
                    F_theta = F_theta + (R_FU{i,j}*obj.theta_v{j} + R_FOmega{i,j}*obj.theta_omega{j});
                    T_theta = T_theta + (aux.hat_(obj.x{i} - xB)*R_FU{i,j} + R_TU{i,j})*obj.theta_v{j} + (aux.hat_(obj.x{i} - xB)*R_FOmega{i,j} + R_TOmega{i,j})*obj.theta_omega{j};
                end
            end
            FTtheta = [F_theta; T_theta];
            FTtotal = FText - FTtheta;
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
                obj (1,1) classes.Swimmer % N independent spheres
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
                
                obj (1,1) classes.Swimmer 

                FText (:,1) double = zeros(6,1) % buoyancy force in the lab frame (i.e e_i basis)
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
                obj (1,1) classes.Swimmer
                t_out double
                y_out double
                color double = [0.75, 0.7, 0.65]
                video.option (1,1) string = "no"
                video.filename (1,1) string = "video_cyl"
                video.format (1,1) string = "mp4"
                video.profile (1,1) string = "MPEG-4"
                video.framerate (1,1) {mustBeNumeric} = 25
                video.orient (1,2) double = [-37.5,30]
                video.layout (1,1) string {mustBeMember(video.layout,["fixed", "mobile"])} = "mobile"
            end

            if all(size(color) == [1,3])
                color = color.*ones(obj.N,1);
            end
            
            T = length(t_out); % number of time steps
            [X0,Y0,Z0] = sphere(20);
            [end_i, end_j] = size(X0);

            if video.option == "yes"
                filename = video.filename;
                v = VideoWriter(filename, video.profile);
                v.Quality = 95;
                v.FrameRate = video.framerate;
                open(v);
            end

            for h = 1:T
                x0_now = y_out(h,5:7)'; % current body position
                Rot = aux.quaternion_to_matrice(y_out(h,1:4)); % current rotation
                
                X_new = cell(1,obj.N);
                Y_new = cell(1,obj.N);
                Z_new = cell(1,obj.N);
                
                X_new_sub = zeros(size(X0));
                Y_new_sub = zeros(size(Y0));
                Z_new_sub = zeros(size(Z0));
                
                theta_anim = zeros(1, obj.N - 1);
                for k = 1:obj.N - 1
                    theta_anim(k) = obj.theta_laws{1, k + 1}(t_out(h));
                end
                
                [xx, ~] = classes.Swimmer.flagellaGeometry(zeros(3,1), eye(3), obj.a{1}, obj.alpha, theta_anim, obj.l, obj.fl_idx);
                
                for k = 1:obj.N
                    X1 = obj.a{k}*X0; 
                    Y1 = obj.a{k}*Y0; 
                    Z1 = obj.a{k}*Z0; 
                    for i = 1:end_i
                        for j = 1:end_j
                            new_grid = Rot*([X1(i,j) Y1(i,j) Z1(i,j)]') + x0_now + Rot*xx(:,k);
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
                
                
                drawnow
                for i = 1:obj.N
                    s = surf(X_new{i}, Y_new{i}, Z_new{i});
                    set(s, 'FaceColor', color(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.7)
                    hold on
                end
                hold off
                ax = gca;
                ax.FontSize = 20;
                axis equal
                xlabel('$x$ (mm)', Interpreter='latex', FontSize=27)
                ylabel('$y$ (mm)', Interpreter='latex', FontSize=27)
                zlabel('$z$ (mm)', Interpreter='latex', FontSize=27)
                
                view(video.orient(1),video.orient(2))
                if video.layout == "mobile"
                    scale = 40*obj.a{2};
                    axis_vec = [x0_now(1)-scale x0_now(1)+scale x0_now(2)-scale x0_now(2)+scale x0_now(3)-scale x0_now(3)+scale];
                    axis(axis_vec)
                elseif video.layout == "fixed"
                    axis_vec = [min(y_out(:, 5))-7*obj.a{1} max(y_out(:, 5))+7*obj.a{1}...
                        min(y_out(:, 6))-4*obj.a{1} max(y_out(:, 6))+4*obj.a{1}...
                        min(y_out(:, 7))-3*obj.a{1} max(y_out(:, 7))+3*obj.a{1}];
                    axis(axis_vec)
                end

                if obj.setup == 1
                    hold on
                    lim = axis_vec;
                    if lim(1) > 0
                        lim(1) = 0;
                    elseif lim(2) < 0
                        lim(2) = 0;
                    end
                    [wall_y, wall_z] = meshgrid(linspace(lim(3), lim(4), 21), linspace(lim(5), lim(6), 21)); % Generate x and y data
                    wall_x = zeros(size(wall_y, 1)); % Generate z data
                    surf(wall_x, wall_y, wall_z)
                    hold off
                end
                
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