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
    end

    properties (Access = protected)
        theta_v; 
        theta_omega; 

        ff7;
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
                position.theta (1, :) double = [0.1   0.05 -0.7 -0.65 -0.55 -0.8  -0.3   -0.2   0.1   0.2]
                setup.setup (1, 1) double = 0
            end
            
            if isscalar(Nm) % same number of spheres per flagella
                Nm = Nm*ones(1, Nf);
                if length(position.theta) == Nm(1)
                    theta = position.theta;
                    position.theta = zeros(1, sum(Nm));
                    idx = cumsum(Nm);
                    for i = 1:Nf
                        position.theta(idx(i) - Nm(1) + 1:idx(i)) = (-1)^(-1 + i)*theta;
                    end
                end
            end
            
            N = sum(Nm) + 1;
            xx = zeros(3,N);
            t_i = cell(1,N);
            fl_idx = [cumsum(Nm) - Nm + 2; cumsum(Nm) + 1];
            
            if length(rho) == 2
                rho = [max(rho), min(rho)*ones(1, N - 1)];
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

            xx(:,1) = zeros(3,1); % coordinates of the body
            
            % body frame
            t1B = [1; 0; 0];
            t2B = [0; 1; 0];
            
            % placing spheres
            for i = 1:Nf
                % vector from the center of the body to the attachment
                % point of a given flagella 
                l0 = cos(alpha(i))*t1B + sin(alpha(i))*t2B; 
                for j = fl_idx(1, i):fl_idx(2, i)
                    t_i{j} = cos(alpha(i) + sum(position.theta(fl_idx(1, i) - 1:j - 1)))*t1B...
                        + sin(alpha(i) + sum(position.theta(fl_idx(1, i) - 1:j - 1)))*t2B;
                    sum_ti = zeros(3,1);
                    for k = fl_idx(1, i):j
                        sum_ti = sum_ti + t_i{k};
                    end
                    xx(:,j) = a(1)*l0 + l*sum_ti;
                end
            end
            
            % readjusting position using the body as the reference
            for k = 1:N
                xx(:,k) = position.xB + position.rot*xx(:,k);
            end
            
            obj@classes.ConstrainedAggregate(xx, a, rho, setup.setup);
            
            obj.xB = position.xB;
            obj.Nf = Nf;
            obj.Nm = Nm;
            obj.l = l;
            obj.alpha = alpha;
            obj.theta = position.theta;
            obj.theta_dot = zeros(1, N - 1);
            obj.fl_idx = fl_idx;

            obj.theta_v = cell(1, N);
            obj.theta_v{1} = zeros(3, 1);
            obj.theta_omega = cell(1, N);
            obj.theta_omega{1} = zeros(3, 1);
            
            % first cell contains the rotation matrix i.e. the expression
            % of the mobile frame basis vectors (column 1 is t1G etc ...)
            t_i{1} = position.rot;
            obj.t_i = t_i;
            
        end
        
        function setPosition(obj, xB, rot, t)
            % modify obj properties in order to repostion the swimmer
            % according to its new x0, rotation and swimming law(t) in the
            % lab frame
            
            % properties modified : xB, x, xG, t_i, theta, theta_dot, 
            % theta_v, theta_omega, gd_M
            
            obj.xB = xB;
            obj.x{1} = xB;
            obj.t_i{1} = rot;
            
            for k = 1:obj.N - 1
                obj.theta(k) = obj.theta_laws{1, k + 1}(t);
                obj.theta_dot(k) = obj.theta_laws{2, k + 1}(t);
            end
            
            theta_om = zeros(1, obj.N);

            for i = 1:obj.Nf
                
                l0 = cos(obj.alpha(i))*rot(:,1) + sin(obj.alpha(i))*rot(:,2);
                theta_om(obj.fl_idx(1, i):obj.fl_idx(2, i)) = cumsum(obj.theta_dot(obj.fl_idx(1, i) - 1:obj.fl_idx(2, i) - 1));
                
                for j = obj.fl_idx(1, i):obj.fl_idx(2, i) % spheres indexes of flagellum i
                    obj.t_i{j} = cos(obj.alpha(i) + sum(obj.theta(obj.fl_idx(1, i) - 1:j - 1)))*rot(:,1)...
                        + sin(obj.alpha(i) + sum(obj.theta(obj.fl_idx(1, i) - 1:j - 1)))*rot(:,2);
                    
                    sum_ts = zeros(3,1);
                    for k = obj.fl_idx(1, i):j
                        sum_ts = sum_ts + obj.t_i{k};
                    end
                    
                    obj.x{j} = obj.xB + obj.a{1}*l0 + obj.l*sum_ts;
                end

            end
            obj.xG = sum([obj.m{:}].*cell2mat(obj.x)./obj.totalMass, 2);

            for i = 1:obj.N
                obj.theta_omega{i} = theta_om(i)*rot(:,3);
            end
            
            for i = 1:obj.Nf
                for j = obj.fl_idx(1, i):obj.fl_idx(2, i) % spheres indexes of flagellum i
                    sum_thetatj = zeros(3,1);
                    for k = obj.fl_idx(1, i):j
                        sum_thetatj = sum_thetatj + cross(obj.theta_omega{k}, obj.t_i{k});
                    end
                    obj.theta_v{j} = obj.l*sum_thetatj;
                end
            end
            
            % updated grand resistance matrix
            DurlofskyTensors(obj); 
            if obj.setup == 0
                Calculation_gdM(obj); 
            elseif obj.setup == 1
                WallTensors(obj);
                Calculation_gdM_Wall(obj); 
            end

        end
        
        function [uB, omegaB] = getVelocities(obj, FText, xB, qB, t, flow) % flagella mass neglected

            arguments % everything given in the lab frame
                obj (1,1) classes.Swimmer % N interdependant spheres
                FText (6,1) double % external forces and torques acting on the swimmer
                xB (3,1) double % body position of the swimmer
                qB (4,1) double % quaternions
                t (1,1) double % time (to get the right angles evolution)
                flow.mu (1,1) double = 1 
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)
            end
            % ambiant flow (lab frame)
            u_inf = @(x) flow.U0 + cross(flow.Omega_inf, x)...
                + aux.recTransfo(flow.E_inf)*x; 
            
            % Update the positions of the spheres 
            Rot = quaternion_to_matrice(qB); % updated rotation state
            setPosition(obj, xB, Rot, t)
            
            % extra terms for forces and torques arising from swimming
            F_theta = zeros(3,1);
            T_theta = zeros(3,1);
            [~, ~, ~, ~, ~, R_FU, R_FOmega, R_TU, R_TOmega] = getResistanceTensors(obj, xB);
            small_res_mat = [obj.A, obj.Bt; obj.B, obj.C];
            
            for i = 1:obj.N
                for j = 1:obj.N
                    F_theta = F_theta + (R_FU{i,j}*obj.theta_v{j} + R_FOmega{i,j}*obj.theta_omega{j});
                    T_theta = T_theta + (hat_(obj.x{i} - xB)*R_FU{i,j} + R_TU{i,j})*obj.theta_v{j} + (hat_(obj.x{i} - xB)*R_FOmega{i,j} + R_TOmega{i,j})*obj.theta_omega{j};
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
            %------------------------------------%
        end
        function [U, V, W, uB, omegaB] = getFlowField(obj, t, FText, xB, qB, X, Y, flow, option) 
            
            arguments % everything given in the lab frame
                obj (1,1) classes.Swimmer % N independant spheres
                t (1,1) double % time (to get the right angles evolution)
                FText (6,1) double  % buoyancy force in the lab frame (i.e e_i basis) 6x1 column vector
                xB (3,1) double = zeros(3,1) % body position of the swimmer
                qB (4,1) double = [1 0 0 0] % quaternions
                X double = NaN
                Y double = NaN
                flow.mu (1,1) double = 1 % in non inertial cases, mu is out of the equations
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)
                option.margin (1,1) double = 1.1
            end
            if isnan(X)
                [X, Y] = meshgrid(-10*obj.a{1} + xB(1) : obj.a{1}/2 : 10*obj.a{1} + xB(1), -10*obj.a{1} + xB(2) : obj.a{1}/2 : 10*obj.a{1} + xB(2));
            end
            
            [uB, omegaB] = getVelocities(obj, FText, xB, qB, t,...
                    "mu", flow.mu, "U0", flow.U0,...
                    "Omega_inf", flow.Omega_inf, "E_inf", flow.E_inf);
            getFTS(obj, uB, omegaB)
            
            [ni, nj] = size(X);
            U = zeros(ni, nj);
            V = zeros(ni, nj);
            W = zeros(ni, nj);
            % if there is a wall, cancel the flow on the other side 
            if obj.setup == 1
                other_side_of_the_wall_condition = gt(obj.x{1}(1),0);
            else
                other_side_of_the_wall_condition = NaN;
            end
            for i = 1:ni
                for j = 1:nj
                    x_c = [X(i,j); Y(i,j); xB(3)];
                    overlap_condition = zeros(1, obj.N);
                    % cancels the flow at the location of the spheres
                    for k = 1:obj.N
                        overlap_condition(k) = norm(x_c - obj.x{k})^2 <= (option.margin*obj.a{k})^2; % marge de 10%
                    end
                    if any(overlap_condition)
                        U(i,j) = NaN;
                        V(i,j) = NaN;
                        W(i,j) = NaN;
                    elseif and(other_side_of_the_wall_condition == true, lt(X(i,j), 0))
                        % if the swimmer is on the "right" side of the wall
                        U(i,j) = NaN;
                        V(i,j) = NaN;
                        W(i,j) = NaN;
                    elseif and(other_side_of_the_wall_condition == false, gt(X(i,j), 0))
                        % if the swimmer is on the "left" side of the wall
                        U(i,j) = NaN;
                        V(i,j) = NaN;
                        W(i,j) = NaN;
                    else
                        [U(i,j), V(i,j), W(i,j)] = FlowField(x_c, obj,...
                            "mu", flow.mu, "U0", flow.U0,...
                            "Omega_inf", flow.Omega_inf, "E_inf", flow.E_inf);
                    end
                end
            end
        end
        function getFTS(obj, uB , omegaB) %%% pb: comment sont déterminées les theta_v et theta_omega sans l'info du temps ? 
            % disturbance flow
            
            % % ambiant flow (lab frame)
            % u_inf = @(x) U0 + cross(Omega_inf, x) + SphereAssembly.RecTransfo(E_inf)*x; % E vector into matrix  % lab frame
            
            gdU = zeros(3*obj.N,1);
            gdOmega = zeros(3*obj.N,1);
            gdE = zeros(9*obj.N,1);
            
            % spheres velocity
            for i=1:obj.N
                obj.u{i} = uB - hat_(obj.x{i})*omegaB + obj.theta_v{i};
                obj.omega{i} = omegaB + obj.theta_omega{i};
            end 

            for i = 1:obj.N
                gdU(3*i-2:3*i) = obj.u{i};
                gdOmega(3*i-2:3*i) = obj.omega{i};
            end
        
            res = obj.gdM\[gdU;gdOmega;gdE];
   
            for i = 1:obj.N
                obj.Force{i} = res(3*i - 2:3*i);
                obj.Torque{i} = res(3*obj.N + 3*i - 2:3*obj.N + 3*i);
                obj.Stresslet{i} = res(6*obj.N + 9*i - 8:6*obj.N + 9*i);
            end 
        end

        function P = PowerOutput(obj, t, xB, qB) % WIP, need FText in input arg + updated getVelocities arg
            if and(nargin == 2, obj.setup == 0)
                xB = obj.xB;
                qB = matrice_to_quaternion(obj.t_i{1});
            end
            [uB, omegaB] = getVelocities(obj, t, xB, qB);
            getFTS(obj, uB, omegaB)
            
            hydro_forces = -cell2mat(obj.Force);
            hydro_torques = -cell2mat(obj.Torque);

            reaction_forces = zeros(3, obj.N);
            % segment_torques = zeros(3, obj.N);
            torque_t1 = zeros(3, obj.N);
            torque_t2 = zeros(3, obj.N);
            torque_t3 = zeros(3, obj.N);
            crosslf = zeros(3, obj.N);
            crosslsumf = zeros(3, obj.N);
            
            for i = 2:obj.N
                crosslf(:, i) = cross(obj.l*obj.t_i{i}, hydro_forces(:, i));
            end
            
            for n = 1:obj.Nf
                start_idx = obj.fl_idx(1, n);
                end_idx = obj.fl_idx(2, n);
                for i = start_idx:end_idx
                    reaction_forces(:, i) = - sum(hydro_forces(:, i:end_idx), 2);
                    torque_t1(:, i) = - sum(hydro_torques(:, i:end_idx), 2);
                    torque_t2(:, i) = - sum(crosslf(:, i:end_idx), 2);
                end
                for i = start_idx:end_idx - 1
                    crosslsumf(:, i) = cross(obj.l*obj.t_i{i}, sum(reaction_forces(:, i + 1:end_idx),2));
                end
            end
            
            % for i = 2:obj.N - 1
            %     crosslsumf(:, i) = cross(obj.t_i{i}, reaction_forces(:, i + 1:end));
            % end

            for n = 1:obj.Nf
                start_idx = obj.fl_idx(1, n);
                end_idx = obj.fl_idx(2, n);
                for i = start_idx:end_idx
                    torque_t3(:, i) = sum(crosslsumf(:, i:end_idx), 2);
                end
            end
            
            segment_torques = torque_t1 + torque_t2 + torque_t3;

            cdottt = zeros(1, obj.N);

            for i = 2:obj.N
                cdottt(i) = dot(segment_torques(:, i), obj.theta_dot(i - 1)*obj.t_i{1}(:,3));
            end
            P = sum(cdottt);
            
        end

        function [t_out,y_out] = Dynamics(obj, FText, flow, CI, options)
            
            arguments % everything given in the lab frame
                
                obj (1,1) classes.Swimmer 
                % Resistance tensors precalculated in the reference frame (Rot_init = eye(3))
                % Sphere positions are evaluated thanks to X_ID and Rot(t)

                FText (:,1) double = zeros(6,1) % buoyancy force in the lab frame (i.e e_i basis)
                flow.mu (1,1) double = 1
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)

                CI.qB (4,1) double = matrice_to_quaternion(obj.Rot_init);
                % DO NOT forget to change the Rot_init property if the
                % inital configuration != reference configuration
                % or assign the corresponding q0
                
                CI.xB (3,1) double = obj.xB
                
                options.tspan double = 0:0.04:100
                % options.compute_angles (4,1) double = [0;0;0;0]
                % % 1st component is 0 ("no") or 1 ("yes")
                % % the three remaining components are Euler angles : [phi; theta; psi]
                % % in the initial configuration
                % % if compute_angles is set to "yes", q0 is automatically
                % % calculted from Euler angles and so does not require to be
                % % set (will be rewritten otherwise)
            end
            
            CI_ = [CI.qB; CI.xB];
            opts = odeset('RelTol',10^(-6),'AbsTol',10^(-6));
            
            % [t_out,y_out] = ode45(@(t,y) Swimmer_Forces(obj, t, y, f, flow.mu, flow.U0, flow.Omega_inf, flow.E_inf), options.tspan, CI_, opts);
            [t_out,y_out] = ode45(@(t,y) Stokes_Forces(obj, t, y, FText,...
                flow.mu, flow.U0, flow.Omega_inf, flow.E_inf),...
                options.tspan, CI_, opts);
            
        end

        function y_dot = Stokes_Forces(obj, t, y, FText, mu, U0, Omega_inf, E_inf)
            % differential equation to solve the dynamics
            % y is a column vector
            % y(1:4) quaternion qB
            % y(5:7) coordinates of the body xB
            % if computed_angles == "yes", y(8:10) Euler angles
            % y(8) = phi, y(9) = theta, y(10) = psi
            
            [uB, omegaB] = getVelocities(obj, FText, y(5:7), y(1:4), t,...
                    "mu", mu, "U0", U0, "Omega_inf", Omega_inf, "E_inf", E_inf);

            % outputs
            dq_dt = quaternion_dot(y(1:4), omegaB, "lab"); % new rotation matrix
            dr_dt = uB; % translational velocity in the lab frame
            
            % OmegaG = Rot'*omegaG;
            %-------------------------------------
            
            y_dot = [dq_dt; dr_dt];

            display(t)
            % display([t, obj.theta(1), obj.theta_dot(1)])
        end
        
        function Animation(obj, t_out, y_out, C, Video) % 
            arguments
                obj (1,1) classes.Swimmer
                t_out double
                y_out double
                C.color double = []
                Video.option (1,1) string = "no"
                Video.filename (1,1) string = "video_cyl"
                Video.format (1,1) string = "mp4"
                Video.profile (1,1) string = "MPEG-4"
                Video.framerate (1,1) {mustBeNumeric} = 25
                Video.orient (1,2) double = [-37.5,30]
                Video.layout (1,1) string {mustBeMember(Video.layout,["fixed", "mobile"])} = "mobile"
            end

            if isempty(C.color)
                C.color = [1 1 1].*ones(obj.N,1);
            end
            
            T = length(t_out); % number of time steps
            [X0,Y0,Z0] = sphere(10);
            [end_i, end_j] = size(X0);

            if Video.option == "yes"
                filename = Video.filename;
                video = VideoWriter(filename, Video.profile);
                video.Quality = 95;
                video.FrameRate = Video.framerate;
                %video.FileFormat = Video.format;
                open(video);
            end

            for h = 1:T
                x0_now = y_out(h,5:7)'; % current body position
                Rot = quaternion_to_matrice(y_out(h,1:4)); % current rotation
                
                % X_new = zeros([obj.N, size(X0)]);
                % Y_new = zeros([obj.N, size(Y0)]);
                % Z_new = zeros([obj.N, size(Z0)]);
                
                X_new = cell(1,obj.N);
                Y_new = cell(1,obj.N);
                Z_new = cell(1,obj.N);
                
                X_new_sub = zeros(size(X0));
                Y_new_sub = zeros(size(Y0));
                Z_new_sub = zeros(size(Z0));
                
                theta_anim = zeros(1, obj.N - 1);
                theta_dot_anim = zeros(1, obj.N - 1);
                for k = 1:obj.N - 1
                    theta_anim(k) = obj.theta_laws{1, k + 1}(t_out(h));
                    theta_dot_anim(k) = obj.theta_laws{2, k + 1}(t_out(h));
                end
                
                t1G = [1; 0; 0];
                t2G = [0; 1; 0];
                ti = cell(1, obj.N);
                xx = zeros(3, obj.N);

                for i = 1:obj.Nf
                
                    l0 = cos(obj.alpha(i))*t1G + sin(obj.alpha(i))*t2G;
                    
                    for j = obj.fl_idx(1, i):obj.fl_idx(2, i)
                        ti{j} = cos(obj.alpha(i) + sum(theta_anim(obj.fl_idx(1, i) - 1:j - 1)))*t1G...
                            + sin(obj.alpha(i) + sum(theta_anim(obj.fl_idx(1, i) - 1:j - 1)))*t2G;
                        sum_ti = zeros(3,1);
                        for k = obj.fl_idx(1, i):j
                            sum_ti = sum_ti + ti{k};
                        end
                        % relative position from the body in the mobile frame
                        xx(:,j) = obj.a{1}*l0 + obj.l*sum_ti;
                    end
    
                end
                
                
                for k = 1:obj.N
                    X1 = obj.a{k}*X0; %+ obj.x{k}(1);
                    Y1 = obj.a{k}*Y0; %+ obj.x{k}(2);
                    Z1 = obj.a{k}*Z0; %+ obj.x{k}(3);
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
                    c = colormap(C.color(i,:));
                    set(s, 'FaceColor', c)
                    % s.EdgeColor = 'none';
                    hold on
                end
                for i = 1:obj.Nf
                    l0 = obj.a{1}*(cos(obj.alpha(i))*t1G + sin(obj.alpha(i))*t2G);
                    l0_lab = Rot*l0 + x0_now + Rot*xx(:,1);
                    markersize = 35; % 100 for mobile, 10 with a wall in the lab
                    scatter3(l0_lab(1), l0_lab(2), l0_lab(3), markersize, 'k', 'filled')
                end
                hold off
                xticks([0 0.02 0.04])
                ax = gca;
                ax.FontSize = 20;
                axis equal
                xlabel('$x$ (mm)', Interpreter='latex', FontSize=27)
                ylabel('$y$ (mm)', Interpreter='latex', FontSize=27)
                zlabel('z')
                
                
                view(Video.orient(1),Video.orient(2))
                % axis ? difficile de penser à une bonne facon d'afficher
                % de manière générale
                if Video.layout == "mobile"
                    scale = 40*obj.a{2};
                    axis_vec = [x0_now(1)-scale x0_now(1)+scale x0_now(2)-scale x0_now(2)+scale x0_now(3)-scale x0_now(3)+scale];
                    axis(axis_vec)
                elseif Video.layout == "fixed"
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
                
                disp(h)
                
                if Video.option == "yes"
                    frame = getframe(gcf);
                    writeVideo(video, frame);
                end

            end

            if Video.option == "yes"
                close(video);
            end

        end

        
        
    end

end