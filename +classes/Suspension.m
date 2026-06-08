classdef Suspension < classes.SphereAssembly
    
    methods 
        
        function obj = Suspension(xx, a, rho, setup)
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

            obj@classes.SphereAssembly(xx, a, rho, setup);
        end

        function [t_out,y_out] = Dynamics(obj, f, flow, CI, options)

            arguments % everything given in the lab frame
                obj (1,1) Suspension % N spheres in free motion (no virtual links)
                f (:,1) double = zeros(6,1) % buoyancy force in the lab frame (i.e e_i basis) 6N column vector
                flow.mu (1,1) double = 0.01
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)
                CI.x0 (:,1) double = cell2mat(obj.x') % 3*N position column vector
                options.tspan double = 0:0.04:200
            end

            for i = 1:obj.N
                % all the spheres have the same orientation at the start 
                q0(4*i - 3:4*i) = matrice_to_quaternion(eye(3)); % 4*N rotation column vector
            end

            CI_ = [q0'; CI.x0];

            opts = odeset('RelTol',10^(-6),'AbsTol',10^(-6));
            [t_out,y_out] = ode45(@(t,y) Stokes_Forces(obj, t, y, f, flow.mu, flow.U0, flow.Omega_inf, flow.E_inf), options.tspan, CI_, opts);
        end

        function [u, omega] = getVelocities(obj, f, x_info, flow) % non-inertial
            % returns the velocities of N free spheres (line vectors)
            % immersed in a linear flow in the LAB frame coordinates
            % for t time steps (t >= 1)
            % for dynamics resolution: t = 1
            % for retroactive resolution of the velocities: t > 1

            arguments % everything given in the lab frame
                obj (1,1) Suspension % N independant spheres
                f (:,1) double  % buoyancy force in the lab frame (i.e e_i basis) 6*N column vector
                x_info double % list with the positions of the N spheres at t time steps (3*N columns)
                flow.mu (1,1) double = 0.01
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)
            end

            [t, ~] = size(x_info); % number of time steps
            u = zeros(t, 3*obj.N);
            omega = zeros(t, 3*obj.N);

            % ambiant flow (lab frame)
            u_inf = @(x) flow.U0 + cross(flow.Omega_inf, x) + SphereAssembly.RecTransfo(flow.E_inf)*x; % E vector into matrix  % lab frame
            E_vec = zeros(9*obj.N,1);
            for i = 1:obj.N
                E_vec(9*i - 8:9*i) = flow.E_inf;
            end

            for time = 1:t % going through time steps

                for i = 1:obj.N
                    obj.x{i} = x_info(3*i - 2 : 3*i)'; % updated centers of mass for each sphere
                end
                % -> obj properties modifications : x
                
                DurlofskyTensors(obj); % -> SA properties modifications : J1, J2, J3, J4_93, R0, R3, R4_93, K1_39, K2_39, K3_39 & K4_99
                if obj.setup == 0
                    Calculation_gd_M(obj); % -> SA properties modifications : gd_M
                elseif obj.setup == 1
                    WallTensors(obj);% -> SA properties modifications : J1star, J2star, J3star, J4star93, R1star, R2star, R3star, R4star93, K1star39, K2star39, K3star39 & K4star99
                    Calculation_gd_M_Wall(obj); % -> SA properties modifications : gd_M
                end

                gdR = inv(obj.gdM);

                R_FU = gdR(1:6*obj.N, 1:6*obj.N);
                R_FE = gdR(1:6*obj.N, 6*obj.N+1:15*obj.N);

                U_vec = (1/flow.mu)*(R_FU\(f + R_FE*E_vec));

                % M_UF = obj.gdM(1:6*obj.N, 1:6*obj.N);
                % M_US = obj.gdM(1:6*obj.N, 6*obj.N+1:15*obj.N);
                % M_EF = obj.gdM(6*obj.N+1:15*obj.N, 1:6*obj.N);
                % M_ES = obj.gdM(6*obj.N+1:15*obj.N, 6*obj.N+1:15*obj.N);
                % 
                % U_vec = (1/flow.mu)*((M_UF - M_US*(M_ES\M_EF))*f - M_US/M_ES*E_vec);

                for i = 1:obj.N
                    u(time, 3*i - 2:3*i) = U_vec(3*i - 2:3*i) + u_inf(obj.x{i});
                    omega(time, 3*i - 2:3*i) = U_vec(3*obj.N + 3*i - 2:3*obj.N + 3*i) + flow.Omega_inf;
                end
            end
        end
        
        function Animation(obj, t_out, y_out, C, Video) % free motion spheres
            arguments
                obj (1,1) Suspension
                t_out double
                y_out double
                C.color double = []
                Video.option (1,1) string = "no"
                Video.filename (1,1) string = "video_cyl"
                Video.format (1,1) string = "mp4"
                Video.framerate (1,1) {mustBeNumeric} = 50
                Video.orient (1,2) double = [-37.5,30]
            end

            if isempty(C.color)
                C.color = [1 1 1].*ones(obj.N,1);
            end

            T = length(t_out); % number of time steps
            [X0,Y0,Z0] = sphere(10);

            [end_i, end_j] = size(X0);

            if Video.option == "yes"
                filename = Video.filename;
                video = VideoWriter(filename);
                video.Quality = 95;
                video.FrameRate = Video.framerate;
                % video.FileFormat = Video.format;
                open(video);
            end

            X_new = cell(1, obj.N);
            Y_new = cell(1, obj.N);
            Z_new = cell(1, obj.N);

            % x_tot = zeros(3,1);

            for h = 1:T
                for n = 1:obj.N
                    X_new{n} = zeros(end_i, end_j);
                    Y_new{n} = zeros(end_i, end_j); 
                    Z_new{n} = zeros(end_i, end_j);

                    xn = y_out(h, 4*obj.N + 3*n - 2 : 4*obj.N + 3*n)'; % position of sphere n
                    Rot = quaternion_to_matrice(y_out(h, 4*n - 3:4*n)); % orientation of sphere n

                    % reshape the sphere radius
                    X1 = obj.a{n}*X0; 
                    Y1 = obj.a{n}*Y0; 
                    Z1 = obj.a{n}*Z0; 
                    for i = 1:end_i
                        for j = 1:end_j
                            new_grid = Rot*([X1(i,j) Y1(i,j) Z1(i,j)]') + xn;
                            % ^ first rotate the sphere pattern (mesh) then
                            % place at the sphere position
                            X_new{n}(i,j) = new_grid(1);
                            Y_new{n}(i,j) = new_grid(2);
                            Z_new{n}(i,j) = new_grid(3);
                        end
                    end
                end

                drawnow
                for i = 1:obj.N
                    s = surf(X_new{i}, Y_new{i}, Z_new{i});
                    c = colormap(C.color(i,:));
                    set(s, 'FaceColor', c)
                    hold on
                end
                hold off
                axis equal
                xlabel('x')
                ylabel('y')
                zlabel('z')

                view(Video.orient(1),Video.orient(2))
                % axis ?

                if obj.setup == 1
                    hold on
                    lim = axis;
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

    methods (Access = protected)

        function y_dot = Stokes_Forces(obj, t, y, f, mu, U0, Omega_inf, E_inf) 
            % differential equation to solve the dynamics
            % y is a column vector
            % 4 components per spheres for the quaternion q
            % then 3 components per spheres for the center of mass r
            % y = (q1; q2; ...; qn; r1; r2; ...; rn);

            dq_dt = zeros(4*obj.N, 1);

            [u, omega] = getVelocities(obj, f, y(4*obj.N + 1:7*obj.N)',...
                "mu", mu, "U0", U0, "Omega_inf", Omega_inf, "E_inf", E_inf);

            for i = 1:obj.N
                dq_dt(4*i - 3 : 4*i, 1) = quaternion_dot(y(4*i - 3 : 4*i), omega(3*i - 2:3*i), "lab"); % new rotation matrix
                dr_dt = u'; % translational velocity in the lab frame
            end
            y_dot = [dq_dt; dr_dt];

            disp("t = " + num2str(t))
        end
        
    end
    
end