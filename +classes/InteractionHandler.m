classdef InteractionHandler < classes.SphereAssembly
    
    properties
        Np;             % number of particles
        particles;      
        idx_particles;
        
        xres;
        % resistance tensors
        A; 
        B;
        Bt;
        C; 
        G;
        H; 
        Gt;
        Ht;
        M;
    end
    
    methods
    
        function system = InteractionHandler(varargin)
            % input arguments are ConstrainedAggregate objects
            % last input argument can be the setup 
            % if no setup specified, setup = 0 (no wall)

            Np = nargin; % number of particles (i.e ConstrainedAggregate objects)
            
            % checking the presence of a setup argument
            if ~isa(varargin{end}, 'classes.ConstrainedAggregate') 
                Np = Np - 1;
                setup = varargin{end};
            else
                setup = 0;
            end
            % checking the minimal required number of particles
            if Np < 2
                error('Minimum two objects are necessary to consider interactions.')
            end
            
            %%% system construction
            particles = cell(1, Np);
            idx = zeros(2, Np);
            N = 0;
            for i = 1:Np
                if ~isa(varargin{i}, 'classes.ConstrainedAggregate')
                    error('Particles must be ConstrainedAggregate objects');
                end
                particles{i} = varargin{i};
                idx(1, i) = N + 1; 
                N = N + particles{i}.N;
                idx(2, i) = N;
            end

            xx = zeros(3, N);
            a = zeros(3, N);
            rho = zeros(3, N);

            for i = 1:Np
                xx(:, idx(1, i):idx(2, i)) = cell2mat(particles{i}.x);
                a(idx(1, i):idx(2, i)) = cell2mat(particles{i}.a);
                rho(idx(1, i):idx(2, i)) = cell2mat(particles{i}.a);
            end

            % system instantiation
            system@classes.SphereAssembly(xx, a, rho, setup);
            system.Np = Np;
            system.particles = particles;
            system.idx_particles = idx;
            
        end

        function getResistanceTensors(system, xr)

            arguments
                system (1, 1) classes.InteractionHandler
                xr (3, :) double = NaN
            end

            if isnan(xr)
                xr = zeros(3, system.Np);
                for i = 1:system.Np
                    % reference point for RigidBody is the center of mass
                    if isa(system.particles{i}, 'classes.RigidBody')
                        xr(:, i) = system.particles{i}.xG;

                    % reference point for Swimmer is the body coordinates
                    elseif isa(system.particles{i}, 'classes.Swimmer')
                        xr(:, i) = system.particles{i}.xB;
                    end
                end
            end
            dim_xr = size(xr);
            if not(dim_xr(2) == system.Np)
                error('Reference points xr must contains one coordinates per particle : size(xr) = [3, system.Np]')
            end

            system.xres = mat2cell(xr, 3, ones(1, system.Np));
            Xr = cell(1, system.N);
            for i = 1:system.Np
                for j = system.idx_particles(1, i):system.idx_particles(2, i) % spheres indexes of particle i
                    Xr{j} = system.x{j} - xr(:, i);
                end
            end

            gdR = inv(system.gdM);
            
            % initialisation
            R_FU = cell(system.N);
            R_FOmega = cell(system.N);
            R_FE = cell(system.N);
            R_TU = cell(system.N);
            R_TOmega = cell(system.N);
            R_TE = cell(system.N);
            R_SU = cell(system.N);
            R_SOmega = cell(system.N);
            R_SE = cell(system.N);
            
            A_mat = cell(system.Np);
            Bt_mat = cell(system.Np);
            B_mat = cell(system.Np);
            C_mat = cell(system.Np);
            G_mat = cell(system.Np);
            Gt_mat = cell(system.Np);
            H_mat = cell(system.Np);
            Ht_mat = cell(system.Np);
            M_mat = cell(system.Np);
            
            for i = 1:system.Np
                for j = 1:system.Np
                    A_mat{i,j} = zeros(3);
                    Bt_mat{i,j} = zeros(3);
                    B_mat{i,j} = zeros(3);
                    C_mat{i,j} = zeros(3);
                    G_mat{i,j} = zeros(3,3,3);
                    Gt_mat{i,j} = zeros(3,3,3);
                    H_mat{i,j} = zeros(3,3,3);
                    Ht_mat{i,j} = zeros(3,3,3);
                    M_mat{i,j} = zeros(3,3,3,3);
                end
            end

            % subdivision of the resistance matrix into blocks
            %%% 1st line %%%
            Block_FU = gdR(1:3*system.N, 1:3*system.N);
            Block_FOmega = gdR(1:3*system.N, 3*system.N + 1:6*system.N);
            Block_FE = gdR(1:3*system.N, 6*system.N + 1:15*system.N);
            %%% 2nd line %%%
            Block_TU = gdR(3*system.N + 1:6*system.N, 1:3*system.N);
            Block_TOmega = gdR(3*system.N + 1:6*system.N, 3*system.N + 1:6*system.N);
            Block_TE = gdR(3*system.N + 1:6*system.N, 6*system.N + 1:15*system.N);
            %%% 3rd line %%%
            Block_SU = gdR(6*system.N + 1:15*system.N, 1:3*system.N);
            Block_SOmega = gdR(6*system.N + 1:15*system.N, 3*system.N + 1:6*system.N);
            Block_SE = gdR(6*system.N + 1:15*system.N, 6*system.N + 1:15*system.N);
            
            % subdivision of blocks into sub-blocks
            %%% square blocks (3x3) %%%
            for ni = 1:3:3*system.N
                for nj = 1:3:3*system.N
                    i_idx = (ni - 1)/3 + 1; 
                    j_idx = (nj - 1)/3 + 1;
                    R_FU{i_idx, j_idx} = Block_FU(ni:ni + 2, nj:nj + 2);
                    R_FOmega{i_idx, j_idx} = Block_FOmega(ni:ni + 2, nj:nj + 2);
                    R_TU{i_idx, j_idx} = Block_TU(ni:ni + 2, nj:nj + 2);
                    R_TOmega{i_idx, j_idx} = Block_TOmega(ni:ni + 2, nj:nj + 2);
                end
            end
            %%% rectangular blocks (3x9) %%%
            for ni = 1:3:3*system.N
                for nj = 1:9:9*system.N
                    i_idx = (ni - 1)/3 + 1;
                    j_idx = (nj - 1)/9 + 1;
                    R_FE{i_idx, j_idx} = aux.recTransfo(Block_FE(ni:ni + 2, nj:nj + 8));
                    R_TE{i_idx, j_idx} = aux.recTransfo(Block_TE(ni:ni + 2, nj:nj + 8));
                end
            end
            %%% rectangular blocks (9x3) %%%
            for ni = 1:9:9*system.N
                for nj = 1:3:3*system.N
                    i_idx = (ni - 1)/9 + 1; 
                    j_idx = (nj - 1)/3 + 1;
                    R_SU{i_idx, j_idx} = aux.recTransfo(Block_SU(ni:ni + 8, nj:nj + 2));
                    R_SOmega{i_idx, j_idx} = aux.recTransfo(Block_SOmega(ni:ni + 8, nj:nj + 2));
                end
            end
            %%% square block (9x9) %%%
            for ni = 1:9:9*system.N
                for nj = 1:9:9*system.N
                    i_idx = (ni - 1)/9 + 1;
                    j_idx = (nj - 1)/9 + 1;
                    R_SE{i_idx, j_idx} = aux.recTransfo(Block_SE(ni:ni + 8, nj:nj + 8));
                end
            end
            
            for i = 1:system.Np
                for j = 1:system.Np
                    % k_alpha count the spheres belonging to the particle i
                    k_alpha = system.idx_particles(1, i) - 1; % counting starts from the last index of the previous particle. 
                    for alpha = 1:system.particles{i}.N
                        % k_beta count the spheres belonging to the particle j
                        k_beta = system.idx_particles(1, j) - 1; % same
                        k_alpha = k_alpha + 1;
                        for beta = 1:system.particles{j}.N
                            k_beta = k_beta + 1;
                            A_mat{i,j} = A_mat{i,j} + R_FU{k_alpha, k_beta}; 
                            Bt_mat{i,j} = Bt_mat{i,j} + (R_FOmega{k_alpha, k_beta} - R_FU{k_alpha, k_beta}*aux.hat_(Xr{k_beta})); 
                            B_mat{i,j} = B_mat{i,j} + R_TU{k_alpha, k_beta} + aux.hat_(Xr{k_alpha})*R_FU{k_alpha, k_beta}; 
                            C_mat{i,j} = C_mat{i,j} + R_TOmega{k_alpha, k_beta} - R_TU{k_alpha, k_beta}*aux.hat_(Xr{k_beta}) ...
                                + aux.hat_(Xr{k_alpha})*R_FOmega{k_alpha, k_beta} - aux.hat_(Xr{k_alpha})*R_FU{k_alpha, k_beta}*aux.hat_(Xr{k_beta});
                            Gt_mat{i,j} = Gt_mat{i,j} + R_FE{k_alpha, k_beta} + tensorprod(R_FU{k_alpha, k_beta},Xr{k_beta});
                            G_mat{i,j} = G_mat{i,j} + R_SU{k_alpha, k_beta} + (1/2)*permute(tensorprod(R_FU{k_alpha, k_beta},Xr{k_alpha}), [3,1,2])...
                                + (1/2)*permute(tensorprod(R_FU{k_alpha, k_beta},Xr{k_alpha}),[1,3,2])...
                                - (1/3)*tensorprod(eye(3), (Xr{k_alpha}'*R_FU{k_alpha, k_beta})'); 
                            H_mat{i,j} = H_mat{i,j} + R_SOmega{k_alpha, k_beta} - tensorprod(R_SU{k_alpha, k_beta}, aux.hat_(Xr{k_beta}), 3, 1)...
                                + (1/2)*(permute(tensorprod(R_FOmega{k_alpha, k_beta}, Xr{k_alpha}), [3,1,2]) ...
                                - permute(tensorprod(R_FU{k_alpha, k_beta}*aux.hat_(Xr{k_beta}), Xr{k_alpha}), [3,1,2]))...
                                + (1/2)*(permute(tensorprod(R_FOmega{k_alpha, k_beta},Xr{k_alpha}),[1,3,2]) ...
                                - permute(tensorprod(R_FU{k_alpha, k_beta}*aux.hat_(Xr{k_beta}), Xr{k_alpha}), [1,3,2]))...
                                - (1/3)*(tensorprod(eye(3), (Xr{k_alpha}'*R_FOmega{k_alpha, k_beta})') ...
                                - tensorprod(eye(3), (Xr{k_alpha}'*R_FU{k_alpha, k_beta}*aux.hat_(Xr{k_beta}))')); 
                            Ht_mat{i,j} = Ht_mat{i,j} + R_TE{k_alpha, k_beta} + tensorprod(R_TU{k_alpha, k_beta}, Xr{k_beta})...
                                + tensorprod(aux.hat_(Xr{k_alpha}),R_FE{k_alpha, k_beta},2,1)...
                                + tensorprod(aux.hat_(Xr{k_alpha}),tensorprod(R_FU{k_alpha, k_beta},Xr{k_beta}),2,1);
                            M_mat{i,j} = M_mat{i,j} + R_SE{k_alpha, k_beta} + tensorprod(R_SU{k_alpha, k_beta}, Xr{k_beta})...
                                    + (1/2)*permute(tensorprod(R_FE{k_alpha, k_beta}, Xr{k_alpha}), [4,1,2,3])...
                                    + (1/2)*permute(tensorprod(tensorprod(R_FU{k_alpha, k_beta}, Xr{k_beta}), Xr{k_alpha}), [4,1,2,3])...
                                    + (1/2)*permute(tensorprod(R_FE{k_alpha, k_beta}, Xr{k_alpha}), [1,4,2,3])...
                                    + (1/2)*permute(tensorprod(tensorprod(R_FU{k_alpha, k_beta}, Xr{k_beta}), Xr{k_alpha}), [1,4,2,3])...
                                    - (1/3)*tensorprod(eye(3), tensorprod(R_FE{k_alpha, k_beta}, Xr{k_alpha}, 1))...
                                    - (1/3)*tensorprod(eye(3), tensorprod(tensorprod(R_FU{k_alpha, k_beta}, Xr{k_beta}), Xr{k_alpha}, 1));
                        end
                    end
                end
            end

            system.A = A_mat;
            system.B = B_mat;
            system.Bt = Bt_mat;
            system.C = C_mat;
            system.Gt = Gt_mat;
            system.Ht = Ht_mat;
            system.G = G_mat;
            system.H = H_mat;
            system.M = M_mat;

        end

        function setPosition(system, xrefs, rots, times) % good for RB, need to check for Swimmer
            % xrefs: cell vector with xref property of each particle inside
            %        the system
            % rots:  cell vector with the rotation matrix of each particle
            % times: cell vector with time info if necessary
            
            % ajdusting input formats (if the input are not in cell format)
            if isa(xrefs, 'double') 
                xrefs = mat2cell(xrefs, 3, ones(1, system.Np));
            end
            if isa(rots, 'double')
                dim_rot = size(rots);
                if dim_rot(1) == 3 % if input are rotation matrices 
                    rots = mat2cell(rots, 3, 3*ones(1, system.Np));
                elseif dim_rot(1) == 4 % if inputs are quaternions
                    rots = mat2cell(rots, 4, ones(1, system.Np));
                end
            end
            % if ~(isa(system.particles, 'classes.RigidBody')) %%%%% to check with the Swimmer
            %     if isa(times, 'double')
            %         times = mat2cell(times, 3, ones(1, system.Np)); 
            %     end
            % end
            
            for i = 1:system.Np
                % if the particle i is a RigidBody, then times{i} won't be
                % called
                system.particles{i}.calledFromIH = true;  
                if isa(system.particles{i}, 'classes.RigidBody')
                    setPosition(system.particles{i}, xrefs{i}, rots{i})
                elseif isa(system.particles{i}, 'classes.Swimmer')
                    setPosition(system.particles{i}, xrefs{i}, rots{i}, times{i})
                end
            end
            
            % updated positions
            for i = 1:system.Np
                idx = system.idx_particles(1, i):system.idx_particles(2, i);
                for j = 1:system.particles{i}.N
                    system.x{idx(j)} = system.particles{i}.x{j};
                end
            end
            
            % updated grand resistance matrix
            DurlofskyTensors(system); 
            if system.setup == 0 
                Calculation_gdM(system); 
            elseif system.setup == 1
                WallTensors(system) 
                Calculation_gdM_Wall(system)
            end
            
        end

        function [urefs, omegarefs] = getVelocities(system, FText, xrefs, quats, times, flow) % WIP, updates for Swimmer needed
            
            arguments % everything given in the lab frame
                system (1, 1) classes.InteractionHandler % N independant spheres
                FText (:, 1) double % external forces [F1; F2; ...] then torques [... T1; T2; ...]
                xrefs (:, 1) double % center of mass of the solid
                quats (:, 1) double % quaternions
                times (1, :) double = NaN % times : unused for the RigidBody class
                flow.mu (1, 1) double = 1
                flow.U0 (3, 1) double = zeros(3,1)
                flow.Omega_inf (3, 1) double = zeros(3,1)
                flow.E_inf (9, 1) double = zeros(9,1)
            end
            
            % ambiant flow (lab frame)
            u_inf = @(x) flow.U0 + cross(flow.Omega_inf, x)...
                + aux.recTransfo(flow.E_inf)*x; 
            E_vec = repmat(flow.E_inf, system.Np, 1);
            
            xrefs_cell = cell(1, system.Np);
            Rot = cell(1, system.Np);
            isRB = false(1, system.Np);
            for i = 1:system.Np
                xrefs_cell{i} = xrefs(3*i - 2 : 3*i); % new positions
                Rot{i} = aux.quaternion_to_matrice(quats(4*i - 3 : 4*i)); % new rotation state
                isRB(i) = isa(system.particles{i}, "classes.RigidBody");
            end
            xrefs = reshape(xrefs,[3,system.Np]);

            % updating positions
            if all(isRB)
                setPosition(system, xrefs_cell, Rot)
            else
                setPosition(system, xrefs_cell, Rot, times)
            end
            getResistanceTensors(system, xrefs) 
            
            A_mat  = zeros(3*system.Np);
            B_mat  = zeros(3*system.Np);
            Bt_mat = zeros(3*system.Np);
            C_mat  = zeros(3*system.Np);
            Gt_mat = zeros(3*system.Np, 9*system.Np);
            Ht_mat = zeros(3*system.Np, 9*system.Np);
            for i = 1:system.Np
                for j = 1:system.Np
                    A_mat(3*i - 2:3*i, 3*j - 2:3*j) = system.A{i,j};
                    B_mat(3*i - 2:3*i, 3*j - 2:3*j) = system.B{i,j};
                    Bt_mat(3*i - 2:3*i, 3*j - 2:3*j) = system.Bt{i,j};
                    C_mat(3*i - 2:3*i, 3*j - 2:3*j) = system.C{i,j};
                    Gt_mat(3*i - 2:3*i, 9*j - 8:9*j) = aux.recTransfo(system.Gt{i,j}, "39");
                    Ht_mat(3*i - 2:3*i, 9*j - 8:9*j) = aux.recTransfo(system.Ht{i,j}, "39");
                end
            end
            small_res_mat = [A_mat, Bt_mat; B_mat, C_mat];
            
            %%%%%%%
            %%% ici pas de termes de forçage liés à la nage par ex
            %%% voir comment faire dans le cas de plusieurs nageurs ?
            %%% pour un nageur et un solide indéformable ?
            ghe = [Gt_mat*E_vec; Ht_mat*E_vec];
            RHS = (1/flow.mu)*FText + ghe; 
            sol = small_res_mat\RHS; 

            % particles velocities in the lab frame
            urefs = cell(1, system.Np);
            omegarefs = cell(1, system.Np);
            for i = 1:system.Np
                urefs{i} = sol(3*i - 2:3*i) + u_inf(xrefs_cell{i}); 
                omegarefs{i} = sol(3*system.Np + 3*i - 2: 3*system.Np + 3*i) + flow.Omega_inf;
            end

        end

        function y_dot = Stokes_Forces(system, t, y, FText, mu, U0, Omega_inf, E_inf)
            % everything is calculated in the lab frame
            % 4 components per particles for the quaternion q
            % then 3 components per particles for the reference point x
            % y = (q1; q2; ...; qn; x1; x2; ...; xn);
            
            [urefs, omegarefs] = getVelocities(system, FText, y(4*system.Np + 1:7*system.Np), ...
                y(1:4*system.Np), t, "mu", mu, "U0", U0, "Omega_inf", Omega_inf, "E_inf", E_inf);
            
            Rot = cell(1, system.Np);
            for i = 1:system.Np
                Rot{i} = aux.quaternion_to_matrice(y(4*i - 3 : 4*i)); % new rotation state
            end

            % outputs
            dq_dt = zeros(4*system.Np, 1);
            dr_dt = zeros(3*system.Np, 1);
            for i = 1:system.Np
                dq_dt(4*i - 3:4*i) = aux.quaternion_dot(y(4*i - 3 : 4*i), omegarefs{i}, "lab"); % new rotation matrix
                dr_dt(3*i - 2:3*i) = urefs{i}; % translational velocity in the lab frame
            end
            
            y_dot = [dq_dt; dr_dt];
            % display(t)
        end

        function [t_out, y_out] = Dynamics(system, FText, flow, CI, options) % WIP (Swimmer)

            arguments % everything given in the lab frame

                system (1,1) classes.InteractionHandler
                % Resistance tensors precalculated in the reference frame (Rot_init = eye(3))
                % Sphere positions are evaluated thanks to X_ID and Rot(t)

                FText (:,1) double = NaN % buoyancy force in the lab frame (i.e e_i basis)
                flow.mu (1,1) double = 1
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)

                % CI.angles (3,1) double = [0; pi/4; 0] % [phi; theta; psi]
                % theta must not be zero (or need to add a new
                % parametrisation)
                % angles MUST be specified if compute_angles = "yes"
                % I've included Euler angles into compute_angles

                CI.q0 (:,1) double = NaN
                % DO NOT forget to change the Rot_init property if the
                % inital configuration != reference configuration
                % or assign the corresponding q0

                CI.xref0 (:,1) double = NaN
                options.tspan double = 0:0.04:100
            end
            
            if isnan(FText)
                FText = zeros(6*system.Np, 1);
            end
            if isnan(CI.q0)
                CI.q0 = zeros(4*system.Np, 1);
                for i = 1:system.Np
                    CI.q0(4*i - 3:4*i) = aux.matrice_to_quaternion(system.particles{i}.Rot_init);
                end
            end
            if isnan(CI.xref0)
                CI.xref0 = zeros(3*system.Np, 1);
                for i = 1:system.Np
                    if isa(system.particles{i}, "classes.RigidBody")
                        CI.xref0(3*i - 2:3*i) = system.particles{i}.xG;
                    elseif isa(system.particles{i}, "classes.Swimmer")
                        CI.xref0(3*i - 2:3*i) = system.particles{i}.xB;
                    end
                end
            end
            
            CI_ = [CI.q0; CI.xref0];
            opts = odeset('RelTol',10^(-6),'AbsTol',10^(-6),'OutputFcn', @aux.checkingProgression);
            
            [t_out, y_out] = ode45(@(t,y) Stokes_Forces(system, t, y, FText,...
                flow.mu, flow.U0, flow.Omega_inf, flow.E_inf),...
                options.tspan, CI_, opts);
        end

        function Animation(system, t_out, y_out, color, video, settings) % WIP (Swimmer)
            arguments
                system (1,1) classes.InteractionHandler
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
                color = [0.75, 0.7, 0.65].*ones(system.Np, 1);
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
            
            xx = cell2mat(system.x);
            xscale = max(max(xx(1,:)) - min(xx(1,:)));
            yscale = max(max(xx(2,:)) - min(xx(2,:)));
            zscale = max(max(xx(3,:)) - min(xx(3,:)));
            scale = norm([xscale; yscale; zscale]);
            
            xref_now = cell(1, system.Np);
            Rot = cell(1, system.Np);
            N = zeros(1, system.Np);

            for i = 1:system.Np
                N(i) = system.particles{i}.N;
                max_x(i) = max(y_out(:, 4*system.Np + 3*i - 2)); 
                max_y(i) = max(y_out(:, 4*system.Np + 3*i - 1));
                max_z(i) = max(y_out(:, 4*system.Np + 3*i));
                min_x(i) = min(y_out(:, 4*system.Np + 3*i - 2)); 
                min_y(i) = min(y_out(:, 4*system.Np + 3*i - 1));
                min_z(i) = min(y_out(:, 4*system.Np + 3*i));
            end
            Nmax = max(N);
            X_new = cell(Nmax, system.Np);
            Y_new = cell(Nmax, system.Np);
            Z_new = cell(Nmax, system.Np);
            
            for h = 1:T
                for p = 1:system.Np
                    xref_now{p} = y_out(h, 4*system.Np + 3*p - 2:4*system.Np + 3*p)'; % current ref positions
                    Rot{p} = aux.quaternion_to_matrice(y_out(h, 4*p - 3:4*p)); % current rotation
                    
                    X_new_sub = zeros(size(X0));
                    Y_new_sub = zeros(size(Y0));
                    Z_new_sub = zeros(size(Z0));

                    for k = 1:system.particles{p}.N
                        X1 = system.particles{p}.a{k}*X0; 
                        Y1 = system.particles{p}.a{k}*Y0; 
                        Z1 = system.particles{p}.a{k}*Z0; 
                        for i = 1:end_i
                            for j = 1:end_j
                                new_grid = Rot{p}*([X1(i,j) Y1(i,j) Z1(i,j)]') + xref_now{p} + Rot{p}*system.particles{p}.X_ID{k};
                                % ^ first rotate the sphere pattern (mesh) then
                                % place at the center of mass 
                                % then place at the right position
                                X_new_sub(i,j) = new_grid(1);
                                Y_new_sub(i,j) = new_grid(2);
                                Z_new_sub(i,j) = new_grid(3);
                            end
                        end
                        X_new{k, p} = X_new_sub;
                        Y_new{k, p} = Y_new_sub;
                        Z_new{k, p} = Z_new_sub;
                    end
                end
                
                %%% Plot frames
                % Delete the contents of the axis instead of creating new 
                % objects to lighten MATLAB memories
                cla; 
                
                xref_vec = zeros(3, system.Np);
                for p = 1:system.Np
                    for k = 1:system.particles{p}.N
                        s = surf(X_new{k, p}, Y_new{k, p}, Z_new{k, p});
                        set(s, 'FaceColor', color(p, :), 'EdgeColor', 'none', 'FaceAlpha', 0.7)
                        hold on
                    end
                    xref_vec(:, p) = xref_now{p};
                end
                
                view(video.orient(1), video.orient(2))
                if video.layout == "mobile"
                    axis_vec = [mean(xref_vec(1, :))-scale mean(xref_vec(1, :))+scale...
                        mean(xref_vec(2, :))-scale mean(xref_vec(2, :))+scale...
                        mean(xref_vec(3, :))-scale mean(xref_vec(3, :))+scale];
                elseif video.layout == "fixed"
                    axis_vec = [min(min_x)-scale max(max_x)+scale...
                        min(min_y)-scale max(max_y)+scale...
                        min(min_z)-scale max(max_z)+scale];
                end
                

                if system.setup == 1
                    lim = axis_vec;
                    if lim(1) > 0
                        lim(1) = 0;
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