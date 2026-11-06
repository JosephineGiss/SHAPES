classdef (Abstract) ConstrainedAggregate < classes.SphereAssembly
    
    properties
        xG; % position of the centre of mass

        xres; % point from which resistance tensors are computed
        %%% Resistance tensors;
        A; 
        B;
        Bt;
        C; 
        G;
        H; 
        Gt;
        Ht;
        M;

        xmob; % point from which mobility tensors are computed
        %%% Mobility tensors
        a_mob;
        b_mob;
        c_mob;

        Ug; % velocity of the center of mass
        Omega_g; % angular velocity of the solid (does it work for swimmer ?)

        Rot_init = eye(3); % rotation matrix relative to the center of mass at t = 0
    end

    properties (Access = protected)
        totalMass; % (only used to get xG so useless for Suspension)
    end

    properties (Hidden, Access = public)
        calledFromIH = false;
    end
    
    methods (Static)

        function current_xG = getCenterOfMass(xx, a, rho)
            N = size(xx, 2);
            if isscalar(rho)
                rho = rho.*ones(N,1);
            end
            if isscalar(a)
                a = a.*ones(N,1);
            end
            if iscolumn(rho)
                rho = rho';
            end
            if iscolumn(a)
                a = a';
            end

            m = rho.*(4/3).*pi.*a.^3;
            total_mass = sum(m);
            if total_mass == 0
                current_xG = [0; 0; 0];
            else
                current_xG = sum(m.*xx./total_mass, 2);
            end
        end

    end
    
    methods (Abstract)
        
        setPosition(obj, xref, rot, t)
        % Defines the new position of the aggregate then update its grand
        % resistance matrix 
        % Inputs:
        %   xref - Reference position [3x1] 
        %          RigidBody: center of mass
        %          Swimmer: body position
        %   rot  - Rotation matrix [3x3]
        %   t    - Time (scalar, if necessary)

        getVelocities(obj, FText, xref, quat, time, flow)

        Dynamics(obj, FText, flow, CI, options) 
        Stokes_Forces(obj, t, y, FText, mu, U0, Omega_inf, E_inf)

        getFTS(obj, u_ref , omega_ref, flow)
        % Update the force, torque and stresslet exerted by each sphere on
        % the fluid (Force, Torque and Stresslet properties)
    end

    methods
        
        function obj = ConstrainedAggregate(xx, a, rho, setup)
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
            obj.totalMass = sum([obj.m{:}]);
            if obj.totalMass == 0
                obj.xG = [0; 0; 0];
            else
                obj.xG = sum([obj.m{:}].*cell2mat(obj.x)./obj.totalMass, 2);
            end
        end

        function varargout = getResistanceTensors(obj, xr)

            arguments
                obj (1, 1) classes.ConstrainedAggregate
                xr (3, 1) double = NaN
            end

            if isnan(xr)
                xr = obj.xG;
            end

            obj.xres = xr;
            Xr = cell(1,obj.N);
            for i = 1:obj.N
                Xr{i} = obj.x{i} - xr;
            end
            
            gdR = inv(obj.gdM);
            
            % initialization
            R_FU = cell(obj.N);
            R_FOmega = cell(obj.N);
            R_FE = cell(obj.N);
            R_TU = cell(obj.N);
            R_TOmega = cell(obj.N);
            R_TE = cell(obj.N);
            R_SU = cell(obj.N);
            R_SOmega = cell(obj.N);
            R_SE = cell(obj.N);

            A_mat = zeros(3);
            Bt_mat = zeros(3);
            B_mat = zeros(3);
            C_mat = zeros(3);
            G_mat = zeros(3,3,3);
            Gt_mat = zeros(3,3,3);
            H_mat = zeros(3,3,3);
            Ht_mat = zeros(3,3,3);
            M_mat = zeros(3,3,3,3);

            % subdivision of the resistance matrix into blocks
            %%% 1st line %%%
            Block_FU = gdR(1:3*obj.N, 1:3*obj.N);
            Block_FOmega = gdR(1:3*obj.N, 3*obj.N + 1:6*obj.N);
            Block_FE = gdR(1:3*obj.N, 6*obj.N + 1:15*obj.N);
            %%% 2nd line %%%
            Block_TU = gdR(3*obj.N + 1:6*obj.N, 1:3*obj.N);
            Block_TOmega = gdR(3*obj.N + 1:6*obj.N, 3*obj.N + 1:6*obj.N);
            Block_TE = gdR(3*obj.N + 1:6*obj.N, 6*obj.N + 1:15*obj.N);
            %%% 3rd line %%%
            Block_SU = gdR(6*obj.N + 1:15*obj.N, 1:3*obj.N);
            Block_SOmega = gdR(6*obj.N + 1:15*obj.N, 3*obj.N + 1:6*obj.N);
            Block_SE = gdR(6*obj.N + 1:15*obj.N, 6*obj.N + 1:15*obj.N);

            % subdivision of blocks into sub-blocks
            %%% square blocks (3x3) %%%
            for ni = 1:3:3*obj.N
                for nj = 1:3:3*obj.N
                    i_idx = (ni - 1)/3 + 1; 
                    j_idx = (nj - 1)/3 + 1;
                    R_FU{i_idx, j_idx} = Block_FU(ni:ni + 2, nj:nj + 2);
                    R_FOmega{i_idx, j_idx} = Block_FOmega(ni:ni + 2, nj:nj + 2);
                    R_TU{i_idx, j_idx} = Block_TU(ni:ni + 2, nj:nj + 2);
                    R_TOmega{i_idx, j_idx} = Block_TOmega(ni:ni + 2, nj:nj + 2);
                end
            end
            %%% rectangular blocks (3x9) %%%
            for ni = 1:3:3*obj.N
                for nj = 1:9:9*obj.N
                    i_idx = (ni - 1)/3 + 1;
                    j_idx = (nj - 1)/9 + 1;
                    R_FE{i_idx, j_idx} = aux.recTransfo(Block_FE(ni:ni + 2, nj:nj + 8));
                    R_TE{i_idx, j_idx} = aux.recTransfo(Block_TE(ni:ni + 2, nj:nj + 8));
                end
            end
            %%% rectangular blocks (9x3) %%%
            for ni = 1:9:9*obj.N
                for nj = 1:3:3*obj.N
                    i_idx = (ni - 1)/9 + 1; 
                    j_idx = (nj - 1)/3 + 1;
                    R_SU{i_idx, j_idx} = aux.recTransfo(Block_SU(ni:ni + 8, nj:nj + 2));
                    R_SOmega{i_idx, j_idx} = aux.recTransfo(Block_SOmega(ni:ni + 8, nj:nj + 2));
                end
            end
            %%% square block (9x9) %%%
            for ni = 1:9:9*obj.N
                for nj = 1:9:9*obj.N
                    i_idx = (ni - 1)/9 + 1;
                    j_idx = (nj - 1)/9 + 1;
                    R_SE{i_idx, j_idx} = aux.recTransfo(Block_SE(ni:ni + 8, nj:nj + 8));
                end
            end
            
            for alpha = 1:obj.N
                for beta = 1:obj.N
                    A_mat = A_mat + R_FU{alpha,beta};
                    Bt_mat = Bt_mat + (R_FOmega{alpha,beta} - R_FU{alpha,beta}*aux.hat_(Xr{beta}));
                    B_mat = B_mat + R_TU{alpha,beta} + aux.hat_(Xr{alpha})*R_FU{alpha,beta};
                    C_mat = C_mat + R_TOmega{alpha,beta} - R_TU{alpha,beta}*aux.hat_(Xr{beta}) + aux.hat_(Xr{alpha})*R_FOmega{alpha,beta}...
                        - aux.hat_(Xr{alpha})*R_FU{alpha,beta}*aux.hat_(Xr{beta});
                    Gt_mat = Gt_mat + R_FE{alpha, beta} + tensorprod(R_FU{alpha, beta},Xr{beta});
                    G_mat = G_mat + R_SU{alpha, beta} + (1/2)*permute(tensorprod(R_FU{alpha,beta},Xr{alpha}),[3,1,2])...
                        + (1/2)*permute(tensorprod(R_FU{alpha,beta},Xr{alpha}),[1,3,2])...
                        - (1/3)*tensorprod(eye(3),(Xr{alpha}'*R_FU{alpha,beta})'); 
                    Ht_mat = Ht_mat + R_TE{alpha, beta} + tensorprod(R_TU{alpha, beta},Xr{beta}) + tensorprod(aux.hat_(Xr{alpha}),R_FE{alpha,beta},2,1)...
                        + tensorprod(aux.hat_(Xr{alpha}),tensorprod(R_FU{alpha, beta},Xr{beta}),2,1);
                    H_mat = H_mat + R_SOmega{alpha, beta} - tensorprod(R_SU{alpha, beta}, aux.hat_(Xr{beta}),3,1)...
                        + (1/2)*(permute(tensorprod(R_FOmega{alpha,beta},Xr{alpha}),[3,1,2]) ...
                        - permute(tensorprod(R_FU{alpha,beta}*aux.hat_(Xr{beta}),Xr{alpha}),[3,1,2]))...
                        + (1/2)*(permute(tensorprod(R_FOmega{alpha,beta},Xr{alpha}),[1,3,2]) ...
                        - permute(tensorprod(R_FU{alpha,beta}*aux.hat_(Xr{beta}),Xr{alpha}),[1,3,2]))...
                        - (1/3)*(tensorprod(eye(3),(Xr{alpha}'*R_FOmega{alpha,beta})') ...
                        - tensorprod(eye(3),(Xr{alpha}'*R_FU{alpha,beta}*aux.hat_(Xr{beta}))'));  
                    M_mat = M_mat + R_SE{alpha, beta} + tensorprod(R_SU{alpha, beta}, Xr{beta})...
                        + (1/2)*permute(tensorprod(R_FE{alpha, beta}, Xr{alpha}),[4,1,2,3])...
                        + (1/2)*permute(tensorprod(tensorprod(R_FU{alpha, beta}, Xr{beta}), Xr{alpha}),[4,1,2,3])...
                        + (1/2)*permute(tensorprod(R_FE{alpha, beta}, Xr{alpha}),[1,4,2,3])...
                        + (1/2)*permute(tensorprod(tensorprod(R_FU{alpha, beta}, Xr{beta}), Xr{alpha}),[1,4,2,3])...
                        - (1/3)*tensorprod(eye(3),tensorprod(R_FE{alpha, beta}, Xr{alpha}, 1))...
                        - (1/3)*tensorprod(eye(3), tensorprod(tensorprod(R_FU{alpha, beta}, Xr{beta}), Xr{alpha}, 1));
                end
            end

            Gt_sym = 0.5*(Gt_mat + permute(Gt_mat, [1 3 2]));
            Ht_sym = 0.5*(Ht_mat + permute(Ht_mat, [1 3 2]));
            M_sym = 0.5*(M_mat + permute(M_mat, [1 2 4 3]));
            trace_Gt = zeros(3,1);
            trace_Ht = zeros(3,1);
            trace_M = zeros(3,3);
            for i = 1:3
                trace_M = trace_M + M_mat(:,:,i,i);
                for l = 1:3
                    trace_Gt(i) = trace_Gt(i) + Gt_mat(i,l,l);
                    trace_Ht(i) = trace_Ht(i) + Ht_mat(i,l,l);
                end
            end
            delta = eye(3);
            traceless_correction_Gt = (1/3)*trace_Gt.*reshape(delta, [1,3,3]);
            traceless_correction_Ht = (1/3)*trace_Ht.*reshape(delta, [1,3,3]);
            traceless_correction_M = (1/3)*trace_M.*reshape(delta, [1,1,3,3]);
            Gt_final = Gt_sym - traceless_correction_Gt;
            Ht_final = Ht_sym - traceless_correction_Ht;
            M_final = M_sym - traceless_correction_M;

            obj.A = A_mat;
            obj.B = B_mat;
            obj.Bt = Bt_mat;
            obj.C = C_mat;
            obj.Gt = Gt_final;
            obj.Ht = Ht_final;
            obj.G = G_mat;
            obj.H = H_mat;
            obj.M = M_final;

            if nargout > 0
                varargout{1} = obj.A;
                varargout{2} = obj.B;
                varargout{3} = obj.C;
                varargout{4} = obj.G;
                varargout{5} = obj.H;
                varargout{6} = R_FU;
                varargout{7} = R_FOmega;
                varargout{8} = R_TU;
                varargout{9} = R_TOmega;
            end
            
        end
        function varargout = getMobilityTensors(obj)
            obj.a_mob = inv(obj.A - obj.Bt*(obj.C\obj.B)); % 5.22 Kim & Karrila
            obj.b_mob = - (obj.C\obj.B)/(obj.A - obj.Bt*(obj.C\obj.B)); % 5.23 Kim & Karrila
            obj.c_mob = inv(obj.C - obj.B*(obj.A\obj.Bt)); % 5.24 Kim & Karrila
            obj.xmob = obj.xres;

            if nargout > 0
                varargout{1} = obj.a_mob;
                varargout{2} = obj.b_mob;
                varargout{3} = obj.c_mob;
            end
        end
        function [xr, varargout] = CenterOfResistance(obj)
            % from Happel & Brenner 5-4.15
            % by default, the function uses the center of mass as its
            % reference

            % Levi-Civita symbol
            epsilon = zeros(3,3,3);
            epsilon(1,2,3) = 1;
            epsilon(2,3,1) = 1;
            epsilon(3,1,2) = 1;
            epsilon(3,2,1) = -1;
            epsilon(2,1,3) = -1;
            epsilon(1,3,2) = -1;

            %%%
            if isempty(obj.xres)
                getResistanceTensors(obj);
            elseif any(not(obj.xres == obj.xG))
                getResistanceTensors(obj);
            end
            if nargout > 1
                varargout{1} = obj.A;
                varargout{2} = obj.B;
                varargout{3} = obj.C;
                varargout{4} = obj.G;
                varargout{5} = obj.H;
            end
            r_or = -inv(trace(obj.A)*eye(3) - obj.A)*tensorprod(epsilon,obj.B,[2 3], [1 2]);
            xr = obj.xG + r_or;
        end
        function [xm, varargout] = CenterOfMobility(obj)
            % by default, the function uses the center of mass as its
            % reference

            % Levi-Civita symbol
            epsilon = zeros(3,3,3);
            epsilon(1,2,3) = 1;
            epsilon(2,3,1) = 1;
            epsilon(3,1,2) = 1;
            epsilon(3,2,1) = -1;
            epsilon(2,1,3) = -1;
            epsilon(1,3,2) = -1;

            %%%
            if isempty(obj.xres)
                getResistanceTensors(obj);
            elseif any(not(obj.xres == obj.xG))
                getResistanceTensors(obj);
            end
            if isempty(obj.xmob)
                getMobilityTensors(obj);
            elseif any(not(obj.xmob == obj.xG))
                getMobilityTensors(obj);
            end
            if nargout > 1 
                varargout{1} = obj.a_mob;
                varargout{2} = obj.b_mob;
                varargout{3} = obj.c_mob;
            end

            r_om = -inv(trace(obj.c_mob)*eye(3) - obj.c_mob)*tensorprod(epsilon, obj.b_mob', [2 3], [1 2]);
            xm = obj.xG + r_om;
        end
        function varargout = translationTheorem_RT(obj, x_new)
            % from Happel & Brenner 1983 
            % and from Kim & Karrila 1991 (5.2.3)

            % x_new : point from which the function is going to compute the
            % resistance tensors
            
            x_ref = obj.xres; % point from which resistance tensors are currently computed
            vec_x = x_new - x_ref;
            
            obj.C = obj.C - aux.hat_(vec_x)*obj.A*aux.hat_(vec_x) + obj.B*aux.hat_(vec_x) - aux.hat_(vec_x)*obj.Bt; % (5-4.10 HB)
            obj.B = obj.B - aux.hat_(vec_x)*obj.A; % (5-4.12 HB)
            obj.Bt = obj.B';
            
            obj.xres = x_new;
            
            if nargout > 0
                varargout{1} = obj.A;
                varargout{2} = obj.B;
                varargout{3} = obj.C;
            end
            
        end
        function varargout = translationTheorem_MT(obj, x_new)
            % from Kim & Karrila 1991 (5.3.1)
            
            % x_new : point from which the function is going to compute the
            % mobility tensors
            
            x_ref = obj.xmob; % point from which mobility tensors are currently computed
            vec_x = x_new - x_ref;
            
            obj.a_mob = obj.a_mob - aux.hat_(vec_x)*obj.c_mob*aux.hat_(vec_x) - aux.hat_(vec_x)*obj.b_mob + obj.b_mob'*aux.hat_(vec_x);
            obj.b_mob = obj.b_mob + obj.c_mob*aux.hat_(vec_x);
            
            obj.xmob = x_new;
            
            if nargout > 0
                varargout{1} = obj.a_mob;
                varargout{2} = obj.b_mob;
                varargout{3} = obj.c_mob;
            end
        end
        
        function [U, V, W, uB, omegaB] = getFlowField(obj, t, FText, xref, qref, X, Y, flow, option)
            % Velocity field on a plane in space (laboratory coordinate system).
            % Available planes (option.plane): "xy", "xz", "yz" 
            
            arguments % everything given in the lab frame
                obj (1,1) classes.ConstrainedAggregate % N independent spheres
                t (1,1) double % time (to get the right angles evolution)
                FText (6,1) double % buoyancy force in the lab frame (i.e e_i basis)
                xref (3,1) double = zeros(3,1) % reference position of the object
                qref (4,1) double = [1 0 0 0] % reference quaternions
                X double = NaN
                Y double = NaN
                flow.mu (1,1) double = 1
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)
                option.margin (1,1) double = 1.1
                option.plane (1,1) string {mustBeMember(option.plane, ["xy","xz","yz"])} = "xy"
                option.height (1,1) double = 0
            end
            
            if isscalar(X) && isnan(X)
                error("Choose a grid [X,Y].")
            end
            
            % Velocities and forces/torques
            [uB, omegaB] = getVelocities(obj, FText, xref, qref, t, ...
                "mu", flow.mu, "U0", flow.U0, ...
                "Omega_inf", flow.Omega_inf, "E_inf", flow.E_inf);
            getFTS(obj, uB, omegaB, ...
                "mu", flow.mu, "U0", flow.U0, ...
                "Omega_inf", flow.Omega_inf, "E_inf", flow.E_inf)
            
            [ni, nj] = size(X);
            U = zeros(ni, nj);
            V = zeros(ni, nj);
            W = zeros(ni, nj);
            
            % If there is a wall, cancel the flow on the other side 
            hasWall = (obj.setup == 1);
            if hasWall
                sphereRight = obj.x{1}(1) > 0;
            end
            
            for i = 1:ni
                for j = 1:nj
                    if option.plane == "xy"
                        x_c = [X(i,j); Y(i,j); option.height];
                    elseif option.plane == "xz"
                        x_c = [X(i,j); option.height; Y(i,j)];
                    elseif option.plane == "yz"
                        x_c = [option.height; X(i,j); Y(i,j)];
                    end
                    % cancels the flow at the location of the spheres
                    overlap = false;
                    for k = 1:obj.N
                        if norm(x_c - obj.x{k})^2 <= (option.margin*obj.a{k})^2
                            overlap = true;
                            break
                        end
                    end
                    
                    wrongSide = false;
                    if hasWall
                        if sphereRight && x_c(1) < 0
                            wrongSide = true;
                        elseif ~sphereRight && x_c(1) > 0
                            wrongSide = true;
                        end
                    end
                    
                    if overlap || wrongSide
                        U(i,j) = NaN; V(i,j) = NaN; W(i,j) = NaN;
                    else
                        [U(i,j), V(i,j), W(i,j)] = FlowField(x_c, obj, ...
                            "mu", flow.mu, "U0", flow.U0, ...
                            "Omega_inf", flow.Omega_inf, "E_inf", flow.E_inf);
                    end
                end
            end
        end
        
        function Show(obj, color, settings)
            
            arguments % everything given in the lab frame
                obj (1,1) classes.ConstrainedAggregate % N interdependent spheres
                color (:, 3) double = [0.75, 0.7, 0.65]
                settings.preset string {mustBeMember(settings.preset, ["default", "sketch"])} = "default"
                settings.center_of_mass string {mustBeMember(settings.center_of_mass, ["on", "off"])} = "off"
            end

            if nargin == 1
                color = [0.75, 0.7, 0.65].*ones(obj.N,1);
            end

            Show@classes.SphereAssembly(obj, color, preset = settings.preset)
            
            hold on
            if settings.center_of_mass == "on"
                plot3(obj.xG(1), obj.xG(2), obj.xG(3), 'kx', 'MarkerSize', 12, 'LineWidth', 3)
            end
            hold off

        end

    end

end