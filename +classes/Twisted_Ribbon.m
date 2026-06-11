classdef Twisted_Ribbon < classes.RigidBody

    properties
        chirality; 
        L; % total length
        s; % distance along the axis of the helix for it to rotate by 2 pi
        W; % width
        thetaf; % twist angle
        lambda;
    end
    
    methods
        
        function obj = Twisted_Ribbon(chirality, a, L, W, s, lambda, rho, position, setup)
            
            arguments
                chirality (1,1) double {mustBeMember(chirality, [-1, 1])} = 1
                a (1,1) double = 0.00025
                L (1,1) double = 0.0075
                W (1,1) double = 0.004
                s (1,1) double = 0.02
                lambda (1,1) double = 2.5
                rho (1,1) double = 1030
                position.xG (3,1) double = [0;0;0] % mobile frame origin
                position.rot (3,3) double = eye(3)
                setup.setup (1,1) double = 0
            end
            
            dAB = lambda*a;
            thetaf = (L*2*pi)/s;
            
            % dAB must to be choose here as a multiple of the length
            Nb_rows = floor(L/dAB) + 1; % numbers of rows along the length
            new_lambda = L/((Nb_rows - 1)*a);
            if Nb_rows == 1
                new_lambda = lambda;
            end
            real_dAB = new_lambda*a;
            
            N_row = floor(W/real_dAB + 1); % number of spheres along the row
            N = N_row*Nb_rows;
            xx = zeros(3, N);
            
            x = 0;
            y = linspace(-(W/2), W/2, N_row);
            
            for i = 1:Nb_rows
                L_row = (i - 1)*real_dAB; % position of the row along the length
                k = 0;
                for j = (i - 1)*N_row + 1:i*N_row % count the spheres on the row
                    k = k + 1;
                    xx(:, j) = aux.rot_u([0;0;1], chirality*(2*pi*((i-1)*real_dAB - L/2)/s))*[x; y(k); L_row - (L/2)];
                end
            end
            
            % readjusting position
            current_xG = classes.ConstrainedAggregate.getCenterOfMass(xx, a, rho);
            X_ID = mat2cell(xx - current_xG, 3, ones(1, N)); 
            for k = 1:N
                xx(: ,k) = position.rot*(xx(:, k) - current_xG) + position.xG;
            end
            
            % object initialization
            obj@classes.RigidBody(xx, a, rho, setup.setup); 
            obj.L = L;
            obj.s = s; 
            obj.W = W;
            obj.thetaf = thetaf; 
            obj.lambda = new_lambda;
            obj.X_ID = X_ID;
            obj.Rot_init = position.rot;
            
        end
        
    end

end