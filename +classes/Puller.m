classdef Puller < classes.Swimmer

    properties (Constant)
        % properties preallocated to reproduce the swim motion of the green
        % algae Chlamydomonas reinhardtii
        % using Nm = 10, and alpha = [3*pi/7, 4*pi/7].
        preset_thetas = ...
            [0.1   0.05 -0.7  -0.65 -0.55 -0.8  -0.3   -0.2   0.1   0.2;  
            -0.3   0.15  0.02 -0.1  -0.7  -0.65 -0.65  -0.65 -0.3   0;    
            -0.6   0.1   0.08 -0.05 -0.1  -0.4  -0.6   -0.7  -0.7  -0.1;  
            -0.75  0.05  0.05 -0.   -0.05 -0.05 -0.4   -0.5  -0.5  -0.35; 
            -0.95  0     0.    0.   -0.05  0.08 -0.1   -0.2  -0.3  -0.4;  
            -1.03 -0.15  0.05 -0.    0.    0.    0.    -0.05 -0.1  -0.25; 
            -1.05 -0.25 -0.1   0.1   0.05 -0.   -0.05  -0.    0.   -0.15; 
            -1.08 -0.28 -0.2   0.    0.23  0.    0.    -0.05 -0.1  -0.18; 
            -1.10 -0.32 -0.35  0.05  0.26  0.1   0.    -0.0  -0.15 -0.25; 
            -1.05 -0.45 -0.4  -0.15  0.35  0.23  0.1   -0.05 -0.1  -0.3;  
            -1.03 -0.65 -0.45 -0.25  0.4   0.3   0.2   -0.   -0.1 -0.15;  
            -1.   -0.5  -0.6  -0.5   0.05  0.4   0.4    0.15 -0.   -0.1;  
            -0.4  -0.8  -0.7  -0.6  -0.3   0.2   0.35   0.2   0.2   0.1;  
            -0.1  -0.4  -0.8  -0.7  -0.7  -0.2   0.1    0.2   0.3   0.1   
            ];
    end

    methods 
        
        function obj = Puller(Nm, alpha, a, l, rho, position, setup)
            arguments
                Nm (1, :) double = 10
                alpha (1, 2) double = [3*pi/7, 4*pi/7]
                a (1, :) double = [0.0035, 0.00035];
                l (1, :) double = 0.0010 
                rho (1, :) double = [1030, 0]
                position.xB (3, 1) double = [0; 0; 0]
                position.rot (3, 3) double = eye(3)
                position.mirror (1,:) double = [1 -1]
                position.theta (1, :) double = classes.Puller.preset_thetas(1,:)
                setup.setup (1, 1) double = 0
            end
            Nf = numel(alpha);
            obj@classes.Swimmer(Nf, Nm, alpha, a, l, rho, xB = position.xB,...
                rot = position.rot, mirror = position.mirror,...
                theta = position.theta, setup = setup.setup);
        end
        
        function getSwimmingLaws(obj, option)

            arguments
                obj (1,1) classes.Puller
                option.T (1,:) double = 0.5
                option.data cell = {}
                option.phase (1,:) double = 0
                option.nHarm (1,1) double = Inf
            end

            Nf = obj.Nf;
            T = option.T.* ones(1, Nf);
            phase = option.phase.* ones(1, Nf);
            data = option.data;
            if isempty(data)
                assert(all(obj.Nm == 10), "Default law requires 10 spheres per flagellum.")
                data = repmat({classes.Puller.preset_thetas}, 1, Nf);
            end
            for i = 1:Nf
                D = obj.mirror(i)*data{i};
                [nt, nm] = size(D);
                assert(nm == obj.Nm(i), "data{%d} must have %d columns", i, obj.Nm(i))
                w = 2*pi/T(i);
                C = fft(D, [], 1)/nt;
                % nh = min(option.nHarm, floor((nt - 1)/2));
                nh = min(option.nHarm, floor(nt/2));
                k = (1:nh)';
                a0 = real(C(1,:));
                ak = 2*real(C(2:nh+1,:));
                bk = -2*imag(C(2:nh+1,:));

                if mod(nt,2) == 0 && nh == nt/2 % cosine only
                    ak(end,:) = real(C(nh+1,:));
                    bk(end,:) = 0;
                end

                for j = 1:nm
                    idx = obj.fl_idx(1,i) + j - 1;
                    obj.theta_laws{1,idx} = @(t) a0(j) + ak(:,j).'*cos(k*w*(t - phase(i))) ...
                        + bk(:,j).'*sin(k*w*(t - phase(i)));
                    obj.theta_laws{2,idx} = @(t) w*( -(k.*ak(:,j)).'*sin(k*w*(t - phase(i))) ...
                        + (k.*bk(:,j)).'*cos(k*w*(t - phase(i))) );
                end
                
            end

        end

    end

end