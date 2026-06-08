classdef Puller < classes.Swimmer

    properties
        % properties preallocated to reproduce the swim motion of the green
        % algae Chlamydomonas reinhardtii
        preset_alphas = [3*pi/7, 4*pi/7];
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
                position.theta (1, :) double = [0.1   0.05 -0.7 -0.65 -0.55 -0.8  -0.3   -0.2   0.1   0.2]
                setup.setup (1, 1) double = 0
            end
            Nf = 2;
            obj@classes.Swimmer(Nf, Nm, alpha, a, l, rho, xB = position.xB,...
                rot = position.rot, theta = position.theta, setup = setup.setup);
        end
        
        function getSwimmingLaws(obj, option) 
            % première ligne : theta
            % deuxième ligne : theta_dot
            
            arguments 
                obj (1,1) classes.Puller
                option.T (1,1) double = 2 % period (irl T = 0.02s)
                option.preset = 1
            end
            
            if option.preset == 1
                % Fourier fit of the theta angles based on the experimental
                % paper Geyer et. al 2013
                if obj.Nm ~= 10
                    error('The swimmer must have 2 flaggela of 10 spheres.')
                end
                thetas = obj.preset_thetas;
                f = cell(1,obj.Nm(1));
                for i = 1:obj.Nm(1)
                    yfit = [thetas(:,i); thetas(:,i); thetas(1,i)];
                    xfit = linspace(0, 2*option.T, length(yfit));
                    [f7,~] = fit(xfit',yfit,"fourier7");
                    f{i} = f7;
                end
                
                for i = 1:obj.Nf
                    for j = 1:obj.Nm(i)
                        coeffs_and_freq = coeffvalues(f{j});
                        coeffs = (-1)^(1 + i)*coeffs_and_freq(1:end - 1);
                        w = coeffs_and_freq(end);
                        obj.theta_laws{1, obj.fl_idx(1, i) + j - 1} = @(t) coeffs(1) + coeffs(2)*cos(w*t)...
                            + coeffs(3)*sin(w*t) + coeffs(4)*cos(2*w*t) + coeffs(5)*sin(2*w*t) + coeffs(6)*cos(3*w*t) + coeffs(7)*sin(3*w*t)...
                            + coeffs(8)*cos(4*w*t) + coeffs(9)*sin(4*w*t) + coeffs(10)*cos(5*w*t) + coeffs(11)*sin(5*w*t)...
                            + coeffs(12)*cos(6*w*t) + coeffs(13)*sin(6*w*t) + coeffs(14)*cos(7*w*t) + coeffs(15)*sin(7*w*t);
                        obj.theta_laws{2, obj.fl_idx(1, i) + j - 1} = @(t) -w*coeffs(2)*sin(w*t)...
                            + w*coeffs(3)*cos(w*t) - 2*w*coeffs(4)*sin(2*w*t) + 2*w*coeffs(5)*cos(2*w*t) - 3*w*coeffs(6)*sin(3*w*t) + 3*w*coeffs(7)*cos(3*w*t)...
                            - 4*w*coeffs(8)*sin(4*w*t) + 4*w*coeffs(9)*cos(4*w*t) - 5*w*coeffs(10)*sin(5*w*t) + 5*w*coeffs(11)*cos(5*w*t)...
                            - 6*w*coeffs(12)*sin(6*w*t) + 6*w*coeffs(13)*cos(6*w*t) - 7*w*coeffs(14)*sin(7*w*t) + 7*w*coeffs(15)*cos(7*w*t);
                    end
                end
                obj.ff7 = f;
            end
            
        end

    end

end