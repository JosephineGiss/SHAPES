classdef Pusher < Swimmer

    properties
        preset_thetas = [0.1   0.05 -0.7 -0.65 -0.55 -0.8  -0.3   -0.2   0.1   0.2;  %1
            -0.3   0.15  0.02 -0.1  -0.7  -0.65 -0.65  -0.65 -0.3   0;    %2
            -0.6   0.1   0.08 -0.05 -0.1  -0.4  -0.6   -0.7  -0.7  -0.1;  %3
            -0.75  0.05  0.05 -0.   -0.05 -0.05 -0.4   -0.5  -0.5  -0.35; %4
            -0.95  0     0.    0.   -0.05  0.08 -0.1   -0.2  -0.3  -0.4;  %5
            -1.03 -0.15  0.05 -0.    0.    0.    0.    -0.05 -0.1  -0.25; %6
            -1.05 -0.25 -0.1   0.1   0.05 -0.   -0.05  -0.    0.   -0.15; %7
            -1.08 -0.28 -0.2   0.    0.23  0.    0.    -0.05 -0.1  -0.18; %8
            -1.10 -0.32 -0.35  0.05  0.26  0.1   0.    -0.0  -0.15 -0.25; %9
            -1.05 -0.45 -0.4  -0.15  0.35  0.23  0.1   -0.05 -0.1  -0.3;  %10
            -1.03 -0.65 -0.45 -0.25  0.4   0.3   0.2   -0.   -0.1 -0.15;  %11
            -1.   -0.5  -0.6  -0.5   0.05  0.4   0.4    0.15 -0.   -0.1;  %12
            -0.4  -0.8  -0.7  -0.6  -0.3   0.2   0.35   0.2   0.2   0.1;  %13
            -0.1  -0.4  -0.8  -0.7  -0.7  -0.2   0.1    0.2   0.3   0.1   %14
            ];
        lamdba = 1.5; % longueur d'onde
    end

    methods 
        
        function obj = Pusher(Nm, alpha, a, l, rho, position, setup)
            arguments
                Nm (1, :) double = 10
                alpha (1, 1) double = 0
                a (1, :) double = [0.0025, 0.0005];
                l (1, :) double = 0.0015 
                rho (1, :) double = [1030, 0]
                position.xB (3, 1) double = [0; 0; 0]
                position.rot (3, 3) double = eye(3)
                position.theta (1, :) double = [0.1   0.05 -0.7 -0.65 -0.55 -0.8  -0.3   -0.2   0.1   0.2]
                setup.setup (1, 1) double = 0
            end
            Nf = 1;
            obj@Swimmer(Nf, Nm, alpha, a, l, rho, xB = position.xB,...
                rot = position.rot, theta = position.theta, setup = setup.setup);
        end
        
        function getSwimmingLaws(obj, option) 
            % première ligne : theta
            % deuxième ligne : theta_dot
            
            arguments 
                obj (1,1) Pusher
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
            elseif option.preset == 2
                for i = 1:obj.Nm
                    obj.theta_laws{1, 1 + i} = @(t) func_SwimmingLaw(obj, i, t, option.T);
                end
            end
            
        end
        function theta_new = func_SwimmingLaw(obj, i, t, T)
            % theta_new = zeros(1, obj.Nm);
            k = 2*pi/obj.lamdba;
            omega = 2*pi/T;
            y = @(x,t)  0.5*(exp(x)-1).*sin(k*x - omega*t); % omega : pulsation
            % y_dot = @(x,t) - 0.5*omega*(exp(x)-1).*cos(k*x - omega*t);
            theta_test = obj.theta(i);
            %%%% attention, ça ne prend en compte que le mouvement plan
            %%%% (x,y) 
            options = optimoptions('fsolve', 'Display', 'none', 'OptimalityTolerance', 10^(-10)); 
            fun = @(theta) y(obj.x{i + 1}(1) + obj.l*cos(theta), t) - (obj.x{i + 1}(2) + obj.l*sin(theta));
            theta_new = fsolve(fun, theta_test, options);
            % fun_dot = @(theta) y_dot(obj.x{i + 1} + obj.l*cos(theta), t);
            % comment calculer theta_dot ?
        end

    end

end