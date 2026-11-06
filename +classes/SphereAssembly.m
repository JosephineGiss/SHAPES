classdef SphereAssembly < handle
    
    
    properties (Access = public)
        N; % number of spheres
        x; % positions of the  spheres
        a;
        rho; % density of the spheres;
        m; % mass of the spheres;
        
        gdM; % grand mobility matrix
        
        u; % velocities of the spheres
        omega; % angular velocities of the spheres
        Force;
        Torque;
        Stresslet;
        
        setup; % determines the environment & boundaries of the SphereAssembly
        % O (default) : SphereAssembly in an unbounded fluid
        % 1 : SphereAssembly is near a plane wall placed at x = 0
    end 
    
    properties (Access = protected)
        
        % Durlofsky tensors
        J1;
        J2;
        J3;
        J4_93;
        R0;
        R3;
        R4_93;
        K1_39;
        K2_39;
        K3_39;
        K4_99;
        
        % wall tensors
        J1star;
        J2star;
        J3star;
        J4star93;
        R1star;
        R2star;
        R3star;
        R4star93;
        K1star39;
        K2star39;
        K3star39;
        K4star99;

    end 

    methods (Static)
        
        %%% Durlofsky tensors %%%
        function [res] = tens_J0(x) 
            r=norm(x);
            res=((1/r)*eye(3)+(1/r^3)*(x*x'));
        end
        function [res] = tens_J1(x,beta) 
            r=norm(x);
            res=((1/r+beta^2/(3*r^3))*eye(3)+(1/r^3-beta^2/r^5)*(x*x'));
        end
        function [res] = tens_J2(x,alpha,beta) 
            r=norm(x);
            res=((1/r+(alpha^2+beta^2)/(3*r^3))*eye(3)+(1/r^3-(alpha^2+beta^2)/r^5)*(x*x'));
        end
        function [res] = tens_J3(x) 
            r=norm(x);
            res=(1/r^3)*[[0,x(3),-x(2)];[-x(3),0,x(1)];[x(2),-x(1),0]];
        end       
        function [res] = tens_J4_93(x1,alpha,beta)
            r=norm(x1);
            x=x1(1); y=x1(2); z=x1(3);
            r3=r^3; r5=r^5; r7=r^7;
            x2=x^2; y2=y^2; z2=z^2;
            x3=x^3; y3=y^3; z3=z^3;
            xzy=x*z*y; 
            xy2=x*y2; 
            xz2=x*z2;
            yx2=y*x2;
            yz2=y*z2;
            zy2=z*y2;
            zx2=z*x2;
            beta2=beta^2; alpha2=alpha^2;
            alpha2_r5=alpha2/r5; beta2_r5=beta2/r5; 
            alpha2_r7=alpha2/r7; beta2_r7=beta2/r7;
            
            J4_121 = -3*yx2/r5 - y*beta2_r5 - 3/5*y*alpha2_r5 + 3*yx2*alpha2_r7 + 5*yx2*beta2_r7;
            J4_122 = -3*xy2/r5 - x*beta2_r5 - 3/5*x*alpha2_r5 + 3*xy2*alpha2_r7 + 5*xy2*beta2_r7;
            J4_123 = 3*xzy*alpha2_r7 - 3*xzy/r5 + 5*xzy*beta2_r7;

            J4_131 = -3*zx2/r5 - z*beta2_r5 - 3/5*z*alpha2_r5 + 3*zx2*alpha2_r7 + 5*zx2*beta2_r7;
            J4_133 = -3*xz2/r5 - x*beta2_r5 - 3/5*x*alpha2_r5 + 3*xz2*alpha2_r7 + 5*xz2*beta2_r7;

            J4_232 = -3*zy2/r5 - z*beta2_r5 - 3/5*z*alpha2_r5 + 3*zy2*alpha2_r7 + 5*zy2*beta2_r7;
            J4_233 = -3*yz2/r5 - y*beta2_r5 - 3/5*y*alpha2_r5 + 3*yz2*alpha2_r7 + 5*yz2*beta2_r7;
            
            res=[[x/r3 - 3*x3/r5 - 3*x*beta2_r5 - 9/5*x*alpha2_r5 + 3*x3*alpha2_r7 + 5*x3*beta2_r7,...
                y/r3 - 3*yx2/r5 - y*beta2_r5 - 3/5*y*alpha2_r5 + 3*yx2*alpha2_r7 + 5*yx2*beta2_r7,...
                z/r3 - 3*zx2/r5 - z*beta2_r5 - 3/5*z*alpha2_r5 + 3*zx2*alpha2_r7 + 5*zx2*beta2_r7];...
                [x/r3 - 3*xy2/r5 - x*beta2_r5 - 3/5*x*alpha2_r5 + 3*xy2*alpha2_r7 + 5*xy2*beta2_r7,...
                y/r3 - 3*y3/r5 - 3*y*beta2_r5 - 9/5*y*alpha2_r5 + 3*y3*alpha2_r7 + 5*y3*beta2_r7,...
                z/r3 - 3*zy2/r5 - z*beta2_r5 - 3/5*z*alpha2_r5 + 3*zy2*alpha2_r7 + 5*zy2*beta2_r7];...
                [x/r3 - 3*xz2/r5 - x*beta2_r5 - 3/5*x*alpha2_r5 + 3*xz2*alpha2_r7 + 5*xz2*beta2_r7,...
                y/r3 - 3*yz2/r5 - y*beta2_r5 - 3/5*y*alpha2_r5 + 3*yz2*alpha2_r7 + 5*yz2*beta2_r7,...
                z/r3 - 3*z3/r5 - 3*z*beta2_r5 - 9/5*z*alpha2_r5 + 3*z3*alpha2_r7 + 5*z3*beta2_r7];...
                [J4_121,...
                J4_122,...
                J4_123];...
                [J4_131,...
                J4_123,...
                J4_133];...
                [J4_123,...
                J4_232,...
                J4_233];...
                [J4_121,...
                J4_122,...
                J4_123];...
                [J4_131,...
                J4_123,...
                J4_133];...
                [J4_123,...
                J4_232,...
                J4_233]];
        end 
       
        function [res] = tens_R0(x) 
            r=norm(x);
            res=(1/r^3)*[[0,x(3),-x(2)];[-x(3),0,x(1)];[x(2),-x(1),0]];
        end      
        function [res] = tens_R3(x1)
            r=norm(x1);
            x=x1(1); y=x1(2); z=x1(3);

            res= [[-1/2*1/r^3 + 3/2*x^2/r^5, 3/2*x*y/r^5, 3/2*x*z/r^5];...
                [3/2*x*y/r^5, -1/2*1/r^3 + 3/2*y^2/r^5, 3/2*y*z/r^5];...
                [3/2*x*z/r^5, 3/2*y*z/r^5, -1/2*1/r^3 + 3/2*z^2/r^5]];
        end          
        function [res] = tens_R4_93(x1)
            r=norm(x1);
            x=x1(1); y=x1(2); z=x1(3);

            res= [[0, -3*x*z/r^5, 3*x*y/r^5];...
                [3*y*z/r^5, 0, -3*x*y/r^5];...
                [-3*y*z/r^5, 3*x*z/r^5, 0];...
                [3/2*x*z/r^5, -3/2*y*z/r^5, -3/2*x^2/r^5 + 3/2*y^2/r^5];...
                [-3/2*x*y/r^5, 3/2*x^2/r^5 - 3/2*z^2/r^5, 3/2*y*z/r^5];...
                [-3/2*y^2/r^5 + 3/2*z^2/r^5, 3/2*x*y/r^5, -3/2*x*z/r^5];...
                [3/2*x*z/r^5, -3/2*y*z/r^5, -3/2*x^2/r^5 + 3/2*y^2/r^5];...
                [-3/2*x*y/r^5, 3/2*x^2/r^5 - 3/2*z^2/r^5, 3/2*y*z/r^5];...
                [-3/2*y^2/r^5 + 3/2*z^2/r^5, 3/2*x*y/r^5, -3/2*x*z/r^5]];
        end
 
        function [res] = tens_K1_39(x1,beta)
            r=norm(x1);
            x=x1(1); y=x1(2); z=x1(3);

            res=[[-9/5*x*beta^2/r^5 + 3*x^3*beta^2/r^7 + x/r^3 - 3*x^3/r^5,...
                -3/5*x*beta^2/r^5 + 3*x*beta^2*y^2/r^7 + x/r^3 - 3*x*y^2/r^5,...
                -3/5*x*beta^2/r^5 + 3*x*beta^2*z^2/r^7 + x/r^3 - 3*x*z^2/r^5,...
                3*y*beta^2*x^2/r^7 - 3/5*y*beta^2/r^5 - 3*y*x^2/r^5,...
                3*z*beta^2*x^2/r^7 - 3/5*z*beta^2/r^5 - 3*z*x^2/r^5,...
                3*x*z*y*beta^2/r^7 - 3*x*z*y/r^5,...
                3*y*beta^2*x^2/r^7 - 3/5*y*beta^2/r^5 - 3*y*x^2/r^5,...
                3*z*beta^2*x^2/r^7 - 3/5*z*beta^2/r^5 - 3*z*x^2/r^5,...
                3*x*z*y*beta^2/r^7 - 3*x*z*y/r^5];...
                [3*y*beta^2*x^2/r^7 - 3/5*y*beta^2/r^5 + y/r^3 - 3*y*x^2/r^5,...
                -9/5*y*beta^2/r^5 + 3*y^3*beta^2/r^7 + y/r^3 - 3*y^3/r^5,...
                -3/5*y*beta^2/r^5 + 3*y*beta^2*z^2/r^7 + y/r^3 - 3*y*z^2/r^5,...
                -3/5*x*beta^2/r^5 + 3*x*beta^2*y^2/r^7 - 3*x*y^2/r^5,...
                3*x*z*y*beta^2/r^7 - 3*x*z*y/r^5,...
                3*z*beta^2*y^2/r^7 - 3/5*z*beta^2/r^5 - 3*z*y^2/r^5,...
                -3/5*x*beta^2/r^5 + 3*x*beta^2*y^2/r^7 - 3*x*y^2/r^5,...
                3*x*z*y*beta^2/r^7 - 3*x*z*y/r^5,...
                3*z*beta^2*y^2/r^7 - 3/5*z*beta^2/r^5 - 3*z*y^2/r^5];...
                [3*z*beta^2*x^2/r^7 - 3/5*z*beta^2/r^5 + z/r^3 - 3*z*x^2/r^5,...
                3*z*beta^2*y^2/r^7 - 3/5*z*beta^2/r^5 + z/r^3 - 3*z*y^2/r^5,...
                -9/5*z*beta^2/r^5 + 3*z^3*beta^2/r^7 + z/r^3 - 3*z^3/r^5,...
                3*x*z*y*beta^2/r^7 - 3*x*z*y/r^5,...
                -3/5*x*beta^2/r^5 + 3*x*beta^2*z^2/r^7 - 3*x*z^2/r^5,...
                -3/5*y*beta^2/r^5 + 3*y*beta^2*z^2/r^7 - 3*y*z^2/r^5,...
                3*x*z*y*beta^2/r^7 - 3*x*z*y/r^5,...
                -3/5*x*beta^2/r^5 + 3*x*beta^2*z^2/r^7 - 3*x*z^2/r^5,...
                -3/5*y*beta^2/r^5 + 3*y*beta^2*z^2/r^7 - 3*y*z^2/r^5]];
        end       
        function [res] = tens_K2_39(x1,alpha,beta)
            r=norm(x1);
            x=x1(1); y=x1(2); z=x1(3);
            r3=r^3; r5=r^5; r7=r^7;
            x2=x^2; y2=y^2; z2=z^2;
            x3=x^3; y3=y^3; z3=z^3;
            xzy=x*z*y;
            xz2=x*z2;
            xy2=x*y2;
            yx2=y*x2;
            yz2=y*z2;
            zy2=z*y2;
            zx2=z*x2;
            alpha2=alpha^2; beta2=beta^2;
            alpha2_r5=alpha2/r5; beta2_r5=beta2/r5;
            alpha2_r7=alpha2/r7; beta2_r7=beta2/r7;
            
            K2_112 = -3*yx2/r5 - 3/5*y*beta2_r5 - y*alpha2_r5 + 3*yx2*beta2_r7 + 5*yx2*alpha2_r7;
            K2_113 = -3*zx2/r5 - 3/5*z*beta2_r5 - z*alpha2_r5 + 3*zx2*beta2_r7 + 5*zx2*alpha2_r7;
            K2_123 = 5*xzy*alpha2_r7 + 3*xzy*beta2_r7 - 3*xzy/r5;

            K2_212 = -3*xy2/r5 - 3/5*x*beta2_r5 - x*alpha2_r5 + 3*xy2*beta2_r7 + 5*xy2*alpha2_r7;
            K2_223 = -3*zy2/r5 - 3/5*z*beta2_r5 - z*alpha2_r5 + 3*zy2*beta2_r7 + 5*zy2*alpha2_r7;

            K2_313 = -3*xz2/r5 - 3/5*x*beta2_r5 - x*alpha2_r5 + 3*xz2*beta2_r7 + 5*xz2*alpha2_r7;
            K2_323 = -3*yz2/r5 - 3/5*y*beta2_r5 - y*alpha2_r5 + 3*yz2*beta2_r7 + 5*yz2*alpha2_r7;
            
            res=[[x/r3 - 3*x3/r5 - 9/5*x*beta2_r5 - 3*x*alpha2_r5 + 5*x3*alpha2_r7 + 3*x3*beta2_r7,...
                x/r3 - 3/5*x*beta2_r5 - x*alpha2_r5 - 3*xy2/r5 + 5*xy2*alpha2_r7 + 3*xy2*beta2_r7,...
                x/r3 - 3/5*x*beta2_r5 - x*alpha2_r5 - 3*xz2/r5 + 5*xz2*alpha2_r7 + 3*xz2*beta2_r7,...
                K2_112,...
                K2_113,...
                K2_123,...
                K2_112,...
                K2_113,...
                K2_123];...
                [y/r3 - 3/5*y*beta2_r5 - y*alpha2_r5 - 3*yx2/r5 + 5*yx2*alpha2_r7 + 3*yx2*beta2_r7,...
                y/r3 - 9/5*y*beta2_r5 - 3*y*alpha2_r5 - 3*y3/r5 + 5*y3*alpha2_r7 + 3*y3*beta2_r7,...
                y/r3 - 3/5*y*beta2_r5 - y*alpha2_r5 - 3*yz2/r5 + 5*yz2*alpha2_r7 + 3*yz2*beta2_r7,...
                K2_212,...
                K2_123,...
                K2_223,...
                K2_212,...
                K2_123,...
                K2_223];...
                [z/r3 - 3/5*z*beta2_r5 - z*alpha2_r5 - 3*zx2/r5 + 5*zx2*alpha2_r7 + 3*zx2*beta2_r7,...
                z/r3 - 3/5*z*beta2_r5 - z*alpha2_r5 - 3*zy2/r5 + 5*zy2*alpha2_r7 + 3*zy2*beta2_r7,...
                z/r3 - 9/5*z*beta2_r5 - 3*z*alpha2_r5 - 3*z3/r5 + 5*z3*alpha2_r7 + 3*z3*beta2_r7,...
                K2_123,...
                K2_313,...
                K2_323,...
                K2_123,...
                K2_313,...
                K2_323]];
        end        
        function [res] = tens_K3_39(x1)
            r=norm(x1);
            x=x1(1); y=x1(2); z=x1(3);

            res=[[0, -3*y*z/r^5, 3*y*z/r^5, -3/2*x*z/r^5, 3/2*x*y/r^5, 3/2*y^2/r^5 - 3/2*z^2/r^5, -3/2*x*z/r^5, 3/2*x*y/r^5, 3/2*y^2/r^5 - 3/2*z^2/r^5];...
                [3*x*z/r^5, 0, -3*x*z/r^5, 3/2*y*z/r^5, -3/2*x^2/r^5 + 3/2*z^2/r^5, -3/2*x*y/r^5, 3/2*y*z/r^5, -3/2*x^2/r^5 + 3/2*z^2/r^5, -3/2*x*y/r^5];...
                [-3*x*y/r^5, 3*x*y/r^5, 0, 3/2*x^2/r^5 - 3/2*y^2/r^5, -3/2*y*z/r^5, 3/2*x*z/r^5, 3/2*x^2/r^5 - 3/2*y^2/r^5, -3/2*y*z/r^5, 3/2*x*z/r^5]];
        end        
        function [res] = tens_K4_99(x1,alpha,beta)
            r=norm(x1);
            x=x1(1); y=x1(2); z=x1(3);
            
            a2 = alpha^2; b2 = beta^2;
            a2b2 = a2*b2;
            
            r3 = r^3; r5 = r^5; r7 = r^7; r9 = r^9; r11 = r^11;
            
            x2 = x^2; y2 = y^2; z2 = z^2;
            x3 = x^3; y3 = y^3; z3 = z^3;
            x4 = x2^2; y4 = y2^2; z4 = z2^2;
            
            xy = x*y; xz = x*z; yz = y*z;
            x2y2 = x2*y2; y2z2 = y2*z2; x2z2 = x2*z2;
            
            xy3 = x*y^3;
            x3y = x3*y;
            yz3 = y*z3;
            xz3 = x*z3;
            x3z = x3*z;
            
            a2b2x2 = a2b2*x2; a2b2y2 = a2b2*y2; a2b2z2 = a2b2*z2;
            a2_r9 = a2/r9; b2_r9 = b2/r9;
            
            x2_r7 = x2/r7; y2_r7 = y2/r7; z2_r7 = z2/r7;
            x2_r9 = x2/r9; y2_r9 = y2/r9; z2_r9 = z2/r9;
            xy_r5 = xy/r5;
            
            r11v20 = r11*20;
            r11v10 = r11*10;
            r5v5 = r5*5;
            a2b2v189 = a2b2*189;
            a2_r7 = a2/r7;
            
            a2x2_r7 = a2*x2_r7; b2x2_r7 = b2*x2_r7; 
            a2y2_r7 = a2*y2_r7; b2y2_r7 = b2*y2_r7;
            a2z2_r7 = a2*z2_r7; b2z2_r7 = b2*z2_r7;
            a2x2_r9 = a2*x2_r9; b2x2_r9 = b2*x2_r9;
            
            xyb2 = xy*b2;
            xzb2 = xz*b2;
            xya2 = xy*a2;
            xza2 = xz*a2; 
            
            xyb2_r7 = xyb2/r7;
            xya2_r7 = xya2/r7;
            x2y2_r9 = x2y2/r9;
            
            K4_1122 = -(a2b2v189)/(20*r7) + (189*a2b2x2)/(10*r9) + (189*a2b2y2)/(10*r9) - (a2b2v189*x4)/(r11v20) - (189*a2b2x2*y2)/(r11v10) - (a2b2v189*y4)/(r11v20) + (a2b2v189*z4)/(r11v20) - (3*a2)/(r5v5) + 3*a2x2_r7 + 3*a2y2_r7 - 21*a2*x2y2_r9 - (3*b2)/(r5v5) + 3*b2x2_r7 + 3*b2y2_r7 - 21*b2*x2y2_r9 + 1/r3 - 3*x2/r5 - 3*y2/r5 + 15*x2y2/r7;
            K4_1133 = -(a2b2v189)/(20*r7) + (189*a2b2x2)/(10*r9) + (189*a2b2z2)/(10*r9) - (a2b2v189*x4)/(r11v20) - (189*a2b2x2*z2)/(r11v10) + (a2b2v189*y4)/(r11v20) - (a2b2v189*z4)/(r11v20) - (3*a2)/(r5v5) + 3*a2x2_r7 + 3*a2z2_r7 - 21*a2*x2z2/r9 - (3*b2)/(r5v5) + 3*b2x2_r7 + 3*b2z2_r7 - 21*b2*x2z2/r9 + 1/r3 - 3*x2/r5 - 3*z2/r5 + 15*x2z2/r7;
            K4_1112 = -6*xy_r5 + 15*x3y/r7 + 9*xyb2_r7 + 9*xya2_r7 - 21*x3y*a2_r9 - 21*x3y*b2_r9;
            K4_1113 = -6*xz/r5 + 15*x3z/r7 + 9*xzb2/r7 + 9*xza2/r7 - 21*x3z*a2_r9 - 21*x3z*b2_r9;
            K4_1123 = 3*yz*a2_r7 - 21*yz*a2x2_r9 + 3*yz*b2/r7 - 21*yz*b2x2_r9 - 3*yz/r5 + 15*yz*x2_r7;
            %
            K4_2233 = -(a2b2v189)/(20*r7) + (189*a2b2y2)/(10*r9) + (189*a2b2z2)/(10*r9) + (a2b2v189*x4)/(r11v20) - (a2b2v189*y4)/(r11v20) - (189*y2z2*b2*a2)/(r11v10) - (a2b2v189*z4)/(r11v20) - (3*a2)/(r5v5) + 3*a2y2_r7 + 3*a2z2_r7 - 21*a2*y2z2/r9 - (3*b2)/(r5v5) + 3*b2y2_r7 + 3*b2z2_r7 - 21*b2*y2z2/r9 + 1/r3 - 3*y2/r5 - 3*z2/r5 + 15*y2z2/r7;
            K4_2212 = -6*xy_r5 + 15*xy3/r7 + 9*xyb2_r7 + 9*xya2_r7 - 21*xy3*a2_r9 - 21*xy3*b2_r9;
            K4_2213 = 3*xza2/r7 - 21*xza2*y2_r9 + 3*xzb2/r7 - 21*xzb2*y2_r9 - 3*xz/r5 + 15*xz*y2_r7;
            K4_2223 = -6*yz/r5 + 15*z*y3/r7 + 9*yz*b2/r7 + 9*yz*a2_r7 - 21*z*y3*a2_r9 - 21*z*y3*b2_r9;
            %
            K4_3312 = 3*xya2_r7 - 21*xya2*z2_r9 + 3*xyb2_r7 - 21*xyb2*z2_r9 - 3*xy_r5 + 15*xy*z2_r7;
            K4_3313 = -6*xz/r5 + 15*xz3/r7 + 9*xzb2/r7 + 9*xza2/r7 - 21*xz3*a2_r9 - 21*xz3*b2_r9;
            K4_3323 = -6*yz/r5 + 15*yz3/r7 + 9*yz*b2/r7 + 9*yz*a2_r7 - 21*yz3*a2_r9 - 21*yz3*b2_r9;
            %
            K4_1212 = -(3*a2)/(r5v5) + 3*a2x2_r7 + 3*a2y2_r7 - 21*a2*x2y2_r9 - (3*b2)/(r5v5) + 3*b2x2_r7 + 3*b2y2_r7 - 21*b2*x2y2_r9 - (3*x2)/(2*r5) - (3*y2)/(2*r5) + 15*x2y2/r7;
            K4_1213 = 3*yz*a2_r7 - 21*yz*a2x2_r9 + 3*yz*b2/r7 - 21*yz*b2x2_r9 - (3*yz)/(2*r5) + 15*yz*x2_r7;
            K4_1223 = 3*xza2/r7 - 21*xza2*y2_r9 + 3*xzb2/r7 - 21*xzb2*y2_r9 - (3*xz)/(2*r5) + 15*xz*y2_r7;
            %
            K4_1313 = -(3*a2)/(r5v5) + 3*a2x2_r7 + 3*a2z2_r7 - 21*a2*x2z2/r9 - (3*b2)/(r5v5) + 3*b2x2_r7 + 3*b2z2_r7 - 21*b2*x2z2/r9 - (3*x2)/(2*r5) - (3*z2)/(2*r5) + 15*x2z2/r7;
            K4_1323 = 3*xya2_r7 - 21*xya2*z2_r9 + 3*xyb2_r7 - 21*xyb2*z2_r9 - (3*xy)/(2*r5) + 15*xy*z2_r7;
            %
            K4_2323 = -(3*a2)/(r5v5) + 3*a2y2_r7 + 3*a2z2_r7 - 21*a2*y2z2/r9 - (3*b2)/(r5v5) + 3*b2y2_r7 + 3*b2z2_r7 - 21*b2*y2z2/r9 - (3*y2)/(2*r5) - (3*z2)/(2*r5) + 15*y2z2/r7;
            
            res=[[-(9*a2)/(r5v5) + 18*a2x2_r7 - 21*a2*x4/r9 - (9*b2)/(r5v5) + 18*b2x2_r7 - 21*b2*x4/r9 + 1/r3 - 12*x2/r5 + 15*x4/r7,...
                K4_1122,...
                K4_1133,...
                K4_1112,...
                K4_1113,...
                K4_1123,...
                K4_1112,...
                K4_1113,...
                K4_1123];...
                % 2nd line
                [K4_1122,...
                -(9*a2)/(r5v5) + 18*a2y2_r7 - 21*a2*y4/r9 - (9*b2)/(r5v5) + 18*b2y2_r7 - 21*b2*y4/r9 + 1/r3 - 12*y2/r5 + 15*y4/r7,...
                K4_2233,...
                K4_2212,...
                K4_2213,...
                K4_2223,...
                K4_2212,...
                K4_2213,...
                K4_2223];...
                % 3rd line
                [K4_1133,...
                K4_2233,...
                -(9*a2)/(r5v5) + 18*a2z2_r7 - 21*a2*z4/r9 - (9*b2)/(r5v5) + 18*b2z2_r7 - 21*b2*z4/r9 + 1/r3 - 12*z2/r5 + 15*z4/r7,...
                K4_3312,...
                K4_3313,...
                K4_3323,...
                K4_3312,...
                K4_3313,...
                K4_3323];...
                % 4th line
                [K4_1112,...
                K4_2212,...
                K4_3312,...
                K4_1212,...
                K4_1213,...
                K4_1223,...
                K4_1212,...
                K4_1213,...
                K4_1223];...
                % 5th line
                [K4_1113,...
                K4_2213,...
                K4_3313,...
                K4_1213,...
                K4_1313,...
                K4_1323,...
                K4_1213,...
                K4_1313,...
                K4_1323];...
                % 6th line
                [K4_1123,...
                K4_2223,...
                K4_3323,...
                K4_1223,...
                K4_1323,...
                K4_2323,...
                K4_1223,...
                K4_1323,...
                K4_2323];...
                % 7th line
                [K4_1112,...
                K4_2212,...
                K4_3312,...
                K4_1212,...
                K4_1213,...
                K4_1223,...
                K4_1212,...
                K4_1213,...
                K4_1223];...
                % 8th line
                [K4_1113,...
                K4_2213,...
                K4_3313,...
                K4_1213,...
                K4_1313,...
                K4_1323,...
                K4_1213,...
                K4_1313,...
                K4_1323];...
                % 9th line
                [K4_1123,...
                K4_2223,...
                K4_3323,...
                K4_1223,...
                K4_1323,...
                K4_2323,...
                K4_1223,...
                K4_1323,...
                K4_2323]];
        end 
        
        %%% wall tensors %%%
        function [res] = tens_J1star(x1,x_alpha, beta) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
            
            rs3 = rs^3; rs5 = rs^5; rs7 = rs^7;
            b2 = beta^2;

            J1star23 = -ys*zs/rs3 + b2*ys*zs/rs5 - 6*x^2*ys*zs/rs5 + 10*b2*x*xs*ys*zs/rs7 - 6*x*xs*ys*zs/rs5;
            
            res = [-1/rs - b2/(3*rs3) - xs^2/rs3 + b2*xs^2/rs5 - 2*x^2/rs3 + 6*x^2*xs^2/rs5 - 2*x*xs/rs3 + 6*b2*x*xs/rs5 - 10*b2*x*xs^3/rs7 + 6*x*xs^3/rs5, ... 
                -xs*ys/rs3 + b2*xs*ys/rs5 + 6*x^2*xs*ys/rs5 - 2*x*ys/rs3 + 2*b2*x*ys/rs5 - 10*b2*x*xs^2*ys/rs7 + 6*x*xs^2*ys/rs5,... 
                -xs*zs/rs3 + b2*xs*zs/rs5 + 6*x^2*xs*zs/rs5 - 2*x*zs/rs3 + 2*b2*x*zs/rs5 - 10*b2*x*xs^2*zs/rs7 + 6*x*xs^2*zs/rs5;... 
                % 2nd line
                -xs*ys/rs3 + b2*xs*ys/rs5 - 6*x^2*xs*ys/rs5 - 2*x*ys/rs3 - 2*b2*x*ys/rs5 + 10*b2*x*xs^2*ys/rs7 - 6*x*xs^2*ys/rs5,... 
                -1/rs - b2/(3*rs3) - ys^2/rs3 + b2*ys^2/rs5 + 2*x^2/rs3 - 6*x^2*ys^2/rs5 + 2*x*xs/rs3 - 2*b2*x*xs/rs5 + 10*b2*x*xs*ys^2/rs7 - 6*x*xs*ys^2/rs5,...
                J1star23;...
                % 3rd line
                -xs*zs/rs3 + b2*xs*zs/rs5 - 6*x^2*xs*zs/rs5 - 2*x*zs/rs3 - 2*b2*x*zs/rs5 + 10*b2*x*xs^2*zs/rs7 - 6*x*xs^2*zs/rs5,... 
                J1star23,...
                -1/rs - b2/(3*rs3) - zs^2/rs3 + b2*zs^2/rs5 + 2*x^2/rs3 - 6*x^2*zs^2/rs5 + 2*x*xs/rs3 - 2*b2*x*xs/rs5 + 10*b2*x*xs*zs^2/rs7 - 6*x*xs*zs^2/rs5];
        end
        function [res] = tens_J2star(x1,x_alpha, alpha, beta) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a; 
            
            a2 = alpha^2; b2 = beta^2; a2b2 = a2*b2;
            rs3 = rs^3; rs5 = rs^5; rs7 = rs^7; rs9 = rs^9;
            xs2 = xs^2; ys2 = ys^2; zs2 = zs^2; x2 = x^2;
            xs3 = xs^3;
            xs4 = xs^4;
            xys = xs*ys; xzs = xs*zs; yzs = ys*zs;
        
        
            J2star23 = -6*xs*x*ys*zs/rs5 - yzs/rs3 + (a2 + b2)*yzs/rs5 - 6*x2*yzs/rs5 - 10*a2*yzs*xs2/rs7 + (70*a2b2*yzs*xs2)/(3*rs9) - (10*a2b2*yzs)/(3*rs7) - 10*(a2 - b2)*x*xs*ys*zs/rs7;
            
            res = [10*a2*xs4/rs7 + (b2 - 5*a2)*xs2/rs5 + 20*a2b2*xs2/rs7 - 2*a2b2/rs5 - (70*a2b2*xs4)/(3*rs9) - (a2 + b2)/(3*rs3) + 6*x2*xs2/rs5 + 6*x*xs3/rs5 - 2*xs*x/rs3 - 2*x2/rs3 - xs2/rs3 - 1/rs + 10*(a2 - b2)*x*xs3/rs7 - 6*(a2 - b2)*x*xs/rs5,...
                -(70*a2b2*xs3*ys)/(3*rs9) + 10*a2b2*xys/rs7 + 10*a2*xs3*ys/rs7 + (b2 - 3*a2)*xys/rs5 + 6*x2*xys/rs5 + 6*x*xs2*ys/rs5 - xys/rs3 - 2*x*ys/rs3 - 2*(a2 - b2)*x*ys/rs5 + 10*(a2 - b2)*x*xs2*ys/rs7,...
                -(70*a2b2*xs3*zs)/(3*rs9) + 10*a2b2*xzs/rs7 + 10*a2*xs3*zs/rs7 + (b2 - 3*a2)*xzs/rs5 + 6*x2*xzs/rs5 + 6*x*xs2*zs/rs5 - xzs/rs3 - 2*x*zs/rs3 - 2*(a2 - b2)*x*zs/rs5 + 10*(a2 - b2)*x*xs2*zs/rs7;...
                %
                (70*a2b2*xs3*ys)/(3*rs9) - 10*a2b2*xys/rs7 - 10*a2*xs3*ys/rs7 + (a2 + b2)*xys/rs5 - 6*x2*xys/rs5 - 6*x*xs2*ys/rs5 - xys/rs3 - 2*x*ys/rs3 + 2*(a2 - b2)*x*ys/rs5 - 10*(a2 - b2)*x*xs2*ys/rs7,...
                -(a2 + b2)/(3*rs3) - 8*a2b2/(3*rs5) + (70*a2b2*xs2*ys2)/(3*rs9) - 1/rs - 6*x2*ys2/rs5 - 6*xs*x*ys2/rs5 + 2*a2*xs2/rs5 + (a2 + b2)*ys2/rs5 + (10*a2b2*zs2)/(3*rs7) + 2*x2/rs3 + 2*xs*x/rs3 - 10*a2*xs2*ys2/rs7 - ys2/rs3 + 2*(a2 - b2)*x*xs/rs5 - 10*(a2 - b2)*x*xs*ys2/rs7,...
                J2star23;...
                %
                (70*a2b2*xs3*zs)/(3*rs9) - 10*a2b2*xzs/rs7 - 10*a2*xs3*zs/rs7 + (a2 + b2)*xzs/rs5 - 6*x2*xzs/rs5 - 6*x*xs2*zs/rs5 - xzs/rs3 - 2*x*zs/rs3 + 2*(a2 - b2)*x*zs/rs5 - 10*(a2 - b2)*x*xs2*zs/rs7,...
                J2star23,...
                -(a2 + b2)/(3*rs3) - 8*a2b2/(3*rs5) + (70*a2b2*xs2*zs2)/(3*rs9) - 1/rs - 6*x2*zs2/rs5 - 6*xs*x*zs2/rs5 + 2*a2*xs2/rs5 + (a2 + b2)*zs2/rs5 + (10*a2b2*ys2)/(3*rs7) + 2*x2/rs3 + 2*xs*x/rs3 - 10*a2*xs2*zs2/rs7 - zs2/rs3 + 2*(a2 - b2)*x*xs/rs5 - 10*(a2 - b2)*x*xs*zs2/rs7];
        
        end
        function [res] = tens_J3star(x1,x_alpha, beta) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
            
            b2 = beta^2;
            rs3 = rs^3; rs5 = rs^5; rs7 = rs^7;
            xs2 = xs^2; ys2 = ys^2; zs2 = zs^2;
            xyzs = xs*ys*zs;
            
            J3star22 = 6*xyzs/rs5 - 10*b2*xyzs/rs7 + 6*x*ys*zs/rs5;

            res =[0,...
                -zs/rs3,...
                ys/rs3;...
                %
                zs/rs3 + 2*b2*zs/rs5 + 6*xs2*zs/rs5 - 10*b2*xs2*zs/rs7 + 6*x*xs*zs/rs5,...
                J3star22,...
                6*xs*zs2/rs5 - 10*b2*xs*zs2/rs7 + 6*x*zs2/rs5 - 2*x/rs3 + 2*b2*xs/rs5 - xs/rs3;... 
                %
                -6*xs2*ys/rs5 + 10*b2*xs2*ys/rs7 - ys/rs3 - 2*b2*ys/rs5 - 6*x*xs*ys/rs5,...
                -6*x*ys2/rs5 + 10*b2*ys2*xs/rs7 - 6*ys2*xs/rs5 + 2*x/rs3 - 2*b2*xs/rs5 + xs/rs3,...
                -J3star22];
        
        end
        function [res] = tens_J4star93(x1,x_alpha, alpha, beta) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
            
            a2 = alpha^2; b2 = beta^2; a2b2 = a2*b2;
            rs3 = rs^3; rs5 = rs^5; rs7 = rs^7; rs9 = rs^9; rs11 = rs^11;
        
            x2 = x^2; 
            xs2 = xs^2; ys2 = ys^2; zs2 = zs^2;
            xs3 = xs^3; ys3 = ys^3; zs3 = zs^3;
            xs4 = xs^4; 
            xs5 = xs^5;
        
            xys = xs*ys; xzs = xs*zs; yzs = ys*zs;
            xyzs = xs*ys*zs;
            
            J4star121 = 6*x2*ys/rs5 - 42*a2*xs4*ys/rs9 - 84*a2b2*xs2*ys/rs9 + (15*a2 + 5*b2)*xs2*ys/rs7 - (42*a2 - 70*b2)*x*xs3*ys/rs9 + (18*a2 - 30*b2)*x*xys/rs7 + 126*a2b2*xs4*ys/rs11 + (3*a2 - 5*b2)*ys/(5*rs5) + 6*a2b2*ys/rs7 - 3*xs2*ys/rs5 - 30*x2*xs2*ys/rs7 - 30*x*xs3*ys/rs7; %%%% checked
            J4star122 = -(9*a2 + 5*b2)*xs/(5*rs5) + 6*a2b2*xs/rs7 + 6*a2*xs3/rs7 - 14*a2b2*xs3/rs9 - 3*ys2*xs/rs5 + 6*x2*xs/rs5 - 30*x2*ys2*xs/rs7 - 30*x*xs2*ys2/rs7 - (6*a2 - 10*b2)*x*zs2/rs7 + (9*a2 + 5*b2)*ys2*xs/rs7 - (42*a2 - 70*b2)*x*xs2*ys2/rs9 + (24*a2 - 40*b2)*x/(5*rs5) + 126*a2b2*ys2*xs3/rs11 + 6*x*xs2/rs5 - 42*a2*ys2*xs3/rs9 - 42*a2b2*xs*ys2/rs9; %%%% checked
            J4star123 = -(42*a2 - 70*b2)*x*yzs*xs2/rs9 - 30*x2*xyzs/rs7 - 30*x*yzs*xs2/rs7 - 42*a2b2*xyzs/rs9 - 42*a2*yzs*xs3/rs9 + 126*a2b2*yzs*xs3/rs11 + (6*a2 - 10*b2)*x*yzs/rs7 - 3*xyzs/rs5 + (9*a2 + 5*b2)*xyzs/rs7; %%%% checked
            
            J4star131 = 6*x2*zs/rs5 - 42*a2*xs4*zs/rs9 - 84*a2b2*xs2*zs/rs9 + (15*a2 + 5*b2)*xs2*zs/rs7 - (42*a2 - 70*b2)*x*xs3*zs/rs9 + (18*a2 - 30*b2)*x*xzs/rs7 + 126*a2b2*xs4*zs/rs11 + (3*a2 - 5*b2)*zs/(5*rs5) + 6*a2b2*zs/rs7 - 3*xs2*zs/rs5 - 30*x2*xs2*zs/rs7 - 30*x*xs3*zs/rs7; %%%% checked
            J4star133 = -(9*a2 + 5*b2)*xs/(5*rs5) + 6*a2b2*xs/rs7 + 6*a2*xs3/rs7 - 14*a2b2*xs3/rs9 - 3*xs*zs2/rs5 + 6*x2*xs/rs5 - 30*x2*zs2*xs/rs7 - 30*x*xs2*zs2/rs7 - (6*a2 - 10*b2)*x*ys2/rs7 + (9*a2 + 5*b2)*xs*zs2/rs7 - (42*a2 - 70*b2)*x*xs2*zs2/rs9 + (24*a2 - 40*b2)*x/(5*rs5) + 126*a2b2*zs2*xs3/rs11 + 6*x*xs2/rs5 - 42*a2*zs2*xs3/rs9 - 42*a2b2*xs*zs2/rs9; %%%% checked
        
            J4star231 = (42*a2 - 70*b2)*x*yzs*xs2/rs9 + 30*x2*xyzs/rs7 + 30*x*yzs*xs2/rs7 + 42*a2b2*xyzs/rs9 + 42*a2*yzs*xs3/rs9 - 126*a2b2*yzs*xs3/rs11 + 6*x*yzs/rs5 - (6*a2 - 10*b2)*x*yzs/rs7 + 3*xyzs/rs5 - (3*a2 + 5*b2)*xyzs/rs7; %%%% checked
            J4star232 = -(3*a2 + 5*b2)*ys2*zs/rs7 + (3*a2 + 5*b2)*zs/(5*rs5) + 42*a2*xs2*zs*ys2/rs9 + (42*a2 - 70*b2)*x*xzs*ys2/rs9 + 3*ys2*zs/rs5 - 126*a2b2*xs2*zs*ys2/rs11 + 12*a2b2*zs/rs7 - 14*a2b2*zs3/rs9 - 6*a2*xs2*zs/rs7 - 6*x2*zs/rs5 + 30*x*xzs*ys2/rs7 - 6*x*xzs/rs5 - (6*a2 - 10*b2)*x*xzs/rs7 + 30*x2*zs*ys2/rs7; %%%% checked
            J4star233 = -(3*a2 + 5*b2)*zs2*ys/rs7 + (3*a2 + 5*b2)*ys/(5*rs5) + 42*a2*xs2*ys*zs2/rs9 + (42*a2 - 70*b2)*x*xys*zs2/rs9 + 3*zs2*ys/rs5 - 126*a2b2*xs2*ys*zs2/rs11 + 12*a2b2*ys/rs7 - 14*a2b2*ys3/rs9 - 6*a2*xs2*ys/rs7 - 6*x2*ys/rs5 + 30*x*xys*zs2/rs7 - 6*x*xys/rs5 - (6*a2 - 10*b2)*x*xys/rs7 + 30*x2*ys*zs2/rs7; %%%% checked
            
            res = [(18*a2 - 30*b2)*x/(5*rs5) - 126*a2b2*xs5/rs11 + 3*xs3/rs5 + 42*a2*xs5/rs9 - 2*x/rs3 + 140*a2b2*xs3/rs9 - (33*a2 + 5*b2)*xs3/rs7 + (42*a2 - 70*b2)*x*xs4/rs9 - 12*x*xs2/rs5 - 18*x2*xs/rs5 + 30*x2*xs3/rs7 + 30*x*xs4/rs7 + (-36*a2 + 60*b2)*x*xs2/rs7 - xs/rs3 - 30*a2b2*xs/rs7 + (9*a2 + 15*b2)*xs/(5*rs5),... %%%% checked
                -(21*a2 + 5*b2)*xs2*ys/rs7 - ys/rs3 + (3*a2 + 5*b2)*ys/(5*rs5) - 6*a2b2*ys/rs7 - 126*a2b2*xs4*ys/rs11 + (42*a2 - 70*b2)*x*xs3*ys/rs9 - 6*x*xys/rs5 + (-18*a2 + 30*b2)*x*xys/rs7 + 84*a2b2*xs2*ys/rs9 + 3*xs2*ys/rs5 - 6*x2*ys/rs5 + 42*a2*xs4*ys/rs9 + 30*x2*xs2*ys/rs7 + 30*x*xs3*ys/rs7,... %%%% checked
                -(21*a2 + 5*b2)*xs2*zs/rs7 - zs/rs3 + (3*a2 + 5*b2)*zs/(5*rs5) - 6*a2b2*zs/rs7 - 126*a2b2*xs4*zs/rs11 + (42*a2 - 70*b2)*x*xs3*zs/rs9 - 6*x*xzs/rs5 + (-18*a2 + 30*b2)*x*xzs/rs7 + 84*a2b2*xs2*zs/rs9 + 3*xs2*zs/rs5 - 6*x2*zs/rs5 + 42*a2*xs4*zs/rs9 + 30*x2*xs2*zs/rs7 + 30*x*xs3*zs/rs7;... %%%% checked
                % 2nd line
                (6*a2 - 10*b2)*x*zs2/rs7 - (3*a2 + 5*b2)*ys2*xs/rs7 + 30*x2*ys2*xs/rs7 + 30*x*xs2*ys2/rs7 - (24*a2 - 40*b2)*x/(5*rs5) + 6*x*ys2/rs5 - 2*x/rs3 - 126*a2b2*ys2*xs3/rs11 + 3*ys2*xs/rs5 - 6*x*xs2/rs5 - 6*x2*xs/rs5 + (42*a2 - 70*b2)*x*xs2*ys2/rs9 + 42*a2b2*xs*ys2/rs9 - xs/rs3 + 42*a2*ys2*xs3/rs9 + 14*a2b2*xs3/rs9 - 6*a2*xs3/rs7 - 6*a2b2*xs/rs7 + (3*a2 + 5*b2)*xs/(5*rs5),... %%%% checked
                -18*x2*ys/rs5 + 42*a2*xs2*ys3/rs9 + 42*a2b2*xs2*ys/rs9 - 126*a2b2*xs2*ys3/rs11 - 18*a2*xs2*ys/rs7 + 30*x2*ys3/rs7 - 18*x*xys/rs5 - (18*a2 - 30*b2)*x*xys/rs7 + (42*a2 - 70*b2)*xs*x*ys3/rs9 + 30*xs*x*ys3/rs7 + 3*ys3/rs5 - (3*a2 + 5*b2)*ys3/rs7 + 14*a2b2*ys3/rs9 + (9*a2 + 15*b2)*ys/(5*rs5) - 6*a2b2*ys/rs7 - ys/rs3,... %%%% checked
                -6*a2*xs2*zs/rs7 + 30*x*xzs*ys2/rs7 - 6*x*xzs/rs5 - (6*a2 - 10*b2)*x*xzs/rs7 + 3*ys2*zs/rs5 - (3*a2 + 5*b2)*ys2*zs/rs7 + 30*x2*zs*ys2/rs7 + (42*a2 - 70*b2)*x*xzs*ys2/rs9 + 42*a2*xs2*zs*ys2/rs9 - 126*a2b2*xs2*zs*ys2/rs11 + (3*a2 + 5*b2)*zs/(5*rs5) + 12*a2b2*zs/rs7 - 14*a2b2*zs3/rs9 - zs/rs3 - 6*x2*zs/rs5;... %%%% checked
                % 3rd line
                (6*a2 - 10*b2)*x*ys2/rs7 - (3*a2 + 5*b2)*xs*zs2/rs7 + 30*x2*zs2*xs/rs7 + 30*x*xs2*zs2/rs7 - (24*a2 - 40*b2)*x/(5*rs5) + 6*x*zs2/rs5 - 2*x/rs3 - 126*a2b2*zs2*xs3/rs11 + 3*xs*zs2/rs5 - 6*x*xs2/rs5 - 6*x2*xs/rs5 + (42*a2 - 70*b2)*x*xs2*zs2/rs9 + 42*a2b2*xs*zs2/rs9 - xs/rs3 + 42*a2*zs2*xs3/rs9 + 14*a2b2*xs3/rs9 - 6*a2*xs3/rs7 - 6*a2b2*xs/rs7 + (3*a2 + 5*b2)*xs/(5*rs5),... %%%% checked
                -6*a2*xs2*ys/rs7 + 30*x*xys*zs2/rs7 - 6*x*xys/rs5 - (6*a2 - 10*b2)*x*xys/rs7 + 3*zs2*ys/rs5 - (3*a2 + 5*b2)*zs2*ys/rs7 + 30*x2*ys*zs2/rs7 + (42*a2 - 70*b2)*x*xys*zs2/rs9 + 42*a2*xs2*ys*zs2/rs9 - 126*a2b2*xs2*ys*zs2/rs11 + (3*a2 + 5*b2)*ys/(5*rs5) + 12*a2b2*ys/rs7 - 14*a2b2*ys3/rs9 - ys/rs3 - 6*x2*ys/rs5,... %%%% checked
                -18*x2*zs/rs5 + 42*a2*xs2*zs3/rs9 + 42*a2b2*xs2*zs/rs9 - 126*a2b2*xs2*zs3/rs11 - 18*a2*xs2*zs/rs7 + 30*x2*zs3/rs7 - 18*x*xzs/rs5 - (18*a2 - 30*b2)*x*xzs/rs7 + (42*a2 - 70*b2)*xs*x*zs3/rs9 + 30*xs*x*zs3/rs7 + 3*zs3/rs5 - (3*a2 + 5*b2)*zs3/rs7 + 14*a2b2*zs3/rs9 + (9*a2 + 15*b2)*zs/(5*rs5) - 6*a2b2*zs/rs7 - zs/rs3;... %%%% checked
                % 4th line
                J4star121,... %%% see above
                J4star122,... %%% see above
                J4star123;... %%% see above
                % 5th line
                J4star131,... %%% see above
                J4star123,... %%% sym checked
                J4star133;... %%% see above
                % 6th line
                J4star231,... %%% see above
                J4star232,... %%% see above
                J4star233;... %%% see above
                % 7th line
                J4star121,... %%% sym checked
                J4star122,... %%% sym checked
                J4star123;... %%% sym checked
                % 8th line
                J4star131,... %%% sym checked
                J4star123,... %%% sym checked
                J4star133;... %%% sym checked
                % 9th line
                J4star231,... %%% sym checked
                J4star232,... %%% sym checked
                J4star233]; %%% sym checked
        
        end

        function [res] = tens_R1star(x1,x_alpha) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
            
            rs3 = rs^3; rs5 = rs^5;

            res = [0, -zs/rs3 + 6*x*xs*zs/rs5, ys/rs3 - 6*x*ys*xs/rs5;...
                zs/rs3, -6*x*ys*zs/rs5, -xs/rs3 - x*(2/rs3 - 6*ys^2/rs5);...
                -ys/rs3, xs/rs3 - x*(-2/rs3 + 6*zs^2/rs5), 6*x*ys*zs/rs5];
        end
        function [res] = tens_R2star(x1,x_alpha, alpha) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
            
            a2 = alpha^2; % alpha = radius of sphere alpha
            rs3 = rs^3; rs5 = rs^5; rs7 = rs^7;
            
            res = [0, -zs/rs3 + 6*x*xs*zs/rs5 + 10*a2*xs^2*zs/rs7 - 2*a2*zs/rs5, ys/rs3 - 6*x*ys*xs/rs5 - 10*a2*xs^2*ys/rs7 + 2*a2*ys/rs5;...
                zs/rs3, -6*x*ys*zs/rs5 - 10*a2*xs*ys*zs/rs7, -xs/rs3 - 2*x/rs3 + 6*x*ys^2/rs5 + 10*a2*xs*ys^2/rs7 - 2*a2*xs/rs5;...
                -ys/rs3, xs/rs3 + 2*x/rs3 - 6*x*zs^2/rs5 - 10*a2*xs*zs^2/rs7 + 2*a2*xs/rs5, 6*x*ys*zs/rs5 + 10*a2*xs*ys*zs/rs7];
        end
        function [res] = tens_R3star(x1,x_alpha) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
            
            rs3 = rs^3; rs5 = rs^5;
            
            res = [1/(2*rs3) - (3*xs^2)/(2*rs5), -(3*xs*ys)/(2*rs5), -(3*xs*zs)/(2*rs5);...
                (3*xs*ys)/(2*rs5), (9*zs^2)/(2*rs5) - 1/rs3 - (3*xs^2)/(2*rs5), -(9*ys*zs)/(2*rs5);...
                (3*xs*zs)/(2*rs5), -(9*ys*zs)/(2*rs5), (9*ys^2)/(2*rs5) - 1/rs3 - (3*xs^2)/(2*rs5)];
        end
        function [res] = tens_R4star93(x1,x_alpha, alpha) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
            
            a2 = alpha^2;
            rs5 = rs^5; rs7 = rs^7; rs9 = rs^9;
            xs2 = xs^2; ys2 = ys^2; zs2 = zs^2;
            xys = xs*ys; xzs = xs*zs; yzs = ys*zs;
            
            R4star121 = (3*xzs)/(2*rs5); %%% checked
            R4star122 = 6*a2*yzs/rs7 - (3*yzs)/(2*rs5) - 42*a2*yzs*xs2/rs9 - 30*x*ys*xs*zs/rs7; %%% checked
            R4star123 = (3*ys2)/(2*rs5) + 30*x*ys2*xs/rs7 - 6*xs*x/rs5 + 42*a2*ys2*xs2/rs9 - 24*a2/(5*rs5) - (3*xs2)/(2*rs5) + 6*a2*zs2/rs7; %%% checked
            
            R4star131 = -(3*xys)/(2*rs5); %%% checked
            R4star132 = -(3*zs2)/(2*rs5) - 30*x*zs2*xs/rs7 + 6*xs*x/rs5 - 42*a2*zs2*xs2/rs9 + 24*a2/(5*rs5) + (3*xs2)/(2*rs5) - 6*a2*ys2/rs7; %%% checked
            
            R4star231 = -(3*zs2)/(2*rs5) + (3*ys2)/(2*rs5); %%% checked
            R4star232 = -6*x*ys/rs5 - (3*xys)/(2*rs5) + 30*x*zs2*ys/rs7 - 6*a2*xys/rs7 + 42*a2*xys*zs2/rs9; %%% checked
            R4star233 = 6*x*zs/rs5 + (3*xzs)/(2*rs5) - 30*x*ys2*zs/rs7 + 6*a2*xzs/rs7 - 42*a2*xzs*ys2/rs9; %%% checked
        
            res = [0,... %
                3*xzs/rs5 - 6*x*zs/rs5 + 30*x*xs2*zs/rs7 - 18*a2*xzs/rs7 + 42*a2*xs^3*zs/rs9,... %
                -3*xys/rs5 + 6*x*ys/rs5 - 30*x*xs2*ys/rs7 + 18*a2*xys/rs7 - 42*a2*xs^3*ys/rs9;... %
                % 2nd line
                -3*yzs/rs5,... %
                -6*x*zs/rs5 + 30*x*ys2*zs/rs7 - 6*a2*xzs/rs7 + 42*a2*xzs*ys2/rs9,... %
                3*xys/rs5 + 18*x*ys/rs5 - 30*x*ys^3/rs7 + 18*a2*xys/rs7 - 42*a2*xs*ys^3/rs9;... %
                % 3rd line
                3*yzs/rs5,... %
                -3*xzs/rs5 - 18*x*zs/rs5 + 30*x*zs^3/rs7 - 18*a2*xzs/rs7 + 42*a2*xs*zs^3/rs9,... %
                6*x*ys/rs5 - 30*x*zs2*ys/rs7 + 6*a2*xys/rs7 - 42*a2*xys*zs2/rs9;... %
                % 4th line
                R4star121,... %
                R4star122,... %
                R4star123;... %
                % 5th line
                R4star131,... %
                R4star132,... %
                -R4star122;... %%% sym checked
                % 6th line
                R4star231,... %
                R4star232,... %
                R4star233;... %
                % 7th line
                R4star121,... %%% sym checked
                R4star122,... %%% sym checked
                R4star123;... %%% sym checked
                % 8th line
                R4star131,... %%% sym checked
                R4star132,... %%% sym checked
                -R4star122;... %%% sym checked
                % 9th line
                R4star231,... %%% sym checked
                R4star232,... %%% sym checked
                R4star233]; %%% sym checked
        
        end
        
        function [res] = tens_K1star39(x1,x_alpha, beta) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
        
            b2 = beta^2;

            rs3 = rs^3; rs5 = rs^5; rs7 = rs^7; rs9 = rs^9;
            
            x2 = x^2;
            xs2 = xs^2; ys2 = ys^2; zs2 = zs^2;
            xs3 = xs^3; ys3 = ys^3; zs3 = zs^3;
            
            xys = xs*ys; xy2s = xs*ys2; yx2s = ys*xs2;
            xzs = xs*zs; xz2s = xs*zs2; zx2s = zs*xs2;
            yzs = ys*zs; yz2s = ys*zs2; zy2s = zs*ys2;
        
            xyzs = xs*ys*zs; 
            x2yzs = xs2*ys*zs; xy2zs = xs*ys2*zs; xyz2s = xs*ys*zs2;
            
            
            K1star112 = 3*yx2s/rs5 - 3*b2*yx2s/rs7 + (3*b2*ys)/(5*rs5) + 6*x2*ys/rs5 - 30*x2*yx2s/rs7 + 12*x*xys/rs5 - 18*b2*x*xys/rs7 - 30*x*xs3*ys/rs7 + 42*b2*x*xs3*ys/rs9; %%%% checked
            K1star113 = 3*zx2s/rs5 - 3*b2*zx2s/rs7 + (3*b2*zs)/(5*rs5) + 6*x2*zs/rs5 - 30*x2*zx2s/rs7 + 12*x*xzs/rs5 - 18*b2*x*xzs/rs7 - 30*x*xs3*zs/rs7 + 42*b2*x*xs3*zs/rs9; %%%% checked
            K1star123 = 3*xyzs/rs5 - 3*b2*xyzs/rs7 - 30*x2*xyzs/rs7 + 6*x*yzs/rs5 - 6*b2*x*yzs/rs7 - 30*x*x2yzs/rs7 + 42*b2*x*x2yzs/rs9; %%%% checked
            
            K1star212 = 3*xy2s/rs5 - 3*b2*xy2s/rs7 + (3*b2*xs)/(5*rs5) - 6*x2*xs/rs5 + 30*x2*xy2s/rs7 - 6*x*xs2/rs5 + 6*b2*x*xs2/rs7 - (6*b2*x)/(5*rs5) + 30*x*xs2*ys2/rs7 - 42*b2*x*xs2*ys2/rs9 + 6*b2*x*ys2/rs7; %%%% checked
            K1star213 = 3*xyzs/rs5 - 3*b2*xyzs/rs7 + 30*x2*xyzs/rs7 + 30*x*x2yzs/rs7 - 42*b2*x*x2yzs/rs9 + 6*b2*x*yzs/rs7; %%%% checked
            K1star223 = 3*zy2s/rs5 - 3*b2*zy2s/rs7 + (3*b2*zs)/(5*rs5) + 30*x2*zy2s/rs7 - 6*x2*zs/rs5 - 6*x*xzs/rs5 + 6*b2*x*xzs/rs7 + 30*x*xy2zs/rs7 - 42*b2*x*xy2zs/rs9; %%%% checked
            
            K1star313 = 3*xz2s/rs5 - 3*b2*xz2s/rs7 + (3*b2*xs)/(5*rs5) - 6*x2*xs/rs5 + 30*x2*xz2s/rs7 - 6*x*xs2/rs5 + 6*b2*x*xs2/rs7 - (6*b2*x)/(5*rs5) + 30*x*xs2*zs2/rs7 - 42*b2*x*xs2*zs2/rs9 + 6*b2*x*zs2/rs7; %%%% checked
            K1star323 = 3*yz2s/rs5 - 3*b2*yz2s/rs7 + (3*b2*ys)/(5*rs5) + 30*x2*yz2s/rs7 - 6*x2*ys/rs5 - 6*x*xys/rs5 + 30*x*xyz2s/rs7 - 42*b2*x*xyz2s/rs9 + 6*b2*x*xys/rs7; %%%% checked
        
            res = [-xs/rs3 + 3*xs3/rs5 - 3*b2*xs3/rs7 + (9*b2*xs)/(5*rs5) - 30*x2*xs3/rs7 + 18*x2*xs/rs5 - 2*x/rs3 + 24*x*xs2/rs5 - 36*b2*x*xs2/rs7 + (18*b2*x)/(5*rs5) - 30*x*xs^4/rs7 + 42*b2*x*xs^4/rs9, ... %%%% checked
                -xs/rs3 + 3*xy2s/rs5 - 3*b2*xy2s/rs7 + (3*b2*xs)/(5*rs5) + 6*x2*xs/rs5 - 30*x2*xy2s/rs7 - 2*x/rs3 + 6*x*ys2/rs5 - 6*b2*x*ys2/rs7 + (6*b2*x)/(5*rs5) + 6*x*xs2/rs5 - 30*x*xs2*ys2/rs7 + 42*b2*x*xs2*ys2/rs9 - 6*b2*x*xs2/rs7,... %%%% checked
                -xs/rs3 + 3*xz2s/rs5 - 3*b2*xz2s/rs7 + (3*b2*xs)/(5*rs5) + 6*x2*xs/rs5 - 30*x2*xz2s/rs7 - 2*x/rs3 + 6*x*zs2/rs5 - 6*b2*x*zs2/rs7 + (6*b2*x)/(5*rs5) + 6*x*xs2/rs5 - 30*x*xs2*zs2/rs7 + 42*b2*x*xs2*zs2/rs9 - 6*b2*x*xs2/rs7,... %%%% checked
                K1star112,... %
                K1star113,... %
                K1star123,... %
                K1star112,... %%% sym checked
                K1star113,... %%% sym checked
                K1star123;... %%% sym checked
                % 2nd line
                -ys/rs3 + 3*yx2s/rs5 - 3*b2*yx2s/rs7 + (3*b2*ys)/(5*rs5) - 6*x2*ys/rs5 + 30*x2*yx2s/rs7 - 6*x*xys/rs5 + 30*x*xs3*ys/rs7 - 42*b2*x*xs3*ys/rs9 + 18*b2*x*xys/rs7,... %%%% checked
                -ys/rs3 + 3*ys3/rs5 - 3*b2*ys3/rs7 + (9*b2*ys)/(5*rs5) + 30*x2*ys3/rs7 - 18*x2*ys/rs5 - 18*x*xys/rs5 + 18*b2*x*xys/rs7 + 30*x*ys3*xs/rs7 - 42*b2*x*ys3*xs/rs9,... %%%% checked
                -ys/rs3 + 3*yz2s/rs5 - 3*b2*yz2s/rs7 + (3*b2*ys)/(5*rs5) - 6*x2*ys/rs5 + 30*x2*yz2s/rs7 - 6*x*xys/rs5 + 30*x*xyz2s/rs7 - 42*b2*x*xyz2s/rs9 + 6*b2*x*xys/rs7,... %%%% checked
                K1star212,... %
                K1star213,... %
                K1star223,... %
                K1star212,... %%% sym checked
                K1star213,... %%% sym checked
                K1star223;... %%% sym checked
                % 3rd line
                -zs/rs3 + 3*zx2s/rs5 - 3*b2*zx2s/rs7 + (3*b2*zs)/(5*rs5) + 30*x2*zx2s/rs7 - 6*x2*zs/rs5 - 6*x*xzs/rs5 + 30*x*xs3*zs/rs7 - 42*b2*x*xs3*zs/rs9 + 18*b2*x*xzs/rs7,... %%%% checked
                -zs/rs3 + 3*zy2s/rs5 - 3*b2*zy2s/rs7 + (3*b2*zs)/(5*rs5) + 30*x2*zy2s/rs7 - 6*x2*zs/rs5 - 6*x*xzs/rs5 + 6*b2*x*xzs/rs7 + 30*x*xy2zs/rs7 - 42*b2*x*xy2zs/rs9,... %%%% checked
                -zs/rs3 + 3*zs3/rs5 - 3*b2*zs3/rs7 + (9*b2*zs)/(5*rs5) + 30*x2*zs3/rs7 - 18*x2*zs/rs5 - 18*x*xzs/rs5 + 18*b2*x*xzs/rs7 + 30*x*zs3*xs/rs7 - 42*b2*x*zs3*xs/rs9,... %%%% checked
                K1star213,... %%% sym checked
                K1star313,... %
                K1star323,... %
                K1star213,... %%% sym checked
                K1star313,... %%% sym checked
                K1star323]; %%% sym checked
        end
        function [res] = tens_K2star39(x1,x_alpha, alpha, beta) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
            
            a2 = alpha^2; b2 = beta^2; a2b2 = a2*b2;
            rs3 = rs^3; rs5 = rs^5; rs7 = rs^7; rs9 = rs^9; rs11 = rs^11; 
            
            x2 = x^2;
            xs2 = xs^2; ys2 = ys^2; zs2 = zs^2;
            xs3 = xs^3; ys3 = ys^3; zs3 = zs^3;
            xs4 = xs^4; 
            xs5 = xs^5;
            
            xys = xs*ys; xy2s = xs*ys2; yx2s = ys*xs2; 
            xzs = xs*zs; xz2s = xs*zs2; zx2s = zs*xs2;
            yzs = ys*zs; yz2s = ys*zs2; zy2s = zs*ys2;
        
            xyzs = xs*ys*zs; 
            x2yzs = xs2*yzs; %xy2zs = xs*ys2*zs; xyz2s = xs*ys*zs2;
            x3yzs = xs3*yzs;
            x2y2s = xs2*ys2; x3y2s = xs3*ys2;
            x2z2s = xs2*zs2; x3z2s = xs3*zs2;
            
            
            K2star112 = -84*a2b2*yx2s/rs9 + (35*a2 - 3*b2)*yx2s/rs7 - 30*x2*yx2s/rs7 + 12*x*xys/rs5 - 30*x*xs3*ys/rs7 + 126*a2b2*ys*xs4/rs11 - 70*a2*ys*xs4/rs9 + 6*a2b2*ys/rs7 - (5*a2 - 3*b2)*ys/(5*rs5) + 3*yx2s/rs5 - (70*a2 - 42*b2)*x*xs3*ys/rs9 + (30*a2 - 18*b2)*x*xys/rs7 + 6*x2*ys/rs5; %%%% checked
            K2star113 = -84*a2b2*zx2s/rs9 + (35*a2 - 3*b2)*zx2s/rs7 - 30*x2*zx2s/rs7 + 12*x*xzs/rs5 - 30*x*xs3*zs/rs7 + 126*a2b2*zs*xs4/rs11 - 70*a2*zs*xs4/rs9 + 6*a2b2*zs/rs7 - (5*a2 - 3*b2)*zs/(5*rs5) + 3*zx2s/rs5 - (70*a2 - 42*b2)*x*xs3*zs/rs9 + (30*a2 - 18*b2)*x*xzs/rs7 + 6*x2*zs/rs5; %%%% checked
            K2star123 = 6*x*yzs/rs5 + (10*a2 - 6*b2)*x*yzs/rs7 + 3*xyzs/rs5 - 42*a2b2*xyzs/rs9 + (15*a2 - 3*b2)*xyzs/rs7 - 70*a2*x3yzs/rs9 + 126*a2b2*x3yzs/rs11 - (70*a2 - 42*b2)*x*x2yzs/rs9 - 30*x2*xyzs/rs7 - 30*x*x2yzs/rs7; %%%% checked
            
            K2star212 = (10*a2 - 6*b2)*x*zs2/rs7 - 6*x2*xs/rs5 - 6*x*xs2/rs5 - (40*a2 - 24*b2)*x/(5*rs5) + 70*a2*x3y2s/rs9 - 126*a2b2*x3y2s/rs11 + 30*x2*xy2s/rs7 + 30*x*x2y2s/rs7 + (70*a2 - 42*b2)*x*x2y2s/rs9 - (15*a2 + 3*b2)*xy2s/rs7 + 3*xy2s/rs5 + 42*a2b2*xy2s/rs9 + (15*a2 + 3*b2)*xs/(5*rs5) - 6*a2b2*xs/rs7 - 10*a2*xs3/rs7 + 14*a2b2*xs3/rs9; %%%% checked
            K2star213 = 30*x2*xyzs/rs7 + 30*x*x2yzs/rs7 + (70*a2 - 42*b2)*x*x2yzs/rs9 - (10*a2 - 6*b2)*x*yzs/rs7 + 70*a2*x3yzs/rs9 - 126*a2b2*x3yzs/rs11 + 3*xyzs/rs5 - (15*a2 + 3*b2)*xyzs/rs7 + 42*a2b2*xyzs/rs9; %%%% checked
            K2star223 = (5*a2 + 3*b2)*zs/(5*rs5) + (70*a2 - 42*b2)*x*ys2*xzs/rs9 + 3*zy2s/rs5 - (5*a2 + 3*b2)*zy2s/rs7 - 126*a2b2*zx2s*ys2/rs11 + 70*a2*zx2s*ys2/rs9 - 6*x2*zs/rs5 + 30*x*ys2*xzs/rs7 - (10*a2 - 6*b2)*x*xzs/rs7 + 30*x2*zy2s/rs7 - 6*x*xzs/rs5 + 12*a2b2*zs/rs7 - 10*a2*zx2s/rs7 - 14*a2b2*zs3/rs9; %%%% checked
            
            K2star313 = (10*a2 - 6*b2)*x*ys2/rs7 - 6*x2*xs/rs5 - 6*x*xs2/rs5 - (40*a2 - 24*b2)*x/(5*rs5) + 70*a2*x3z2s/rs9 - 126*a2b2*x3z2s/rs11 + 30*x2*xz2s/rs7 + 30*x*x2z2s/rs7 + (70*a2 - 42*b2)*x*x2z2s/rs9 - (15*a2 + 3*b2)*xz2s/rs7 + 3*xz2s/rs5 + 42*a2b2*xz2s/rs9 + (15*a2 + 3*b2)*xs/(5*rs5) - 6*a2b2*xs/rs7 - 10*a2*xs3/rs7 + 14*a2b2*xs3/rs9 ; %%%% checked
            K2star323 = (5*a2 + 3*b2)*ys/(5*rs5) + (70*a2 - 42*b2)*x*xys*zs2/rs9 + 3*yz2s/rs5 - (5*a2 + 3*b2)*yz2s/rs7 - 126*a2b2*yx2s*zs2/rs11 + 70*a2*yx2s*zs2/rs9 - 6*x2*ys/rs5 + 30*x*xys*zs2/rs7 - (10*a2 - 6*b2)*x*xys/rs7 + 30*x2*yz2s/rs7 - 6*x*xys/rs5 + 12*a2b2*ys/rs7 - 10*a2*yx2s/rs7 - 14*a2b2*ys3/rs9; %%%% checked
        
            res = [30*a2b2*xs/rs7 - (45*a2 - 9*b2)*xs/(5*rs5) + (65*a2 - 3*b2)*xs3/rs7 - xs/rs3 - 140*a2b2*xs3/rs9 + 3*xs3/rs5 - 70*a2*xs5/rs9 + 126*a2b2*xs5/rs11 - (70*a2 - 42*b2)*x*xs4/rs9 - 30*x2*xs3/rs7 + 18*x2*xs/rs5 + 24*x*xs2/rs5 - 30*x*xs4/rs7 - (30*a2 - 18*b2)*x/(5*rs5) - 2*x/rs3 + (60*a2 - 36*b2)*x*xs2/rs7, ... %%%% checked
                -30*x2*xy2s/rs7 - 30*x*x2y2s/rs7 + 6*a2b2*xs/rs7 - 6*x*zs2/rs5 + 126*a2b2*x3y2s/rs11 - 70*a2*x3y2s/rs9 - xs/rs3 - 14*a2b2*xs3/rs9 + 10*a2*xs3/rs7 + 3*xy2s/rs5 + 6*x2*xs/rs5 + 4*x/rs3 - (10*a2 - 6*b2)*x*zs2/rs7 + (40*a2 - 24*b2)*x/(5*rs5) - (15*a2 - 3*b2)*xs/(5*rs5) - (70*a2 - 42*b2)*x*x2y2s/rs9 - 42*a2b2*xy2s/rs9 + (15*a2 - 3*b2)*xy2s/rs7,... %%%% checked
                -30*x2*xz2s/rs7 - 30*x*x2z2s/rs7 + 6*a2b2*xs/rs7 - 6*x*ys2/rs5 + 126*a2b2*x3z2s/rs11 - 70*a2*x3z2s/rs9 - xs/rs3 - 14*a2b2*xs3/rs9 + 10*a2*xs3/rs7 + 3*xz2s/rs5 + 6*x2*xs/rs5 + 4*x/rs3 - (10*a2 - 6*b2)*x*ys2/rs7 + (40*a2 - 24*b2)*x/(5*rs5) - (15*a2 - 3*b2)*xs/(5*rs5) - (70*a2 - 42*b2)*x*x2z2s/rs9 - 42*a2b2*xz2s/rs9 + (15*a2 - 3*b2)*xz2s/rs7,... %%%% checked
                K2star112,... %%% see above
                K2star113,... %%% see above
                K2star123,... %%% see above
                K2star112,... %%% sym checked
                K2star113,... %%% sym checked
                K2star123;... %%% sym checked
                % 2nd line
                30*x2*yx2s/rs7 - 6*x*xys/rs5 + 30*x*xs3*ys/rs7 - ys/rs3 - 126*a2b2*ys*xs4/rs11 + 70*a2*ys*xs4/rs9 - 6*a2b2*ys/rs7 + (5*a2 + 3*b2)*ys/(5*rs5) + 84*a2b2*yx2s/rs9 - (35*a2 + 3*b2)*yx2s/rs7 + 3*yx2s/rs5 + (70*a2 - 42*b2)*x*xs3*ys/rs9 - (30*a2 - 18*b2)*x*xys/rs7 - 6*x2*ys/rs5,... %%%% checked
                14*a2b2*ys3/rs9 + 30*x2*ys3/rs7 - 18*x2*ys/rs5 + (70*a2 - 42*b2)*x*ys3*xs/rs9 - 18*x*xys/rs5 + 42*a2b2*yx2s/rs9 - 126*a2b2*ys3*xs2/rs11 - 30*a2*yx2s/rs7 + 70*a2*ys3*xs2/rs9 - ys/rs3 - (5*a2 + 3*b2)*ys3/rs7 + 3*ys3/rs5 + 30*x*ys3*xs/rs7 - (30*a2 - 18*b2)*x*xys/rs7 + (15*a2 + 9*b2)*ys/(5*rs5) - 6*a2b2*ys/rs7,... %%%% checked 
                30*x*xys*zs2/rs7 - (5*a2 + 3*b2)*yz2s/rs7 + 14*a2b2*ys3/rs9 + 3*yz2s/rs5 + 70*a2*yx2s*zs2/rs9 - 126*a2b2*yx2s*zs2/rs11 - 6*x2*ys/rs5 - (10*a2 - 6*b2)*x*xys/rs7 - 6*x*xys/rs5 + 12*a2b2*ys/rs7 - 10*a2*yx2s/rs7 - ys/rs3 + 30*x2*yz2s/rs7 + (70*a2 - 42*b2)*x*xys*zs2/rs9 + (5*a2 + 3*b2)*ys/(5*rs5),... %%%% checked
                K2star212,... %%% see above
                K2star213,... %%% see above
                K2star223,... %%% see above
                K2star212,... %%% sym checked
                K2star213,... %%% sym checked
                K2star223;... %%% sym checked
                % 3rd line
                30*x2*zx2s/rs7 - 6*x*xzs/rs5 + 30*x*xs3*zs/rs7 - zs/rs3 - 126*a2b2*zs*xs4/rs11 + 70*a2*zs*xs4/rs9 - 6*a2b2*zs/rs7 + (5*a2 + 3*b2)*zs/(5*rs5) + 84*a2b2*zx2s/rs9 - (35*a2 + 3*b2)*zx2s/rs7 + 3*zx2s/rs5 + (70*a2 - 42*b2)*x*xs3*zs/rs9 - (30*a2 - 18*b2)*x*xzs/rs7 - 6*x2*zs/rs5,... %%%% checked
                30*x*ys2*xzs/rs7 - (5*a2 + 3*b2)*zy2s/rs7 + 14*a2b2*zs3/rs9 + 3*zy2s/rs5 + 70*a2*zx2s*ys2/rs9 - 126*a2b2*zx2s*ys2/rs11 - 6*x2*zs/rs5 - (10*a2 - 6*b2)*x*xzs/rs7 - 6*x*xzs/rs5 + 12*a2b2*zs/rs7 - 10*a2*zx2s/rs7 - zs/rs3 + 30*x2*zy2s/rs7 + (70*a2 - 42*b2)*x*ys2*xzs/rs9 + (5*a2 + 3*b2)*zs/(5*rs5),... %%%% checked
                14*a2b2*zs3/rs9 + 30*x2*zs3/rs7 - 18*x2*zs/rs5 + (70*a2 - 42*b2)*x*zs3*xs/rs9 - 18*x*xzs/rs5 + 42*a2b2*zx2s/rs9 - 126*a2b2*zs3*xs2/rs11 - 30*a2*zx2s/rs7 + 70*a2*zs3*xs2/rs9 - zs/rs3 - (5*a2 + 3*b2)*zs3/rs7 + 3*zs3/rs5 + 30*x*zs3*xs/rs7 - (30*a2 - 18*b2)*x*xzs/rs7 + (15*a2 + 9*b2)*zs/(5*rs5) - 6*a2b2*zs/rs7,... %%%% checked
                K2star213,... %%% sym checked
                K2star313,... %%% see above
                K2star323,... %%% see above
                K2star213,... %%% sym checked
                K2star313,... %%% sym checked
                K2star323];   %%% sym checked
        end
        function [res] = tens_K3star39(x1,x_alpha, beta) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
        
            b2 = beta^2;
            rs5 = rs^5; rs7 = rs^7; rs9 = rs^9; 
            rs5t2 = 2*rs5;
            
            xs2 = xs^2; ys2 = ys^2; zs2 = zs^2;
            xs3 = xs^3; ys3 = ys^3; zs3 = zs^3;
            
            xys = xs*ys; %xy2s = xs*ys2; yx2s = ys*xs2; 
            xzs = xs*zs; %xz2s = xs*zs2; zx2s = zs*xs2;
            yzs = ys*zs; %yz2s = ys*zs2; zy2s = zs*ys2;
        
            xyzs = xs*ys*zs; 
            x2yzs = xs2*yzs; xy2zs = ys2*xzs; xyz2s = xys*zs2;
            x2y2s = xs2*ys2; x2z2s = xs2*zs2; 
            
            K3star112 = (3*xzs)/(rs5t2); %
            K3star113 = -(3*xys)/(rs5t2); %
            K3star123 = (3*zs2)/(rs5t2) - (3*ys2)/(rs5t2); %
            
            K3star212 = -30*x2yzs/rs7 + 42*b2*x2yzs/rs9 - 6*b2*yzs/rs7 - 30*x*xyzs/rs7 + (3*yzs)/(rs5t2); %%%%
            K3star213 = -24*b2/(5*rs5) + (9*xs2)/(rs5t2) + 6*xs*x/rs5 + 6*b2*ys2/rs7 + (3*zs2)/(rs5t2) - 30*xs*x*zs2/rs7 + 42*b2*x2z2s/rs9 - 30*x2z2s/rs7; %%%%
            K3star223 = (9*xys)/(rs5t2) - 30*xyz2s/rs7 - 6*b2*xys/rs7 + 42*b2*xyz2s/rs9 + 6*x*ys/rs5 - 30*x*ys*zs2/rs7; %%%%
            
            K3star312 = 24*b2/(5*rs5) - (9*xs2)/(rs5t2) - 6*xs*x/rs5 - 6*b2*zs2/rs7 - (3*ys2)/(rs5t2) + 30*xs*x*ys2/rs7 - 42*b2*x2y2s/rs9 + 30*x2y2s/rs7; %%%%
            K3star323 = 30*xy2zs/rs7 - 42*b2*xy2zs/rs9 + 6*b2*xzs/rs7 + 30*x*zs*ys2/rs7 - 6*x*zs/rs5 - (9*xzs)/(rs5t2); %%%%
        
            res = [0,... %
                3*yzs/rs5,... %
                -3*yzs/rs5,... %
                K3star112,... %
                K3star113,... %
                K3star123,... %
                K3star112,... %%% sym checked
                K3star113,... %%% sym checked
                K3star123;... %%% sym checked
                % 2nd line
                9*xzs/rs5 - 30*xs3*zs/rs7 + 42*b2*xs3*zs/rs9 - 18*b2*xzs/rs7 + 6*x*zs/rs5 - 30*x*zs*xs2/rs7,... %%%%
                6*xzs/rs5 - 30*xy2zs/rs7 + 42*b2*xy2zs/rs9 - 6*b2*xzs/rs7 + 6*x*zs/rs5 - 30*x*zs*ys2/rs7,... %%%%
                15*xzs/rs5 - 30*zs3*xs/rs7 - 18*b2*xzs/rs7 + 42*b2*zs3*xs/rs9 + 18*x*zs/rs5 - 30*x*zs3/rs7,... %%%%
                K3star212,... %
                K3star213,... %
                K3star223,... %
                K3star212,... %%% sym checked
                K3star213,... %%% sym checked
                K3star223;... %%% sym checked
                % 3rd line
                -9*xys/rs5 + 30*xs3*ys/rs7 - 42*b2*xs3*ys/rs9 + 18*b2*xys/rs7 - 6*x*ys/rs5 + 30*x*ys*xs2/rs7,... %%%%
                -15*xys/rs5 + 30*ys3*xs/rs7 + 18*b2*xys/rs7 - 42*b2*ys3*xs/rs9 - 18*x*ys/rs5 + 30*x*ys3/rs7,... %%%%
                -6*xys/rs5 + 30*xyz2s/rs7 - 42*b2*xyz2s/rs9 + 6*b2*xys/rs7 - 6*x*ys/rs5 + 30*x*ys*zs2/rs7,... %%%%
                K3star312,... %
                -K3star212,... %%% sym checked
                K3star323,... %
                K3star312,... %%% sym checked
                -K3star212,... %%% sym checked
                K3star323]; %%% sym checked
        end
        function [res] = tens_K4star99(x1,x_alpha, alpha, beta) 
            x=x1(1); y=x1(2); z=x1(3);                              % coordinates of a point in space
            x_a = x_alpha(1); y_a = x_alpha(2); z_a = x_alpha(3);   % coordinates of the center of the sphere alpha
            rs = sqrt((x + x_a)^2 + (y - y_a)^2 + (z - z_a)^2);     % star norm
            xs = -(x + x_a); ys = y - y_a; zs = z - z_a;            % star components
            
            a2 = alpha^2; b2 = beta^2; a2b2 = a2*b2;
            rs3 = rs^3; rs5 = rs^5; rs7 = rs^7; rs9 = rs^9; rs11 = rs^11;
            rs5t2 = 2*rs5;
            rs5t5 = 5*rs5; rs7t5 = 5*rs7; rs9t5 = 5*rs9; rs11t5 = 5*rs11; rs13t5 = 5*rs^13;
        
            x2 = x^2;
            xs2 = xs^2; ys2 = ys^2; zs2 = zs^2;
            xs3 = xs^3; ys3 = ys^3; zs3 = zs^3;
            xs4 = xs^4; ys4 = ys^4; zs4 = zs^4;
            xs5 = xs^5;
            xs6 = xs^6;
            
            xys = xs*ys; xy2s = xs*ys2; yx2s = ys*xs2; 
            xzs = xs*zs; xz2s = xs*zs2; zx2s = zs*xs2;
            yzs = ys*zs; %yz2s = ys*zs2; %zy2s = zs*ys2;
        
            xyzs = xs*ys*zs; 
            x2yzs = xs2*yzs; xy2zs = ys2*xzs; xyz2s = xys*zs2;
            x3yzs = xs3*yzs;
            x2y2s = xs2*ys2; x3y2s = xs3*ys2;
            x2z2s = xs2*zs2; x3z2s = xs3*zs2;
            y2z2s = ys2*zs2;
        
            x2y2z2s = xs2*ys2*zs2;
            x4y2s = xs4*ys2; x4z2s = xs4*zs2;
            x3y2zs = xs3*zs*ys2;
            x3z2ys = xs3*ys*zs2;
            
            K4star1122 = -(21*a2 + 3*b2)*xs2/rs7 + (252*a2b2*xs2)/rs9t5 + 2/rs3 + 24*x2/rs5 + 42*a2*xs4/rs9 - (378*a2b2*xs4)/rs11t5 + 3*(a2 + b2)/rs5t5 - 18*a2b2/rs7t5 - 378*a2*x4y2s/rs11 + (4158*a2b2*x4y2s)/rs13t5 - 30*x2*zs2/rs7 + 24*x*xs/rs5 + (126*a2b2*ys2)/rs9t5 - 3*(a2 + b2)*ys2/rs7 - 30*x*xz2s/rs7 - 210*x2*x2y2s/rs9 - 210*x*x3y2s/rs9 - 15*x2y2s/rs7 + (147*a2 + 21*b2)*x2y2s/rs9 - (2268*a2b2*x2y2s)/rs11t5 - 3*zs2/rs5 - 378*(a2 - b2)*x*x3y2s/rs11 + 126*(a2 - b2)*x*xy2s/rs9 + 42*(a2 - b2)*x*xs3/rs9 - 18*(a2 - b2)*x*xs/rs7; %%%% checked
            K4star1133 = -(21*a2 + 3*b2)*xs2/rs7 + (252*a2b2*xs2)/rs9t5 + 2/rs3 + 24*x2/rs5 + 42*a2*xs4/rs9 - (378*a2b2*xs4)/rs11t5 + 3*(a2 + b2)/rs5t5 - 18*a2b2/rs7t5 - 378*a2*x4z2s/rs11 + (4158*a2b2*x4z2s)/rs13t5 - 30*x2*ys2/rs7 + 24*x*xs/rs5 + (126*a2b2*zs2)/rs9t5 - 3*(a2 + b2)*zs2/rs7 - 30*x*xy2s/rs7 - 210*x2*x2z2s/rs9 - 210*x*x3z2s/rs9 - 15*x2z2s/rs7 + (147*a2 + 21*b2)*x2z2s/rs9 - (2268*a2b2*x2z2s)/rs11t5 - 3*ys2/rs5 - 378*(a2 - b2)*x*x3z2s/rs11 + 126*(a2 - b2)*x*xz2s/rs9 + 42*(a2 - b2)*x*xs3/rs9 - 18*(a2 - b2)*x*xs/rs7; %%%% checked
            K4star1112 = 90*x2*xys/rs7 - 210*x2*xs3*ys/rs9 - 210*x*ys*xs4/rs9 - (27*a2 + 9*b2)*xys/rs7 + 126*a2b2*xys/rs9 + 6*xys/rs5 + (273*a2 + 21*b2)*xs3*ys/rs9 - 756*a2b2*xs3*ys/rs11 - 15*xs3*ys/rs7 - 378*a2*xs5*ys/rs11 + (4158*a2b2*xs5*ys)/rs13t5 + 90*x*yx2s/rs7 - 378*(a2 - b2)*x*ys*xs4/rs11 + 252*(a2 - b2)*x*ys*xs2/rs9 - 18*(a2 - b2)*x*ys/rs7; %%%% checked 
            K4star1113 = 90*x2*xzs/rs7 - 210*x2*xs3*zs/rs9 - 210*x*zs*xs4/rs9 - (27*a2 + 9*b2)*xzs/rs7 + 126*a2b2*xzs/rs9 + 6*xzs/rs5 + (273*a2 + 21*b2)*xs3*zs/rs9 - 756*a2b2*xs3*zs/rs11 - 15*xs3*zs/rs7 - 378*a2*xs5*zs/rs11 + (4158*a2b2*xs5*zs)/rs13t5 + 90*x*zx2s/rs7 - 378*(a2 - b2)*x*zs*xs4/rs11 + 252*(a2 - b2)*x*zs*xs2/rs9 - 18*(a2 - b2)*x*zs/rs7; %%%% checked
            K4star1123 = 3*yzs/rs5 - 3*(a2 + b2)*yzs/rs7 + (126*a2b2*yzs)/rs9t5 - 210*x2*x2yzs/rs9 - 210*x*x3yzs/rs9 - 15*x2yzs/rs7 + (147*a2 + 21*b2)*x2yzs/rs9 - (2268*a2b2*x2yzs)/rs11t5 - 378*a2*xs4*yzs/rs11 + (4158*a2b2*xs4*yzs)/rs13t5 + 30*x2*yzs/rs7 + 30*x*xyzs/rs7 + 126*(a2 - b2)*x*xyzs/rs9 - 378*(a2 - b2)*x*x3yzs/rs11; %%%% checked
            
            K4star2233 = -210*x2*y2z2s/rs9 - 378*a2*x2y2z2s/rs11 - 378*a2b2*y2z2s/rs11t5 - 30*x*xs3/rs7 - 42*a2*xs4/rs9 + 36*a2b2/rs7t5 - 30*x2*xs2/rs7 + 24*x2/rs5 - 3*xs2/rs5 - 210*x*xs*y2z2s/rs9 + 4158*a2b2*x2y2z2s/rs13t5 + 2/rs3 - 15*y2z2s/rs7 + 21*(a2 + b2)*y2z2s/rs9 + 24*x*xs/rs5 + 378*a2b2*xs4/rs11t5 - 378*a2b2*xs2/rs9t5 - 12*(a2 + b2)/rs5t5 + (39*a2 + 3*b2)*xs2/rs7 - 378*(a2 - b2)*x*xs*y2z2s/rs11 - 42*(a2 - b2)*x*xs3/rs9 + 36*(a2 - b2)*x*xs/rs7; %%%% checked
            K4star2212 = 90*x2*xys/rs7 + (1134*a2b2*xyz2s)/rs11t5 - (27*a2 + 9*b2)*xys/rs7 + (63*a2 + 21*b2)*ys3*xs/rs9 + 6*xys/rs5 + 126*a2*xs3*ys/rs9 - 378*a2*xs3*ys3/rs11 - (756*a2b2*xys)/rs9t5 + (4158*a2b2*xs3*ys3)/rs13t5 - 15*ys3*xs/rs7 + 90*x*yx2s/rs7 - 210*x2*ys3*xs/rs9 - 210*x*ys3*xs2/rs9 - 378*(a2 - b2)*x*ys3*xs2/rs11 + 126*(a2 - b2)*x*yx2s/rs9 + 42*(a2 - b2)*x*ys3/rs9 - 18*(a2 - b2)*x*ys/rs7; %%%% checked
            K4star2213 = 30*x2*xzs/rs7 + 30*x*zx2s/rs7 + 42*a2*xs3*zs/rs9 - (378*a2b2*xs3*zs)/rs11t5 - (9*a2 + 3*b2)*xzs/rs7 + (126*a2b2*xzs)/rs9t5 + 3*xzs/rs5 + (63*a2 + 21*b2)*xy2zs/rs9 - (1134*a2b2*xy2zs)/rs11t5 - 15*xy2zs/rs7 - 378*a2*x3y2zs/rs11 + (4158*a2b2*x3y2zs)/rs13t5 - 210*x2*xy2zs/rs9 - 210*x*zs*x2y2s/rs9 - 378*(a2 - b2)*x*zs*x2y2s/rs11 - 42*(a2 - b2)*x*zs3/rs9 + 36*(a2 - b2)*x*zs/rs7; %%%% checked
            K4star2223 = 90*x2*yzs/rs7 + 126*a2*x2yzs/rs9 - (1134*a2b2*x2yzs)/rs11t5 - 378*a2*xs2*ys3*zs/rs11 + (4158*a2b2*xs2*ys3*zs)/rs13t5 + 6*yzs/rs5 - (9*a2 + 9*b2)*yzs/rs7 + (126*a2b2*yzs)/rs9t5 + (21*a2 + 21*b2)*ys3*zs/rs9 - (378*a2b2*ys3*zs)/rs11t5 - 210*x2*ys3*zs/rs9 + 90*x*xyzs/rs7 - 210*x*ys3*xzs/rs9 - 15*ys3*zs/rs7 + 126*(a2 - b2)*x*xyzs/rs9 - 378*(a2 - b2)*x*xzs*ys3/rs11; %%%% checked
            
            K4star3312 = 30*x2*xys/rs7 + 30*x*yx2s/rs7 + 42*a2*xs3*ys/rs9 - (378*a2b2*xs3*ys)/rs11t5 - (9*a2 + 3*b2)*xys/rs7 + (126*a2b2*xys)/rs9t5 + 3*xys/rs5 + (63*a2 + 21*b2)*xyz2s/rs9 - (1134*a2b2*xyz2s)/rs11t5 - 15*xyz2s/rs7 - 378*a2*x3z2ys/rs11 + (4158*a2b2*x3z2ys)/rs13t5 - 210*x2*xyz2s/rs9 - 210*x*ys*x2z2s/rs9 - 378*(a2 - b2)*x*ys*x2z2s/rs11 - 42*(a2 - b2)*x*ys3/rs9 + 36*(a2 - b2)*x*ys/rs7; %%%% checked
            K4star3313 = 90*x2*xzs/rs7 + (1134*a2b2*xy2zs)/rs11t5 - (27*a2 + 9*b2)*xzs/rs7 + (63*a2 + 21*b2)*zs3*xs/rs9 + 6*xzs/rs5 + 126*a2*xs3*zs/rs9 - 378*a2*xs3*zs3/rs11 - (756*a2b2*xzs)/rs9t5 + (4158*a2b2*xs3*zs3)/rs13t5 - 15*zs3*xs/rs7 + 90*x*zx2s/rs7 - 210*x2*zs3*xs/rs9 - 210*x*zs3*xs2/rs9 - 378*(a2 - b2)*x*zs3*xs2/rs11 + 126*(a2 - b2)*x*zx2s/rs9 + 42*(a2 - b2)*x*zs3/rs9 - 18*(a2 - b2)*x*zs/rs7; %%%% checked
            K4star3323 = 90*x2*yzs/rs7 + 126*a2*x2yzs/rs9 - (1134*a2b2*x2yzs)/rs11t5 - 378*a2*xs2*ys*zs3/rs11 + (4158*a2b2*xs2*ys*zs3)/rs13t5 + 6*yzs/rs5 - (9*a2 + 9*b2)*yzs/rs7 + (126*a2b2*yzs)/rs9t5 + (21*a2 + 21*b2)*ys*zs3/rs9 - (378*a2b2*ys*zs3)/rs11t5 - 210*x2*ys*zs3/rs9 + 90*x*xyzs/rs7 - 210*x*xys*zs3/rs9 - 15*ys*zs3/rs7 + 126*(a2 - b2)*x*xyzs/rs9 - 378*(a2 - b2)*x*xys*zs3/rs11; %%%% checked
            
            K4star1212 = -(3*ys2)/rs5t2 - 24*x2/rs5 + 378*a2*x4y2s/rs11 - (4158*a2b2*x4y2s)/rs13t5 - (147*a2 + 21*b2)*x2y2s/rs9 + (2268*a2b2*x2y2s)/rs11t5 + 15*x2y2s/rs7 - 24*xs*x/rs5 - 42*a2*xs4/rs9 - (3*xs2)/rs5t2 + (21*a2 + 3*b2)*xs2/rs7 - (252*a2b2*xs2)/rs9t5 + 210*x2*x2y2s/rs9 + 210*x*x3y2s/rs9 + (378*a2b2*xs4)/rs11t5 + 30*x2*zs2/rs7 + (3*a2 + 3*b2)*ys2/rs7 - (126*a2b2*ys2)/rs9t5 - (3*a2 + 3*b2)/rs5t5 + 18*a2b2/rs7t5 + 30*x*xz2s/rs7 - 42*(a2 - b2)*x*xs3/rs9 + 18*(a2 - b2)*x*xs/rs7 + 378*(a2 - b2)*x*x3y2s/rs11 - 126*(a2 - b2)*x*xy2s/rs9; %%%% checked
            K4star1213 = -(3*yzs)/rs5t2 + 3*(a2 + b2)*yzs/rs7 - (126*a2b2*yzs)/rs9t5 + 210*x2*x2yzs/rs9 + 210*x*x3yzs/rs9 + 15*x2yzs/rs7 - (147*a2 + 21*b2)*x2yzs/rs9 + (2268*a2b2*x2yzs)/rs11t5 + 378*a2*xs4*yzs/rs11 - (4158*a2b2*xs4*yzs)/rs13t5 - 30*x2*yzs/rs7 - 30*x*xyzs/rs7 + 378*(a2 - b2)*x*x3yzs/rs11 - 126*(a2 - b2)*x*xyzs/rs9; %%%% checked
            K4star1223 = 378*a2*x3y2zs/rs11 - (4158*a2b2*x3y2zs)/rs13t5 + (378*a2b2*xs3*zs)/rs11t5 + (9*a2 + 3*b2)*xzs/rs7 - (126*a2b2*xzs)/rs9t5 - (3*xzs)/rs5t2 - 42*a2*xs3*zs/rs9 - (63*a2 + 21*b2)*xy2zs/rs9 + (1134*a2b2*xy2zs)/rs11t5 + 15*xy2zs/rs7 - 30*x*zx2s/rs7 + 210*x2*xy2zs/rs9 + 210*x*zs*x2y2s/rs9 - 30*x2*xzs/rs7 + 378*(a2 - b2)*x*zs*x2y2s/rs11 + 42*(a2 - b2)*x*zs3/rs9 - 36*(a2 - b2)*x*zs/rs7; %%%% checked
            
            K4star1313 = -(3*zs2)/rs5t2 - 24*x2/rs5 + 378*a2*x4z2s/rs11 - (4158*a2b2*x4z2s)/rs13t5 - (147*a2 + 21*b2)*x2z2s/rs9 + (2268*a2b2*x2z2s)/rs11t5 + 15*x2z2s/rs7 - 24*xs*x/rs5 - 42*a2*xs4/rs9 - (3*xs2)/rs5t2 + (21*a2 + 3*b2)*xs2/rs7 - (252*a2b2*xs2)/rs9t5 + 210*x2*x2z2s/rs9 + 210*x*x3z2s/rs9 + (378*a2b2*xs4)/rs11t5 + 30*x2*ys2/rs7 + (3*a2 + 3*b2)*zs2/rs7 - (126*a2b2*zs2)/rs9t5 - (3*a2 + 3*b2)/rs5t5 + 18*a2b2/rs7t5 + 30*x*xy2s/rs7 - 42*(a2 - b2)*x*xs3/rs9 + 18*(a2 - b2)*x*xs/rs7 + 378*(a2 - b2)*x*x3z2s/rs11 - 126*(a2 - b2)*x*xz2s/rs9; %%%% checked
            K4star1323 = 378*a2*x3z2ys/rs11 - (4158*a2b2*x3z2ys)/rs13t5 + (378*a2b2*xs3*ys)/rs11t5 + (9*a2 + 3*b2)*xys/rs7 - (126*a2b2*xys)/rs9t5 - (3*xys)/rs5t2 - 42*a2*xs3*ys/rs9 - (63*a2 + 21*b2)*xyz2s/rs9 + (1134*a2b2*xyz2s)/rs11t5 + 15*xyz2s/rs7 - 30*x*yx2s/rs7 + 210*x2*xyz2s/rs9 + 210*x*ys*x2z2s/rs9 - 30*x2*xys/rs7 + 378*(a2 - b2)*x*ys*x2z2s/rs11 + 42*(a2 - b2)*x*ys3/rs9 - 36*(a2 - b2)*x*ys/rs7; %%%% checked
            
            K4star2323 = 24*x2/rs5 - 30*x2*xs2/rs7 + 42*a2*x2z2s/rs9 - (378*a2b2*x2z2s)/rs11t5 + 42*a2*x2y2s/rs9 - (378*a2b2*x2y2s)/rs11t5 + 21*(a2 + b2)*y2z2s/rs9 - (378*a2b2*y2z2s)/rs11t5 + (3*zs2)/rs5t2 + (3*ys2)/rs5t2 + 36*a2b2/rs7t5 + 24*x*xs/rs5 - 12*(a2 + b2)/rs5t5 - 210*x*xs*y2z2s/rs9 - 30*x*xs3/rs7 - 378*a2*x2y2z2s/rs11 + (4158*a2b2*x2y2z2s)/rs13t5 - 15*y2z2s/rs7 - 210*x2*y2z2s/rs9 - 378*(a2 - b2)*x*xs*y2z2s/rs11 - 42*(a2 - b2)*x*xs3/rs9 + 36*(a2 - b2)*x*xs/rs7 - 3*(a2 - b2)*xs2/rs7; %%%% checked
            
            res = [-(108*a2 + 18*b2)*xs2/rs7 + 378*a2b2*xs2/rs9 + 12*xs2/rs5 - 18*x2/rs5 + 180*x2*xs2/rs7 - 210*x2*xs4/rs9 - 210*x*xs5/rs9 + (441*a2 + 21*b2)*xs4/rs9 - 1134*a2b2*xs4/rs11 - 15*xs4/rs7 + 9*(a2 + b2)/rs5t5 - 18*a2b2/rs7 + 180*x*xs3/rs7 - 18*xs*x/rs5 - 1/rs3 - 378*a2*xs6/rs11 + (4158*a2b2*xs6)/rs13t5 - 378*(a2 - b2)*x*xs5/rs11 + 420*(a2 - b2)*x*xs3/rs9 - 90*(a2 - b2)*x*xs/rs7, ... %%%% checked
                K4star1122,... %%% see above
                K4star1133,... %%% see above
                K4star1112,... %%% see above
                K4star1113,... %%% see above
                K4star1123,... %%% see above
                K4star1112,... %%% sym checked
                K4star1113,... %%% sym checked
                K4star1123;... %%% sym checked
                % 2nd line
                K4star1122,... %%% sym checked
                -210*xs*x*ys4/rs9 - 18*xs*x/rs5 - (378*a2b2*ys4)/rs11t5 + 21*(a2 + b2)*ys4/rs9 + 12*ys2/rs5 + (252*a2b2*ys2)/rs9t5 - 18*(a2 + b2)*ys2/rs7 + 180*x*xy2s/rs7 + 9*(a2 + b2)/rs5t5 - 18*a2b2/rs7t5 - 18*a2*xs2/rs7 + (126*a2b2*xs2)/rs9t5 + 180*x2*ys2/rs7 - 15*ys4/rs7 + 252*a2*x2y2s/rs9 - (2268*a2b2*x2y2s)/rs11t5 - 378*a2*xs2*ys4/rs11 + (4158*a2b2*xs2*ys4)/rs13t5 - 1/rs3 - 18*x2/rs5 - 210*x2*ys4/rs9 + 252*(a2 - b2)*x*xy2s/rs9 - 378*(a2 - b2)*x*ys4*xs/rs11 - 18*(a2 - b2)*x*xs/rs7,... %%%% checked
                K4star2233,... %%% see above
                K4star2212,... %%% see above
                K4star2213,... %%% see above
                K4star2223,... %%% see above
                K4star2212,... %%% sym checked
                K4star2213,... %%% sym checked
                K4star2223;... %%% sym checked
                % 3rd line
                K4star1133,... %%% sym checked
                K4star2233,... %%% sym checked
                -210*xs*x*zs4/rs9 - 18*xs*x/rs5 - (378*a2b2*zs4)/rs11t5 + 21*(a2 + b2)*zs4/rs9 + 12*zs2/rs5 + (252*a2b2*zs2)/rs9t5 - 18*(a2 + b2)*zs2/rs7 + 180*x*xz2s/rs7 + 9*(a2 + b2)/rs5t5 - 18*a2b2/rs7t5 - 18*a2*xs2/rs7 + (126*a2b2*xs2)/rs9t5 + 180*x2*zs2/rs7 - 15*zs4/rs7 + 252*a2*x2z2s/rs9 - (2268*a2b2*x2z2s)/rs11t5 - 378*a2*xs2*zs4/rs11 + (4158*a2b2*xs2*zs4)/rs13t5 - 1/rs3 - 18*x2/rs5 - 210*x2*zs4/rs9 + 252*(a2 - b2)*x*xz2s/rs9 - 378*(a2 - b2)*x*zs4*xs/rs11 - 18*(a2 - b2)*x*xs/rs7,... %%%% checked
                K4star3312,... %%% see above
                K4star3313,... %%% see above
                K4star3323,... %%% see above
                K4star3312,... %%% sym checked
                K4star3313,... %%% sym checked
                K4star3323;... %%% sym checked
                % 4th line
                -K4star1112,... %%% sym checked
                -K4star2212,... %%% sym checked
                -K4star3312,... %%% sym checked
                K4star1212,... %%% see above
                K4star1213,... %%% see above
                K4star1223,... %%% see above
                K4star1212,... %%% sym checked
                K4star1213,... %%% sym checked
                K4star1223;... %%% sym checked
                % 5th line
                -K4star1113,... %%% sym checked
                -K4star2213,... %%% sym checked
                -K4star3313,... %%% sym checked
                K4star1213,... %%% sym checked
                K4star1313,... %%% see above
                K4star1323,... %%% see above
                K4star1213,... %%% sym checked
                K4star1313,... %%% sym checked
                K4star1323;... %%% sym checked
                % 6th line
                K4star1123,... %%% sym checked
                K4star2223,... %%% sym checked
                K4star3323,... %%% sym checked
                -K4star1223,... %%% sym checked
                -K4star1323,... %%% sym checked
                K4star2323,... %%% see above
                -K4star1223,... %%% sym checked
                -K4star1323,... %%% sym checked
                K4star2323;... %%% sym checked
                % 7th line
                -K4star1112,... %%% sym checked
                -K4star2212,... %%% sym checked
                -K4star3312,... %%% sym checked
                K4star1212,... %%% sym checked
                K4star1213,... %%% sym checked
                K4star1223,... %%% sym checked
                K4star1212,... %%% sym checked
                K4star1213,... %%% sym checked
                K4star1223;... %%% sym checked
                % 8th line
                -K4star1113,... %%% sym checked
                -K4star2213,... %%% sym checked
                -K4star3313,... %%% sym checked
                K4star1213,... %%% sym checked
                K4star1313,... %%% sym checked
                K4star1323,... %%% sym checked
                K4star1213,... %%% sym checked
                K4star1313,... %%% sym checked
                K4star1323;... %%% sym checked
                % 9th line
                K4star1123,... %%% sym checked
                K4star2223,... %%% sym checked
                K4star3323,... %%% sym checked
                -K4star1223,... %%% sym checked
                -K4star1323,... %%% sym checked
                K4star2323,... %%% sym checked
                -K4star1223,... %%% sym checked
                -K4star1323,... %%% sym checked
                K4star2323]; %%% sym checked
        end
        
    end
    
    methods (Static, Hidden)
        function out = skipCalc(newValue)
            persistent flag
            if isempty(flag), flag = false; end
            if nargin > 0, flag = logical(newValue); end
            out = flag;
        end
        
        function skip = consumeSkip()
            % Reads the flag AND sets it to false (single-use).
            skip = classes.SphereAssembly.skipCalc();
            classes.SphereAssembly.skipCalc(false);
        end
        
        function guard = suspendCalc()
            % Only for subclasses that modify the geometry after the parent constructor
            % Call this immediately before obj@Parent(...). 
            % Reactivates the calculation when 'guard' is destroyed.
            classes.SphereAssembly.skipCalc(true);
            guard = onCleanup(@() classes.SphereAssembly.skipCalc(false));
        end
    end

    methods (Access = protected)
        
        function DurlofskyTensors(obj)
            for i = 1:obj.N
                for j = [1:i - 1, i + 1:obj.N]
                    xij = obj.x{i} - obj.x{j};
                    obj.J1{i,j} = classes.SphereAssembly.tens_J1(xij, obj.a{j});    
                    obj.J2{i,j} = classes.SphereAssembly.tens_J2(xij, obj.a{i}, obj.a{j});
                    obj.J3{i,j} = classes.SphereAssembly.tens_J3(xij);
                    obj.J4_93{i,j} = classes.SphereAssembly.tens_J4_93(xij, obj.a{i}, obj.a{j});
                    obj.R0{i,j} = classes.SphereAssembly.tens_R0(xij);
                    obj.R3{i,j} = classes.SphereAssembly.tens_R3(xij);
                    obj.R4_93{i,j} = classes.SphereAssembly.tens_R4_93(xij);
                    obj.K1_39{i,j} = -classes.SphereAssembly.tens_K1_39(xij, obj.a{j});
                    obj.K2_39{i,j} = -classes.SphereAssembly.tens_K2_39(xij, obj.a{i}, obj.a{j});
                    obj.K3_39{i,j} = -classes.SphereAssembly.tens_K3_39(xij);
                    obj.K4_99{i,j} = -classes.SphereAssembly.tens_K4_99(xij, obj.a{i}, obj.a{j});
                end 
            end        
        end 
        % setup = 0 :
        function Calculation_gdM(obj)
            % initialization
            obj.gdM = zeros(obj.N*15);
            LJ = zeros(3, obj.N*3); % bloc J
            LR = zeros(3, obj.N*3); % bloc R
            LK = zeros(3, obj.N*9); % bloc K
            
            %%% 1st Faxén Law %%%
            for i = 1:obj.N
                for j = 1:obj.N
                    if i == j
                        LJ(:, 3*i - 2 : 3*i) = eye(3)/(6*pi*obj.a{i});
                        LR(:, 3*i - 2 : 3*i) = zeros(3);
                        LK(:, 9*i - 8 : 9*i) = zeros(3, 9);
                    else
                        LJ(:, 3*j - 2 : 3*j) = obj.J2{i,j}/(8*pi);
                        LR(:, 3*j - 2 : 3*j) = obj.R0{i,j}/(8*pi);
                        LK(:, 9*j - 8 : 9*j) = obj.K2_39{i,j}/(8*pi);
                    end
                end
                obj.gdM(3*i - 2 : 3*i, :) = [LJ, LR, LK];
            end
            %%% 2nd Faxén Law %%%
            for i = 1:obj.N
                for j = 1:obj.N
                    if i == j
                        LJ(:, 3*i - 2 : 3*i) = zeros(3);
                        LR(:, 3*i - 2 : 3*i) = eye(3)/(8*pi*obj.a{i}^3);
                        LK(:, 9*i - 8 : 9*i) = zeros(3, 9);
                    else
                        LJ(:, 3*j - 2 : 3*j) = obj.J3{i,j}/(8*pi);
                        LR(:, 3*j - 2 : 3*j) = obj.R3{i,j}/(8*pi);
                        LK(:, 9*j - 8 : 9*j) = obj.K3_39{i,j}/(8*pi);
                    end
                end
                obj.gdM(3*obj.N + 3*i - 2 : 3*obj.N + 3*i, :) = [LJ, LR, LK]; 
            end
            %%% 3rd Faxén Law %%%
            LJ = zeros(9, obj.N*3); % bloc J
            LR = zeros(9, obj.N*3); % bloc R
            LK = zeros(9, obj.N*9); % bloc K
            
            for i = 1:obj.N
                for j = 1:obj.N
                    if i == j
                        LJ(:, 3*i - 2 : 3*i) = zeros(9, 3);
                        LR(:, 3*i - 2 : 3*i) = zeros(9, 3);
                        LK(:, 9*i - 8 : 9*i) = eye(9)/(20*pi/3*obj.a{i}^3);
                    else
                        LJ(:, 3*j - 2 : 3*j) = obj.J4_93{i,j}/(8*pi);
                        LR(:, 3*j - 2 : 3*j) = obj.R4_93{i,j}/(8*pi);
                        LK(:, 9*j - 8 : 9*j) = obj.K4_99{i,j}/(8*pi);
                    end
                end
                obj.gdM(6*obj.N + 9*i - 8 : 6*obj.N + 9*i, :) = [LJ, LR, LK];
            end

            clearDurlofskyTensors(obj)
            
        end 
        
        % setup = 1 :
        function WallTensors(obj)
            for i = 1:obj.N
                for j = 1:obj.N
                    obj.J1star{i,j} = classes.SphereAssembly.tens_J1star(obj.x{i}, obj.x{j}, obj.a{j});    
                    obj.J2star{i,j} = classes.SphereAssembly.tens_J2star(obj.x{i}, obj.x{j}, obj.a{i}, obj.a{j});
                    obj.J3star{i,j} = classes.SphereAssembly.tens_J3star(obj.x{i}, obj.x{j}, obj.a{j});
                    obj.J4star93{i,j} = classes.SphereAssembly.tens_J4star93(obj.x{i}, obj.x{j}, obj.a{i}, obj.a{j});
                    obj.R1star{i,j} = classes.SphereAssembly.tens_R1star(obj.x{i}, obj.x{j});
                    obj.R2star{i,j} = classes.SphereAssembly.tens_R2star(obj.x{i}, obj.x{j}, obj.a{i});
                    obj.R3star{i,j} = classes.SphereAssembly.tens_R3star(obj.x{i}, obj.x{j});
                    obj.R4star93{i,j} = classes.SphereAssembly.tens_R4star93(obj.x{i}, obj.x{j}, obj.a{i});
                    obj.K1star39{i,j} = -classes.SphereAssembly.tens_K1star39(obj.x{i}, obj.x{j}, obj.a{j});
                    obj.K2star39{i,j} = -classes.SphereAssembly.tens_K2star39(obj.x{i}, obj.x{j}, obj.a{i}, obj.a{j});
                    obj.K3star39{i,j} = -classes.SphereAssembly.tens_K3star39(obj.x{i}, obj.x{j}, obj.a{j});
                    obj.K4star99{i,j} = -classes.SphereAssembly.tens_K4star99(obj.x{i}, obj.x{j}, obj.a{i}, obj.a{j});
                end
            end
        end
        function Calculation_gdM_Wall(obj)
            % initialization
            obj.gdM = zeros(obj.N*15);
            LJ = zeros(3, obj.N*3); % bloc J
            LR = zeros(3, obj.N*3); % bloc R
            LK = zeros(3, obj.N*9); % bloc K
            
            %%% 1st Faxén Law %%%
            for i = 1:obj.N
                for j = 1:obj.N
                    if i == j
                        LJ(:, 3*i - 2 : 3*i) = eye(3)/(6*pi*obj.a{i}) + obj.J2star{i,j}/(8*pi);
                        LR(:, 3*i - 2 : 3*i) = obj.R2star{i,j}/(8*pi);
                        LK(:, 9*i - 8 : 9*i) = obj.K2star39{i,j}/(8*pi);
                    else
                        LJ(:, 3*j - 2 : 3*j) = (obj.J2{i,j} + obj.J2star{i,j})/(8*pi);
                        LR(:, 3*j - 2 : 3*j) = (obj.R0{i,j} + obj.R2star{i,j})/(8*pi);
                        LK(:, 9*j - 8 : 9*j) = (obj.K2_39{i,j} + obj.K2star39{i,j})/(8*pi);
                    end
                end
                obj.gdM(3*i - 2 : 3*i, :) = [LJ, LR, LK];
            end
            %%% 2nd Faxén Law %%%
            for i = 1:obj.N
                for j = 1:obj.N
                    if i == j
                        LJ(:, 3*i - 2 : 3*i) = obj.J3star{i,j}/(8*pi);
                        LR(:, 3*i - 2 : 3*i) = eye(3)/(8*pi*obj.a{i}^3) + obj.R3star{i,j}/(8*pi);
                        LK(:, 9*i - 8 : 9*i) = obj.K3star39{i,j}/(8*pi);
                    else
                        LJ(:, 3*j - 2 : 3*j) = (obj.J3{i,j} + obj.J3star{i,j})/(8*pi);
                        LR(:, 3*j - 2 : 3*j) = (obj.R3{i,j} + obj.R3star{i,j})/(8*pi);
                        LK(:, 9*j - 8 : 9*j) = (obj.K3_39{i,j} + obj.K3star39{i,j})/(8*pi);
                    end
                end
                obj.gdM(3*obj.N + 3*i - 2 : 3*obj.N + 3*i, :) = [LJ, LR, LK]; 
            end
            %%% 3rd Faxén Law %%%
            LJ = zeros(9, obj.N*3); % bloc J
            LR = zeros(9, obj.N*3); % bloc R
            LK = zeros(9, obj.N*9); % bloc K

            for i = 1:obj.N
                for j = 1:obj.N
                    if i == j
                        LJ(:, 3*i - 2 : 3*i) = obj.J4star93{i,j}/(8*pi);
                        LR(:, 3*i - 2 : 3*i) = obj.R4star93{i,j}/(8*pi);
                        LK(:, 9*i - 8 : 9*i) = eye(9, 9)/(20*pi/3*obj.a{i}^3) + obj.K4star99{i,j}/(8*pi);
                    else
                        LJ(:, 3*j - 2 : 3*j) = (obj.J4_93{i,j} + obj.J4star93{i,j})/(8*pi);
                        LR(:, 3*j - 2 : 3*j) = (obj.R4_93{i,j} + obj.R4star93{i,j})/(8*pi);
                        LK(:, 9*j - 8 : 9*j) = (obj.K4_99{i,j} + obj.K4star99{i,j})/(8*pi);
                    end
                end
                obj.gdM(6*obj.N + 9*i - 8 : 6*obj.N + 9*i, :) = [LJ, LR, LK];
            end

            clearDurlofskyTensors(obj)
            
        end
        
    end

    methods
        
        %%% Constructor %%%
        function obj = SphereAssembly(xx, a, rho, setup)
            % xx : positions of spheres, 3*N matrix (ith column for the ith sphere)
            % a : radius of the spheres, either a scalar (monodisperse) or N components vector (polydisperse)
            % rho : mass density of the spheres, either a scalar (homogenous) or N components vector (heterogeneous)
            % setup : 0 or 1
            
            skip = classes.SphereAssembly.consumeSkip();
            
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
            
            [~, obj.N] = size(xx);
            
            if isscalar(rho)
                rho = rho.*ones(obj.N,1);
            end
            if isscalar(a)
                a = a.*ones(obj.N,1);
            end
            if isrow(rho)
                rho = rho';
            end
            if isrow(a)
                a = a';
            end

            % errors
            if setup == 1 
                % checking all the spheres are on the same side of the wall
                wall_conditions = [nnz(xx(1,:)>0) == obj.N, nnz(xx(1,:)<0) == obj.N];
                if not(or(wall_conditions(1), wall_conditions(2)))
                    % either all x components of the spheres position are positive
                    % or all x components of the spheres position are negative
                    % if not -> error
                    error('All the spheres must be on the same side of the wall.')
                end
                if any(abs(xx(1,:))<a')
                    error('Some spheres are overlapping the wall.')
                end
            elseif and(setup ~= 0 ,setup ~= 1) % checking the setup exists
                error('Choose an existing setup : \n%s, %s', '0', '1')
            end
            
            m = rho.*(4/3).*pi.*a.^3;
            
            obj.x = num2cell(xx, 1);
            obj.m = num2cell(m);
            obj.a = num2cell(a);
            obj.rho = num2cell(rho);
            
            obj.setup = setup;
            
            %%%% only for subclasses that modify the geometry after the parent constructor
            if skip
                return
            end
            %%%%
            
            DurlofskyTensors(obj);
            if setup == 0
                Calculation_gdM(obj);
            elseif setup == 1
                WallTensors(obj)
                Calculation_gdM_Wall(obj)
            end
            
        end
        
        function clearDurlofskyTensors(obj)
            obj.J1 = [];
            obj.J2 = [];
            obj.J3 = [];
            obj.J4_93 = [];
            obj.R0 = [];
            obj.R3 = [];
            obj.R4_93 = [];
            obj.K1_39 = [];
            obj.K2_39 = [];
            obj.K3_39 = [];
            obj.K4_99 = [];
            
            if obj.setup == 1 
                obj.J1star = [];
                obj.J2star = [];
                obj.J3star = [];
                obj.J4star93 = [];
                obj.R1star = [];
                obj.R2star = [];
                obj.R3star = [];
                obj.R4star93 = [];
                obj.K1star39 = [];
                obj.K2star39 = [];
                obj.K3star39 = [];
                obj.K4star99 = [];
            end
        end
        
        function [Uf, varargout] = FlowField(xc, obj, flow)
            % you need to compute FTS first
            
            arguments
                xc (3,1) double % point of interest
                obj classes.SphereAssembly
                flow.mu (1,1) double = 1
                flow.U0 (3,1) double = zeros(3,1)
                flow.Omega_inf (3,1) double = zeros(3,1)
                flow.E_inf (9,1) double = zeros(9,1)
            end
            
            if isempty(obj.Force)
                message = ['You need to compute forces, torques and stresslets exerted by the spheres.' newline 'Try getFTS(obj)'];
                error(message)
            end

            % ambient flow (lab frame)
            u_inf = @(x) flow.U0 + cross(flow.Omega_inf, x)...
                + aux.recTransfo(flow.E_inf)*x; 

            Uf = zeros(3,1);
            
            if obj.setup == 0
                for i = 1:obj.N
                    Uf = Uf + 1/(8*pi*flow.mu)*(classes.SphereAssembly.tens_J1(xc - obj.x{i}, obj.a{i})*obj.Force{i} ...
                        + classes.SphereAssembly.tens_R0(xc - obj.x{i})*obj.Torque{i} ...
                        - classes.SphereAssembly.tens_K1_39(xc - obj.x{i}, obj.a{i})*obj.Stresslet{i});
                end
            elseif obj.setup == 1
                for i = 1:obj.N
                    Uf = Uf + 1/(8*pi*flow.mu)*((classes.SphereAssembly.tens_J1(xc - obj.x{i}, obj.a{i}) + classes.SphereAssembly.tens_J1star(xc, obj.x{i}, obj.a{i}))*obj.Force{i}...
                        + (classes.SphereAssembly.tens_R0(xc - obj.x{i}) + classes.SphereAssembly.tens_R1star(xc, obj.x{i}))*obj.Torque{i}...
                        + (classes.SphereAssembly.tens_K1_39(xc - obj.x{i}, obj.a{i}) + classes.SphereAssembly.tens_K1star39(xc, obj.x{i}, obj.a{i}))*obj.Stresslet{i});
                end
            end
            
            Uf = Uf + u_inf(xc);

            if nargout == 2
                U = Uf;
                Uf = U(1);
                varargout{1} = U(2);
            elseif nargout == 3
                U = Uf;
                Uf = U(1);
                varargout{1} = U(2);
                varargout{2} = U(3);
            end
        end
        
        function Show(obj, color, settings)
            
            arguments 
                obj (1,1) classes.SphereAssembly 
                color (:, 3) double = [0.75, 0.7, 0.65]
                settings.preset string {mustBeMember(settings.preset, ["default", "sketch"])} = "default"
            end

            if or(nargin == 1, all(size(color) == [1,3]))
                color = color.*ones(obj.N, 1); 
            end

            % figure
            hold on
            
            [X0,Y0,Z0] = sphere(40);
            for i = 1:obj.N
                % allocating the correct radius to each sphere
                X1 = obj.a{i}*X0;
                Y1 = obj.a{i}*Y0;
                Z1 = obj.a{i}*Z0;
                % plotting
                s = surf(X1 + obj.x{i}(1), Y1 + obj.x{i}(2), Z1 + obj.x{i}(3)); ...
                set(s, 'FaceColor', color(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.7)
                view(3)
            end
            
            if obj.setup == 1 % wall 
                lim = axis;
                if lim(1) > 0
                    lim(1) = 0;
                elseif lim(2) < 0
                    lim(2) = 0;
                end
                [wall_y, wall_z] = meshgrid(linspace(lim(3), lim(4), 30),  linspace(lim(5), lim(6), 30)); % Generate x and y data
                wall_x = zeros(size(wall_y)); % Generate z data
                surf(wall_x, wall_y, wall_z, ...
                'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            end
            
            axis equal

            if settings.preset == "default"
                camlight headlight
                lighting gouraud
                material dull
                % material shiny
                view(3)
                grid on
                set(gca, 'FontSize', 12, 'LineWidth', 1.2)
                xlabel('$x$', 'FontSize', 15, 'Interpreter', 'latex')
                ylabel('$y$', 'FontSize', 15, 'Interpreter', 'latex')
                zlabel('$z$', 'FontSize', 15, 'Interpreter', 'latex')
            elseif settings.preset == "sketch"
                view(30, -25)
                camlight headlight
                lighting gouraud
                material shiny
                % view(30, 25)
                view(95,10)
                grid off
                box off
                axis off
            end
            
            set(gcf, 'Color', 'w')
            hold off
        end
        
    end
    
end
