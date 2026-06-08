function R = rot_u(u, alpha) 
% rotation of an angle 'alpha' around the 'u' axis
% returns the corresponding rotation matrix
% u must be an unit vector

    if not(norm(u) == 1)
        error('u must be an unit vector.')
    end
    
    c = cos(alpha);
    s = sin(alpha);

    R = [u(1)^2*(1 - c) + c, u(1)*u(2)*(1 - c) - u(3)*s, u(1)*u(3)*(1 - c) + u(2)*s;...
        u(1)*u(2)*(1 - c) + u(3)*s, u(2)^2*(1 - c) + c, u(2)*u(3)*(1 - c) - u(1)*s;...
        u(1)*u(3)*(1 - c) - u(2)*s, u(2)*u(3)*(1 - c) + u(1)*s, u(3)^2*(1 - c) + c];
end