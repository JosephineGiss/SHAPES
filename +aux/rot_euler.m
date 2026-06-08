function R = rot_euler(phi, theta, psi)
% returns the rotation matrix using Euler angles  

    cphi = cos(phi);
    sphi = sin(phi);
    cth = cos(theta);
    sth = sin(theta);
    cpsi = cos(psi);
    spsi = sin(psi);
    
    R = [cphi*cpsi - sphi*cth*spsi, -cphi*spsi - sphi*cth*cpsi, sphi*sth;...
        sphi*cpsi + cphi*cth*spsi, -sphi*spsi + cphi*cth*cpsi, -cphi*sth;...
        sth*spsi, sth*cpsi, cth];
end