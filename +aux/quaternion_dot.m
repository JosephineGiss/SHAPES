function dq = quaternion_dot(Q,Omega, varargin)

% 
if nargin == 2
    varargin{1} = "mobile"; % frame
end

if varargin{1} == "lab" % frame
    omega = Omega;
else
    R = quaternion_to_matrice(Q);
    omega = R*Omega;
end

norme_Q=sqrt(Q(1)^2+Q(2)^2+Q(3)^2+Q(4)^2);
Q(1)=Q(1)/norme_Q;
Q(2)=Q(2)/norme_Q;
Q(3)=Q(3)/norme_Q;
Q(4)=Q(4)/norme_Q;
    

%integration quaternion

A_omega = [       0, -omega(1), -omega(2), -omega(3)
           omega(1),         0, -omega(3),  omega(2)
           omega(2),  omega(3),         0, -omega(1)
           omega(3), -omega(2),  omega(1),        0];


dq  = 1/2*A_omega*Q;

