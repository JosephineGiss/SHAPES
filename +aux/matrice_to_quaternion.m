function [Q] = matrice_to_quaternion(R)

%QUATERNION_TO_MATRICE: Passage quaternion -> matrice:

Q=zeros(4,1);

Q(1)=sqrt(1+trace(R))/2;
Q(2)=(R(3,2)-R(2,3))/(4*Q(1)); 
Q(3)=(R(1,3)-R(3,1))/(4*Q(1)); 
Q(4)=(R(2,1)-R(1,2))/(4*Q(1)); 

