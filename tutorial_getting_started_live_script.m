clear
close all
%[text] # Getting started
%[text] This tutorial introduces the simplest features of the SHAPES program package and focuses on rigid particles. Here, we show how to instantiate a rigid aggregate using the RigidBody class, compute its resistance tensors and simulate its dynamics, either in an unbounded fluid or near a wall.
%%
%[text] ## 1) Rigid-aggregate shapes
%[text] For each aggregate, the user must choose the coordinates xx and radii a of the beads. The i-th column of the matrix xx contains the coordinates of the bead i. The user must also ensure that beads do not touch or overlap. If all the beads have the same radius, the variable a can be a scalar.
  % runs the section %[control:button:5963]{"position":[1,2]}
xx1 = [0; 0; 0];
xx2 = [5; 0; 0];
xx3 = [5; 0; 5];
xx4 = [5; 5; 5];
% Coordinates for an object made of 4 sphere
xx = [xx1, xx2, xx3, xx4];

a = [1, 1.1, 1.4, 1]; % Radii of the spheres
% a = 1; % same radius for all the spheres

obj = classes.RigidBody(xx,a);

figure
Show(obj)
%%
%[text] ## 2) Resistance tensors and mobility tensors
%[text] In SHAPES, one can compute resistance and mobility tensor of a given aggregate. By default, if no centre is specified, resistance tensors are computed from the centre of mass, but it can be computed from any point xr. Here, as an example, we use the centre of resistance as our reference point and display the tensor B in the command window:
  %[control:button:3539]{"position":[1,2]}
xr = CenterOfResistance(obj);
getResistanceTensors(obj, xr)
disp("Tensor B computed from xr =")
disp(obj.B)
%[text] When you compute resistance tensors multiple times, only the last result is kept. You can check from which point you got your resistance tensors with xres property.
disp("The last time you made the calculations, resistance tensors were computed from:")
disp(obj.xres)
%[text] You can also obtain mobility tensors and compute the centre of mobility.
getMobilityTensors(obj)
xm = CenterOfMobility(obj)
%[text] Mobility tensors are always computed from the same points used for resistance tensors and you can double check which point it is by using xmob property. or this aggregate, the centre of mobility is slightly offset from the centre of resistance.
%%
%[text] ## 3) Interactions with a wall
%[text] Hydrodynamic interactions between a particle and a planar wall can also be studied. Creating a RigidBody object near a wall requires specifying four input arguments: the mass density as the third and the setup as the fourth. Mass density is only useful for dynamic simulations. The setup is either 0 (no wall, by default) or 1 (wall present). To instantiate an object near a wall, all the beads must be on the same side of the wall (which is located at x = 0), but the side itself does not matter.
  %[control:button:8e92]{"position":[1,2]}
N = 4; % number of spheres
a = 0.1; % radius of the spheres
% coordinates of the spheres
xx = zeros(3, N);
for i = 1:N
    xx(:, i) = [0.5*i; 0; 0]; 
end
rho = 1030; % mass density of the particle
setup = 1;  % setup

obj_near_a_wall = classes.RigidBody(xx, a, rho, setup);

figure
Show(obj_near_a_wall)
%[text] Resistance tensors now account for the interactions with the wall.
getResistanceTensors(obj_near_a_wall);
disp("Tensor B computed for a rod near a wall =");
disp(obj_near_a_wall.B);
%%
%[text] ## 4) Pre-existing classes
%[text] Several shapes are available in the SHAPES program package, a comprehensive list can be found in the +classes file of the package. Each of them has their own properties. Here, we made a two-bladed propeller (see Happel and Brenner \[27, Figure 5-4.1\] or Kim and Karrila \[28, Figure 5.3\]) using the ScrewPropeller subclass. The main input arguments are:
%[text] - rad, the radius of each disk;
%[text] - h, the distance between the center of the two disks;
%[text] - theta, the relative orientation of the two disks. \
  %[control:button:53c3]{"position":[1,2]}
rad = 12; 
h = 30;   
theta = pi/4; 
prop = classes.ScrewPropeller(rad, h, theta);

figure
Show(prop)
view(95, 15)

getResistanceTensors(prop); % about the centre of mass
disp("Tensor B computed for a screw propeller =");
disp(prop.B);
%[text] Most classes have optional arguments useful to modify the position and orientation of the particle. Each particle has a default position, specified by the center of mass, and a default orientation, specified by a rotation matrix. The package includes auxiliary functions that are not class methods but are useful for performing certain mathematical operations. For instance, two functions are available to compute a rotation matrix:
%[text] - R = aux.rot\_u(u, alpha) returns the rotation matrix corresponding to the rotation of an angle 'alpha' around the 'u' axis,
%[text] - R = aux.rot\_euler(phi, theta, psi) returns the rotation matrix using Euler angles (parametrization by Goldstein). \
%[text] Here, we use the first option:
R = aux.rot_u([1 0 0], pi/2);
prop_horizontal = classes.ScrewPropeller(rad, h, theta, xG = [5;0;0], rot = R);

figure
Show(prop_horizontal)
view(95, 15)

getResistanceTensors(prop_horizontal);
disp("Tensor B computed for a screw propeller (new orientation) =");
disp(prop_horizontal.B);
%%
%[text] ## 5) Aggregate dynamics
%[text] The translational and angular dynamics of an aggregate is calculated in the same way for any class or subclass. The input arguments are the aggregate obj in question, as well as any external forces and torques FText. In addition, one can also choose different linear velocity fields, the dynamic viscosity mu, as well as location and orientation of the aggregate at the beginning of the simulation. Here, the particles are assumed to be non-inertial so fluid and particle densities should ideally be very close. While the internal integration time steps are automatically adapted by the solver to guarantee precision, the user-defined output time steps can be chosen arbitrarily to obtain data or generate animations without affecting the underlying resolution of the dynamics. As an example, consider a sedimenting propeller in a fluid at rest.
  %[control:button:8df1]{"position":[1,2]}
rho_f = 0997;   % density of the fluid
rho_p = 1000;   % density of the particle
FText = [0; 0; prop.N*(4/3)*pi*prop.a{1}^3*(rho_p - rho_f)*(-9.81); 0; 0; 0]; % external forces, torques
q_init = matrice_to_quaternion(aux.rot_u([0 0 1], theta)); % orientation of the particle at the start of the simulation

t_step = 0.05;
t_start = 0;
t_end = 40;
[t_out,y_out] = Dynamics(prop, FText, tspan = t_start:t_step:t_end, mu = 0.1, q0 = q_init);
%[text] We solve for particle dynamics using the Runge–Kutta–Fehlberg numerical scheme. The vector t\_out stores the time steps (which are chosen using the tspan option in the inputs). The first four columns of y\_out store quaternions (parameterisation of the rotation matrix) and the position of the centre of mass is stored in columns 5 through 7. Animation of the motion can be made using space data over time and can be saved as a video:
%%
  %[control:button:6e8f]{"position":[1,2]}
% this section can last a little longer (5-10min)
% each frame of the video is printed at the screen 
figure
Animation(prop, t_out, y_out, video = "yes", filename = "sedimentation_propeller",...
    orient = [90 0], framerate=1/t_step)
%[text] Sometimes, it may be helpful to adjust the frame rate (or the viscosity mu) to speed up or slow down the animations.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":25.1}
%---
%[control:button:5963]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:3539]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:8e92]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:53c3]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:8df1]
%   data: {"label":"Run","run":"Section"}
%---
%[control:button:6e8f]
%   data: {"label":"Run","run":"Section"}
%---
