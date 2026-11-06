clear 
close all

%% Take a look on the swimmer %%
obj = classes.Puller;
figure(1)
Show(obj)
view(2)

%% Motion of the swimmer %%
vid_filename = "puller_swimming_setup0";
period = 0.5;
getSwimmingLaws(obj, T = period) 
timestep = 0.01;
tstart = 0;
tend = 10;
[t_out,y_out] = Dynamics(obj, tspan = tstart:timestep:tend, xB = zeros(3,1));
%%
figure(2)
Animation(obj, t_out, y_out, orient = [0, 90], option = "yes", filename = vid_filename, framerate = 1/timestep, layout="fixed")

%% Flow field and streamlines

% position of the body
xB = zeros(3,1);
qB = [1 0 0 0];
% vectors grid
[X, Y] = meshgrid(-0.02 : 0.02/22 : 0.02, -0.02 : 0.02/22 : 0.02);

% choose the moment you are interested in
t = 0.13;
f = zeros(6,1);
[U, V, ~, ~, ~] = getFlowField(obj, t, f, xB, qB, X, Y);

figure(3)
q = quiver(X, Y, U, V, 1.5);
q.LineWidth = 1.2;
axis equal

hold on
xmesh = [-0.02 -0.016 -0.012 -0.008 -0.0041 -0.000001 0.000001 0.0041 0.008 0.012 0.0200];
ymesh = [-0.02 -0.005 0 0.002 0.005 0.01  0.02];
[startX, startY] = meshgrid(xmesh, ymesh);
verts = stream2(X, Y, U, V, startX, startY);
lineobj = streamline(verts);
for i = 1:length(verts)
    lineobj(i).LineWidth = 2;
end
axis([-0.02 0.02 -0.02 0.02])
ax = gca;
ax.FontSize = 15; 
xticks([-0.015 -0.01 -0.005 0 0.005 0.01 0.015 0.02])
xticklabels([-15 -10 -5 0 5 10 15 20])
yticklabels([-20 -15 -10 -5 0 5 10 15 20])
xlabel('$x$ [$\mu m$]', 'FontSize', 25, 'Interpreter', 'latex', 'FontWeight', 'bold')
ylabel('$y$ [$\mu m$]', 'FontSize', 25, 'Interpreter', 'latex', 'FontWeight', 'bold')

pos = [obj.x{1}(1)-obj.a{1}, obj.x{1}(2)-obj.a{1}, 2*obj.a{1}, 2*obj.a{1}];
rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', [0. 0. 0.], 'FaceAlpha', 0.7)
for i = 2:obj.N
    pos = [obj.x{i}(1)-obj.a{i}, obj.x{i}(2)-obj.a{i}, 2*obj.a{i}, 2*obj.a{i}];
    rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'k')
end

%% Near a plane wall

objwall = classes.Puller(xB = [0.03; 0; 0], setup = 1);
figure(4)
Show(objwall)
%%
vid_filename_wall = "puller_swimming_setup1";
period = 0.5;
getSwimmingLaws(objwall, T = period) 
timestep = 0.01;
tstart = 0;
tend = 35;
Rot = aux.rot_u([0,0,1], pi/18);
qb = aux.matrice_to_quaternion(Rot);
[t_out_wall,y_out_wall] = Dynamics(objwall, tspan = tstart:timestep:tend, xB = [0.02; 0; 0], qB = qb);

%%
figure(5)
Animation(objwall, t_out_wall, y_out_wall, orient = [0, 90], option = "yes", filename = vid_filename_wall, framerate = 1/timestep, layout = "fixed")

%%
% choose the moment you are interested in
t = 3;
index_t = round((t - tstart)/timestep) + 1;
% position of the body
xB = y_out_wall(index_t, 5:7)';
qB = y_out_wall(index_t, 1:4);
% vectors grid
[X, Y] = meshgrid(-10*objwall.a{1} + xB(1) : objwall.a{1}/2 : 10*objwall.a{1} + xB(1), -10*objwall.a{1} + xB(2) : objwall.a{1}/2 : 10*objwall.a{1} + xB(2));

[U, V, ~, ~, ~] = getFlowField(objwall, t, zeros(6,1), xB, qB, X, Y);

figure(6)
quiver(X, Y, U, V, 2)
axis equal
hold on
Show(objwall)
view(2)
[startX, startY] = meshgrid(-10*objwall.a{1} + xB(1) : 2*objwall.a{1} : 10*objwall.a{1} + xB(1), -10*objwall.a{1} + xB(2) : 2*objwall.a{1} : 10*objwall.a{1} + xB(2));
verts_wall = stream2(X, Y, U, V, startX, startY);
lineobj_wall = streamline(verts_wall);
xl = xline(0, 'k');