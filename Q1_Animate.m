clc
clear all

% Potential Issues
% animation Issues


%% Constants

 syms t;
 VIDEO = 1;

grav = 9.81;
diam = 25;
R = diam;
m = 1;


% angles
Th = pi/10; % angle of plane
a_init = -pi/4;
b_init = 0;
g_init = 0;

aa_init = -4;
bb_init = 2;
gg_init = 5;
%% 

%Plane Dimensions
leng = 800;  % length of plane
width = 2000;  % width of plane

%Initial position of center of disk
x_init = 300;
y_init = 2000;


%% Intitial Conditions comply with constraint equations
xx_init = -( - R*gg_init*sin(a_init) - R*bb_init*cos(a_init)*cos(b_init) + R*aa_init*sin(a_init)*sin(b_init));
yy_init = -( + R*gg_init*cos(a_init) - R*aa_init*cos(a_init)*sin(b_init) - R*bb_init*cos(b_init)*sin(a_init));
%% Solve EOM

x_init = [a_init;b_init;g_init;x_init;y_init;aa_init;bb_init;gg_init;xx_init;yy_init];

% IC's for clockwise (negative rotation)
% x_init = [pi/20; pi/11; pi/20;0; -pi/10;0;0;-1000];
tspan = [0 15];                                 % start and finish times
options = odeset('RelTol',1e-7,'AbsTol',1e-7);  % solver options
sol = ode45(@eom, tspan, x_init, options);      % solve the eoms
dt = 0.01;                                      % set time step
t = tspan(1) : dt : tspan(2);                   % create time vector
X = deval(sol,t);                               % evaluate solution

 %% Plotting States
    h = plot(t,X)
    xlabel('time')
    ylabel('states')
    leg = legend('$\alpha$','$\beta$','$\gamma$','$x$','$y$','$\dot{\alpha}$','$\dot{\beta}$','$\dot{\gamma}$','$\dot{x}$','$\dot{y}$');
    set(leg,'Interpreter','latex')
    set(h(1), 'color', 'b');
    set(h(2), 'color', 'r');
    set(h(3), 'color', 'c');
    set(h(4), 'color', 'g');
    set(h(5), 'color', 'm');
    set(h(6), 'color', 'b');
    set(h(6), 'LineStyle', '--');
    set(h(7), 'color', 'r');
    set(h(7), 'LineStyle', '--');
    set(h(8), 'color', 'c');
    set(h(8), 'LineStyle', '--');
    set(h(9), 'color', 'g');
    set(h(9), 'LineStyle', '--');
    set(h(10), 'color', 'm');
    set(h(10), 'LineStyle', '--');
    figure

%% Rotation Matrices

syms a b g x y aa bb xx yy gg

R_01 = [cos(Th) 0 -sin(Th); 0 1 0; sin(Th) 0 cos(Th)];
R_12 = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
R_23 = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
R_34 = [1 0 0; 0 cos(g) -sin(g); 0 sin(g) cos(g)];

R_10 = transpose(R_01);
R_21 = transpose(R_12);
R_32 = transpose(R_23);
R_43 = transpose(R_34);

R_04 = R_01*R_12*R_23*R_34;
R_40 = transpose(R_04);


%% Plane
x = leng*cos(Th);
y = width;
z = leng*sin(Th);

[X_p,Y_p] = meshgrid(1:10:x,1:10:y);
Z_p = X_p*sin(Th);

% %% Disk
% % Convert x & z co-ordinates from generalised co-ordinates to drawing
% % Generlaised co-ordinates; normal vector is along x
% % Drawing; normal vector is along z
normal_itit_4 = [1;0;0];
normal = ((R_04)*normal_itit_4)';

% for 1r_0A
r1_OA = [x;y;diam*cos(b)]
center = (R_01*r1_OA)';


% Checking Vectors
rr1_OC_vec = ones(3,length(t));
other1_vec = ones(1,length(t));
other2_vec = ones(1,length(t));



rr1_OC = [xx - R*gg*sin(a) - R*bb*cos(a)*cos(b) + R*aa*sin(a)*sin(b);...
         yy + R*gg*cos(a) - R*aa*cos(a)*sin(b) - R*bb*cos(b)*sin(a);0];
     
     
rr1_OA = [xx;yy;-R*sin(b)*bb];      
     
    for i = 1:length(t)
    cla

    % EOM to angle
    a = X(1,i);
    b = X(2,i);
    g = X(3,i);
    x = X(4,i);
    y = X(5,i);
    aa = X(6,i);
    bb = X(7,i);
    gg = X(8,i);
    xx = X(9,i);
    yy = X(10,i);
    
    rr1_OC_vec(:,i) = eval(rr1_OC);
    
    

    
    end
     
 plot(10*t/length(t),rr1_OC_vec)
 xlabel('time');
 ylabel('velocity of contact point')
 title('Velocity of contact point')
 legend({'x','y','z'})
 figure
     

%% Animate

    if VIDEO
        fps = 1/dt;
        MyVideo = VideoWriter('RotatingDisk','MPEG-4');
        MyVideo.FrameRate = fps;
        open(MyVideo);
    end

%% Animation Loop
handle = figure;
surf(X_p,Y_p,Z_p);
view(280,45)
hold on
xlabel('x')
ylabel('y')
zlabel('z')



for i = 1:length(t)
    cla

    % EOM to angle
    a = X(1,i);
    b = X(2,i);
    g = X(3,i);
    x = X(4,i);
    y = X(5,i);
    
    R_01 = [cos(Th) 0 -sin(Th); 0 1 0; sin(Th) 0 cos(Th)];
    R_12 = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
    R_23 = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
    R_34 = [1 0 0; 0 cos(g) -sin(g); 0 sin(g) cos(g)];

    R_04 = R_01*R_12*R_23*R_34;
    R_40 = transpose(R_04);
    
    normal_itit_4 = [1;0;0];
    normal = ((R_04)*normal_itit_4)';

    % for 1r_0A
    r1_OA = [x;y;diam*cos(b)];
    r0_OA = R_01*r1_OA;
    center = (r0_OA)';
  
    
    surf(X_p,Y_p,Z_p); hold on

    % Note that in some cases red & yellow will flip randomly due to
    % animation error
    plotDisk3D(center,normal,diam,g)
    %% 
    txt = ['Time: ' num2str(i/length(t)*10) 'seconds'];
    text(x+6*diam,y-6*diam,z+10*diam,txt)
    
    r = 3; % radius of little tracking spere
    [xs,ys,zs] = sphere;
    r0_AC = R_01*R_12*R_23*[0;0;-diam];
    r0_OA = R_01*r1_OA;
    surf(r*xs+r0_OA(1),r*ys+r0_OA(2),r*zs+r0_OA(3))
    %surf(r*xs+r0_AC(1)+r0_OA(1),r*ys+r0_AC(2)+r0_OA(2),r*zs+r0_AC(3)+r0_OA(3))        
    
    lim = 10*diam
    xlim([x-lim x+lim])
    ylim([y-lim y+lim])
    zlim([z-lim z+lim])
  
    % Graphing axis stuff
   
        if VIDEO
            writeVideo(MyVideo,getframe(handle));
        else
            pause(dt)
        end
 end

    if VIDEO
    close(MyVideo)
    end
    % notifies via terminal when recording is finished
fprintf("Animation recorded\n");


% equations of motion taken from our "ExtractEquations.m" file
function xdot = eom(t,X)
% for ode 45
    a = X(1);
    b = X(2);
    g = X(3);
    x = X(4);
    y = X(5);
    aa = X(6);
    bb = X(7);
    gg = X(8);
    xx = X(9);
    yy = X(10);
    
    %Constants = [grav,diam,R,m];
%     grav = Constants(1);
%     diam = Constants(2);
%     R = Constants(3);
%     m = Constants(4);
grav = 9.81;
diam = 25;
R = diam;
m = 1;
Th = pi/10;
% equations of motion

l1 = -(4*grav*m*cos(b)^3*sin(Th) - 5*grav*m*cos(b)*sin(Th) + 8*grav*m*cos(a)^2*cos(b)^3*sin(Th) - 4*grav*m*cos(b)^3*sin(Th)*sin(a)^2 - 10*grav*m*cos(a)^2*cos(b)*sin(Th) + 15*R*aa^2*m*cos(a)*cos(b)*sin(b) + 15*R*bb^2*m*cos(a)*cos(b)*sin(b) - 15*R*aa*gg*m*cos(a)*cos(b) - 15*R*aa^2*m*cos(a)*cos(b)^3*sin(b) + 12*R*bb^2*m*cos(a)*cos(b)^3*sin(b) - 12*grav*m*cos(Th)*cos(a)*cos(b)^2*sin(b) + 18*R*aa*gg*m*cos(a)*cos(b)^3 + 5*R*aa*bb*m*cos(b)^2*sin(a))/(5*cos(b) + 4*cos(b)^3 + 10*cos(a)^2*cos(b) + 10*cos(b)*sin(a)^2 - 4*cos(a)^2*cos(b)^3 - 4*cos(b)^3*sin(a)^2);
l2 =-(15*R*aa^2*m*cos(b)*sin(a)*sin(b) + 15*R*bb^2*m*cos(b)*sin(a)*sin(b) - 10*grav*m*cos(a)*cos(b)*sin(Th)*sin(a) - 15*R*aa*gg*m*cos(b)*sin(a) - 15*R*aa^2*m*cos(b)^3*sin(a)*sin(b) + 12*R*bb^2*m*cos(b)^3*sin(a)*sin(b) + 12*grav*m*cos(a)*cos(b)^3*sin(Th)*sin(a) - 12*grav*m*cos(Th)*cos(b)^2*sin(a)*sin(b) - 5*R*aa*bb*m*cos(a)*cos(b)^2 + 18*R*aa*gg*m*cos(b)^3*sin(a))/(5*cos(b) + 4*cos(b)^3 + 10*cos(a)^2*cos(b) + 10*cos(b)*sin(a)^2 - 4*cos(a)^2*cos(b)^3 - 4*cos(b)^3*sin(a)^2);
 
 
aaa = (2*bb*gg)/cos(b);
bbb = (4*R*m*sin(2*b)*bb^2 + 4*l1*cos(a)*cos(b) + 4*l2*cos(b)*sin(a) - (R*aa^2*m*sin(2*b))/2 - 4*grav*m*cos(Th)*sin(b) + 2*R*aa*gg*m*cos(b))/(R*m*(4*cos(b)^2 - 5));
ggg = ((cos(b)^2 - 2)*(2*l2*cos(a) - 2*l1*sin(a) + R*aa*bb*m*cos(b)))/(R*m*(sin(b)^2 - 1)) - (2*sin(b)*(2*l2*cos(a)*sin(b) - 2*l1*sin(a)*sin(b) - R*bb*gg*m*cos(b) + R*aa*bb*m*cos(b)*sin(b)))/(R*m*cos(b)^2);
xxx = (l1 - grav*m*sin(Th))/m;
yyy = l2/m;

% matrix of xdots to be used
    xdot = [aa;bb;gg;xx;yy;aaa;bbb;ggg;xxx;yyy];
end





