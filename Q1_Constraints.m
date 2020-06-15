%% CORRECT

clc 
clear all

syms a(t) b(t) g(t) Th(t) x(t) y(t) % alpha, beta, gamma, Theta, x, y
syms da db dg  R

%T = 0;  % Angle Theta

aa = diff(a(t),t);
gg = diff(g(t),t);
bb = diff(b(t),t);
xx = diff(x(t),t);
yy = diff(y(t),t);

aaa = diff(a(t),t,t);
ggg = diff(g(t),t,t);
bbb = diff(b(t),t,t);
xxx = diff(x(t),t,t);
yyy = diff(y(t),t,t);

% Rotation Matrices
R_12 = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
R_23 = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
R_34 = [1 0 0 ; 0 cos(g) -sin(g); 0 sin(g) cos(g)];

R_21 = transpose(R_12);
R_32 = transpose(R_23);
R_43 = transpose(R_34);

R_14 = R_12*R_23*R_34;
R_41 = transpose(R_14);

% Angular Velocities
w1_21 = [0;0;aa];
w2_32 = [0;bb;0];
w3_43 = [gg;0;0];

w2_21 = w1_21;
w3_32 = w2_32;
w4_43 = w3_43;

w4_41 = w4_43+ R_43*w3_32 + R_43*R_32*w2_21;

r1_AC = R_12*R_23*[0;0;-R]; % correct
%r1_AC = [0;0;-R*cos(b)];
r1_OA = [x;y;R*cos(b)]; % correct
%r1_OA = [x;y;R)];
syms rx ry rz

r4_AP = [rx;ry;rz];
r4_AC = R_41*r1_AC;   
r4_AC = r4_AC(t); % make sym not symfun

rr1_OA = diff(r1_OA,t);   
rr4_AP = cross(w4_41,r4_AP)
rr4_AC = simplify(subs(rr4_AP,{'rx','ry','rz'},{r4_AC(1),r4_AC(2),r4_AC(3)}))
rr1_AC = R_14*rr4_AC;
rr1_OC = simplify(rr1_OA+rr1_AC)
rr1_OC = subs(rr1_OC,{a(t),b(t),g(t),x(t),y(t),Th(t)...
            diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t)...
            diff(a(t), t, t),diff(b(t), t, t), diff(g(t), t, t),diff(x(t), t, t), diff(y(t), t, t),diff(Th(t), t, t)},...
            {'a','b','g','x','y','Th','aa','bb','gg','xx','yy',0,'aaa','bbb','ggg','xxx','yyy',0}) 