clc 
clear all

syms a(t) b(t) g(t) Th(t) x(t) y(t)  % alpha, beta, gamma, Theta, x, y
syms aa(t) bb(t) gg(t) xx(t) yy(t) 
syms da db dg dx dy R m

aa(t) = diff(a(t),t);
bb(t) = diff(b(t),t);
gg(t) = diff(g(t),t);
xx(t) = diff(x(t),t);
yy(t) = diff(y(t),t);

aaa(t) = diff(a(t),t,t);
bbb(t) = diff(b(t),t,t);
ggg(t) = diff(g(t),t,t);
xxx(t) = diff(x(t),t,t);
yyy(t) = diff(y(t),t,t);

%% Lagrange L = T - V
% create symbols for derivatives
a_ = sym('a');
b_ = sym('b');
g_ = sym('g');
x_ = sym('x');
y_ = sym('y');
Th_ = sym('Th');
aa_ = sym('aa');
bb_ = sym('bb');
gg_ = sym('gg');
xx_ = sym('xx');
yy_ = sym('yy');
TTh_ = sym('TTh');
aaa_ = sym('aaa');
bbb_ = sym('bbb');
ggg_ = sym('ggg');
xxx_ = sym('xxx');
yyy_ = sym('yyy');
TTTh_ = sym('TTTh');

assume(a(t),'real');
assumeAlso(b(t),'real');
assumeAlso(g(t),'real');
assumeAlso(x(t),'real');
assumeAlso(y(t),'real');
assumeAlso(Th(t),'real');
assumeAlso(R,'real');
assumeAlso(m,'real');


assumeAlso(diff(a(t),t),'real');
assumeAlso(diff(b(t),t),'real');
assumeAlso(diff(g(t),t),'real');
assumeAlso(diff(x(t),t),'real');
assumeAlso(diff(y(t),t),'real');


% Rotation Matrices
R_01 = [cos(Th) 0 -sin(Th); 0 1 0; sin(Th) 0 cos(Th)];
R_12 = [cos(a) -sin(a) 0; sin(a) cos(a) 0; 0 0 1];
R_23 = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
R_34 = [1 0 0; 0 cos(g) -sin(g); 0 sin(g) cos(g)];

R_10 = transpose(R_01);
R_21 = transpose(R_12);
R_32 = transpose(R_23);
R_43 = transpose(R_34);

% Angular Velocities
w0_10 = [0;0;0];
w1_21 = [0;0;aa];
w2_32 = [0;bb;00];
w3_43 = [gg;0;0];

w1_10 = w0_10;
w2_21 = w1_21;
w3_32 = w2_32;
w4_43 = w3_43;

w3_30 = simplify(w3_32 + R_32*w2_21 + R_32*R_21*w1_10);
w4_40 = simplify(w4_43 + R_43*w3_32 + R_43*R_32*w2_21 +R_43*R_32*R_21*w1_10);

%% Constraint Forces
syms l1 l2 % lambda 1 * lambda
 
C1 =  l1*( xx - R*gg*sin(a) - R*bb*cos(a)*cos(b) + R*aa*sin(a)*sin(b));
C2 =  l2*( yy + R*gg*cos(a) - R*aa*cos(a)*sin(b) - R*bb*cos(b)*sin(a));


%% Kinetic Energy
r1_OA = [x;y;R*cos(b)];  % Kinetic energy should be in 
rr1_OA = diff(r1_OA,t);
I4_4A = [1/2*m*R^2 0 0; 0 1/4*m*R^2 0; 0 0 1/4*m*R^2];

Tcom = 1/2*m*rr1_OA'*rr1_OA;
Trot = simplify(1/2*w4_40'*I4_4A*w4_40);
T = Tcom + Trot;

%% Potential Energy
syms grav
r0_OA = R_01*r1_OA;

V = r0_OA'*[0;0;m*grav];


L = simplify(expand(T - V));

% L = subs(L,{a(t),b(t),g(t),x(t),y(t),Th(t)...
%     diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t)},...
%      {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0});

%% Get constraints in correct pfaffian form
CC1 = 1/l1*subs(diff(C1,t),{a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t),...
    diff(a(t), t,t),diff(b(t), t,t), diff(g(t), t,t),diff(x(t), t,t), diff(y(t), t,t),diff(Th(t), t,t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0,aaa_,bbb_,ggg_,xxx_,yyy_,0})
CC2 = 1/l2*subs(diff(C2,t),{a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t),...
    diff(a(t), t,t),diff(b(t), t,t), diff(g(t), t,t),diff(x(t), t,t), diff(y(t), t,t),diff(Th(t), t,t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0,aaa_,bbb_,ggg_,xxx_,yyy_,0})
C = subs(C1 + C2,{a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0})

dW_c = diff(C,aa_)*da + diff(C,bb_)*db + diff(C,gg_)*dg + diff(C,xx_)*dx + diff(C,yy_)*dy;

%% Equations of Motion
%% ERRORS HERE

% This line is edited but not complete
% Get eq in proper form
P1_1 = subs(diff(subs(L,diff(a(t),t),aa_),aa_),aa_,diff(a(t), t));  % correct
P2_1 = subs(diff(subs(diff(subs(L,diff(a(t),t),aa_),aa_),aa_,diff(a(t), t)),t),{a(t),diff(a(t), t),diff(a(t),t,t)},{a_,aa_,aaa_});  % correct
P3_1 = - diff(subs(L,a(t),a_),a_);
eq1 = simplify(P2_1 + P3_1)

P1_2 = subs(diff(subs(L,diff(b(t),t),bb_),bb_),bb_,diff(b(t), t));  % correct
P2_2 = subs(diff(subs(diff(subs(L,diff(b(t),t),bb_),bb_),bb_,diff(b(t), t)),t),{b(t),diff(b(t), t),diff(b(t),t,t)},{b_,bb_,bbb_});  % correct
P3_2 = - diff(subs(L,b(t),b_),b_);
eq2 = simplify(P2_2 + P3_2)

P1_3 = subs(diff(subs(L,diff(g(t),t),gg_),gg_),gg_,diff(g(t), t));  % correct
P2_3 = subs(diff(subs(diff(subs(L,diff(g(t),t),gg_),gg_),gg_,diff(g(t), t)),t),{g(t),diff(g(t), t),diff(g(t),t,t)},{g_,gg_,ggg_});  % correct
P3_3 = - diff(subs(L,g(t),g_),g_);
eq3 = simplify(P2_3 + P3_3)

P1_4 = subs(diff(subs(L,diff(x(t),t),xx_),xx_),xx_,diff(x(t), t));  % correct
P2_4 = subs(diff(subs(diff(subs(L,diff(x(t),t),xx_),xx_),xx_,diff(x(t), t)),t),{x(t),diff(x(t), t),diff(x(t),t,t)},{x_,xx_,xxx_});  % correct
P3_4 = - diff(subs(L,x(t),x_),x_);
eq4 = simplify(P2_4 + P3_4)

P1_5 = subs(diff(subs(L,diff(y(t),t),yy_),yy_),yy_,diff(y(t), t));  % correct
P2_5 = subs(diff(subs(diff(subs(L,diff(y(t),t),yy_),yy_),yy_,diff(y(t), t)),t),{y(t),diff(y(t), t),diff(y(t),t,t)},{y_,yy_,yyy_});  % correct
P3_5 = - diff(subs(L,y(t),y_),y_);
eq5 = simplify(P2_5 + P3_5)


eq1 = subs(eq1- diff(dW_c,da),....
    {a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t),...
    diff(a(t), t,t),diff(b(t), t,t), diff(g(t), t,t),diff(x(t), t,t), diff(y(t), t,t),diff(Th(t), t,t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0,aaa_,bbb_,ggg_,xxx_,yyy_,0})

 %Lines below here are not edited 
eq2 = subs(eq2- diff(dW_c,db),....
    {a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t),...
    diff(a(t), t,t),diff(b(t), t,t), diff(g(t), t,t),diff(x(t), t,t), diff(y(t), t,t),diff(Th(t), t,t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0,aaa_,bbb_,ggg_,xxx_,yyy_,0})
 
eq3 = subs(eq3- diff(dW_c,dg),....
    {a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t),...
    diff(a(t), t,t),diff(b(t), t,t), diff(g(t), t,t),diff(x(t), t,t), diff(y(t), t,t),diff(Th(t), t,t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0,aaa_,bbb_,ggg_,xxx_,yyy_,0})

eq4 = subs(eq4- diff(dW_c,dx),....
    {a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t),...
    diff(a(t), t,t),diff(b(t), t,t), diff(g(t), t,t),diff(x(t), t,t), diff(y(t), t,t),diff(Th(t), t,t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0,aaa_,bbb_,ggg_,xxx_,yyy_,0})
 
eq5 = subs(eq5- diff(dW_c,dy),....
    {a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t),...
    diff(a(t), t,t),diff(b(t), t,t), diff(g(t), t,t),diff(x(t), t,t), diff(y(t), t,t),diff(Th(t), t,t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0,aaa_,bbb_,ggg_,xxx_,yyy_,0})
 
EOM = [eq1; eq2; eq3; eq4; eq5];
EOM = EOM(t); % make symbo
var = [aaa_, bbb_, ggg_, xxx_, yyy_];
 % make symbo

% Match Coefficients
[C, B] = equationsToMatrix(EOM,var);
S =  simplify((inv(C)* B));

% Substitute Values & Extract Equations Of Motion

aaa = S(1)
bbb = S(2)
ggg = S(3)
xxx = S(4)
yyy = S(5)

%% Substiture derivatives into differentiated constraint equations

CC1 = eval(CC1(t))
CC2 = eval(CC2(t))

ConstraintE = [CC1;CC2];
%ConstraintE = ConstraintE(t);
varc = [l1 ,l2];

 
% Match Coefficients
[C, B] = equationsToMatrix(ConstraintE,varc);
Sc =  simplify((inv(C)* B));

l1 = subs(Sc(1),{a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t),...
    diff(a(t), t,t),diff(b(t), t,t), diff(g(t), t,t),diff(x(t), t,t), diff(y(t), t,t),diff(Th(t), t,t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0,aaa_,bbb_,ggg_,xxx_,yyy_,0})
l2 = subs(Sc(2),{a(t),b(t),g(t),x(t),y(t),Th(t)...
    diff(a(t), t),diff(b(t), t), diff(g(t), t),diff(x(t), t), diff(y(t), t),diff(Th(t), t),...
    diff(a(t), t,t),diff(b(t), t,t), diff(g(t), t,t),diff(x(t), t,t), diff(y(t), t,t),diff(Th(t), t,t)},...
     {a_,b_,g_,x_,y_,Th_,aa_,bb_,gg_,xx_,yy_,0,aaa_,bbb_,ggg_,xxx_,yyy_,0})