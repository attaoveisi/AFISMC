clc
clear all
close all

%% Test System parameters
% A = [0 1 0 0 0 0;-2 0 2 0 0 0;0 0 0 1 0 0;2 0 -4 -3 2 3;0 0 0 0 0 1;0 0 2 3 -2 -3];
% epss = 0.1;
% ksi = epss/(1+epss);
% DeltaA = [0 0 0 0 0 0;2*ksi 0 -2*ksi 0 0 0;0 0 0 0 0 0;-2*ksi 0 4*ksi 3*ksi -2*ksi -3*ksi;0 0 0 0 0 1;0 0 -2*ksi -3*ksi 2*ksi 3*ksi];
% B = [0 1 0 0 0 0;0 0 0 0 0 1]';
% DeltaB = [0 -ksi 0 0 0 0;0 0 0 0 0 -ksi]';
% H = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]';
% C = rand(6,6);
% C = eye(6);
%% Beam System parameters
load System

A = A(1:4,1:4);
DeltaA = A*0;
B = BB(1:4,1);
DeltaB = B*0;
H = BB(1:4,2);
C = C(1,1:4);
%% Required generalization
n = sqrt(numel(A));
m = numel(B)/n;
q = numel(C)/n;

dum1 = size(H'*H);
n_H = dum1(1,1);

Bp = inv(B'*B)*B';
Gama = eye(n)-B*Bp;

%% Perurbation Parameters
DeltaAm = Bp*DeltaA;
DeltaAu = Gama*DeltaB;
DeltaBm = Bp*DeltaA;
DeltaBu = Gama*DeltaB;
cm = 1e-6;
DeltaC = C*cm;
bm = norm(DeltaBm,2);
bu = norm(DeltaBu,2);
a = norm(DeltaA,2);
au = norm(DeltaAu,2);
am = norm(DeltaAm,2);
gm = 1e-6;

%% Initial Control Feedback Gain
% LMI Config.
setlmis([]);

Xc = lmivar(1,[n 1]);
Qc = lmivar(2,[m n]);
eps1c = lmivar(1,[1 0]);
eps2c = lmivar(1,[1 0]);
eps3c = lmivar(1,[1 0]);
eps4c = lmivar(1,[1 0]);
eps5c = lmivar(1,[1 0]);
eps6c = lmivar(1,[1 0]);
eps7c = lmivar(1,[1 0]);
eps8c = lmivar(1,[1 0]);
eps9c = lmivar(1,[1 0]);
beta = lmivar(1,[1 0]);

%LMI terms
lmiterm([1 1 1 Xc],A,1,'s'); % LMI #1: A*Xc+Xc*A'
lmiterm([1 1 1 Qc],B,1,'s'); % LMI #1: B*Qc+Xc'*B'
lmiterm([1 1 1 eps1c],1,1); % LMI #1: eps1c*I
lmiterm([1 1 1 eps2c],1,1); % LMI #1: eps2c*I
lmiterm([1 1 1 eps3c],1,1); % LMI #1: eps3c*I
lmiterm([1 1 1 eps4c],1,1); % LMI #1: eps4c*I
lmiterm([1 1 1 eps5c],1,1); % LMI #1: eps5c*I
lmiterm([1 1 1 eps6c],1,1); % LMI #1: eps6c*I
lmiterm([1 1 1 eps7c],1,1); % LMI #1: eps7c*I
lmiterm([1 1 1 eps8c],1,1); % LMI #1: eps8c*I
lmiterm([1 1 1 eps9c],1,1); % LMI #1: eps9c*I
lmiterm([1 1 2 0],0); % LMI #1: 0
lmiterm([1 1 3 -Xc],(1+cm),C'); % LMI #1: (1+cm)*Xc'*C'
lmiterm([1 1 4 -Xc],au,1); % LMI #1: au*Xc'
lmiterm([1 1 5 -Xc],gm,Bp'*B'); % LMI #1: gm*Xc'*Bp'*B'
lmiterm([1 1 6 -Xc],am*bu/(1-bm),1); % LMI #1: am*bu/(1-bm)*Xc'
lmiterm([1 1 7 -Xc],gm*bu/(1-bm),Bp'); % LMI #1: gm*bu/(1-bm)*Xc'*Bp'
lmiterm([1 1 8 -Qc],bu/(1-bm),1); % LMI #1: bu/(1-bm)*Qc'
lmiterm([1 1 9 -Xc],gm,1); % LMI #1: gm*Xc'
lmiterm([1 2 2 beta],-1,1); % LMI #1: -gama^2*I --->  beta*I (beta > 0) for minimization
lmiterm([1 2 3 0],0); % LMI #1: 0
lmiterm([1 2 4 0],0); % LMI #1: 0
lmiterm([1 2 5 0],0); % LMI #1: 0
lmiterm([1 2 6 0],0); % LMI #1: 0
lmiterm([1 2 7 0],0); % LMI #1: 0
lmiterm([1 2 8 0],0); % LMI #1: 0
lmiterm([1 2 9 0],0); % LMI #1: 0
lmiterm([1 2 10 0],H'*Bp'*B'); % LMI #1: H'*Bp'*B'
lmiterm([1 2 11 0],bu/(1-bm)*H'*Bp'); % LMI #1: bu/(1-bm)*H'*Bp'
lmiterm([1 2 12 0],H'); % LMI #1: H'
lmiterm([1 3 3 0],-1); % LMI #1: -I
lmiterm([1 4 4 eps1c],-1,1); % LMI #1: -eps1c*I
lmiterm([1 5 5 eps2c],-1,1); % LMI #1: -eps1c*I
lmiterm([1 6 6 eps4c],-1,1); % LMI #1: -eps1c*I
lmiterm([1 7 7 eps5c],-1,1); % LMI #1: -eps1c*I
lmiterm([1 8 8 eps7c],-1,1); % LMI #1: -eps1c*I
lmiterm([1 9 9 eps8c],-1,1); % LMI #1: -eps1c*I
lmiterm([1 10 10 eps3c],-1,1); % LMI #1: -eps1c*I
lmiterm([1 11 11 eps6c],-1,1); % LMI #1: -eps1c*I
lmiterm([1 12 12 eps9c],-1,1); % LMI #1: -eps1c*I

lmiterm([-2 1 1 Xc],1,1); % LMI #2: Xc>0

lmiterm([-3 1 1 eps1c],1,1); % LMI #3: eps1c>0
lmiterm([-4 1 1 eps2c],1,1); % LMI #4: eps2c>0
lmiterm([-5 1 1 eps3c],1,1); % LMI #5: eps3c>0
lmiterm([-6 1 1 eps4c],1,1); % LMI #6: eps4c>0
lmiterm([-7 1 1 eps5c],1,1); % LMI #7: eps5c>0
lmiterm([-8 1 1 eps6c],1,1); % LMI #8: eps6c>0
lmiterm([-9 1 1 eps7c],1,1); % LMI #9: eps7c>0
lmiterm([-10 1 1 eps8c],1,1); % LMI #10: eps8c>0
lmiterm([-11 1 1 eps9c],1,1); % LMI #11: eps9c>0


LMISYS1 = getlmis;

dumc1 = n*(n+1)/2+m*n+9+1;
C_mat1 = zeros(dumc1,1);
C_mat1(end,1) = 1;

% Finding otimal Solution for beta
[~,xc] = mincx(LMISYS1,C_mat1,[0 1000 -1 100 0]);


Xc = dec2mat(LMISYS1,xc,1);
Qc = dec2mat(LMISYS1,xc,2);
eps1c = dec2mat(LMISYS1,xc,3);
eps2c = dec2mat(LMISYS1,xc,4);
eps3c = dec2mat(LMISYS1,xc,5);
eps4c = dec2mat(LMISYS1,xc,6);
eps5c = dec2mat(LMISYS1,xc,7);
eps6c = dec2mat(LMISYS1,xc,8);
eps7c = dec2mat(LMISYS1,xc,9);
eps8c = dec2mat(LMISYS1,xc,10);
eps9c = dec2mat(LMISYS1,xc,11);

Kc = Qc*inv(Xc);

[eig(A+DeltaA)]
[eig(A+DeltaA+(B+DeltaB)*Kc)]

%% Initial observer Feedback Gain with -alfa^2*W'*W + e'*e
% LMI Config.
setlmis([]);

Ro = lmivar(1,[n 1]);
Qo = lmivar(2,[n q]);
eps1o = lmivar(1,[1 0]);
eps2o = lmivar(1,[1 0]);
eps3o = lmivar(1,[1 0]);
alfa = lmivar(1,[1 0]);

%LMI terms
lmiterm([1 1 1 Ro],1,A,'s'); % LMI #1: R*A+A'*R
lmiterm([1 1 1 Qo],1,C,'s'); % LMI #1: Qo*C+C'*Qo'
lmiterm([1 1 1 eps1o],1,1); % LMI #1: eps1o*I
lmiterm([1 1 1 eps2o],1,1); % LMI #1: eps2o*I
lmiterm([1 1 1 eps3o],1,1); % LMI #1: eps3o*I
lmiterm([1 1 1 0],1); % LMI #1: I

lmiterm([1 1 2 0],0); % LMI #1: 0
lmiterm([1 1 3 0],0); % LMI #1: 0
lmiterm([1 1 4 -Ro],a,1); % LMI #1: a*R'
lmiterm([1 1 5 -Qo],cm*C',1); % LMI #1: cm*C'*Qo'
lmiterm([1 2 2 alfa],-1,1); % LMI #1: -alfa^2*I  -->>>  -alfa*I  (alfa > 0) for minimization
lmiterm([1 2 3 -Ro],H',1); % LMI #1: H'*Ro'
lmiterm([1 3 3 eps2o],-1,1); % LMI #1: -eps2o*I
lmiterm([1 4 4 eps1o],-1,1); % LMI #1: -eps1o*I
lmiterm([1 5 5 eps3o],-1,1); % LMI #1: -eps3o*I

lmiterm([-2 1 1 Ro],1,1); % LMI #2: Ro>0

lmiterm([-3 1 1 eps1o],1,1); % LMI #3: eps1o>0
lmiterm([-4 1 1 eps2o],1,1); % LMI #4: eps2o>0
lmiterm([-5 1 1 eps3o],1,1); % LMI #5: eps3o>0

lmiterm([-6 1 1 alfa],1,1); % LMI #5: alfa>0

LMISYS2 = getlmis;

dumc2 = n*(n+1)/2+n*q+3+1;
C_mat2 = zeros(dumc2,1);
C_mat2(end,1) = 1;

% Finding otimal Solution for alfa
[~,xo] = mincx(LMISYS2,C_mat2,[0 1000 -1 100 0]);

Ro = dec2mat(LMISYS2,xo,1);
Qo = dec2mat(LMISYS2,xo,2);
eps1o = dec2mat(LMISYS2,xo,3);
eps2o = dec2mat(LMISYS2,xo,4);
eps3o = dec2mat(LMISYS2,xo,5);
alfa = dec2mat(LMISYS2,xo,6);

Lo = inv(Ro)*Qo;
%% Hybrid system control gain Gama
gama = 1;

%% Iteration Parameters
imax = 2;

%% LMI Solver
for i = 1:imax
setlmis([]);

% LMI Variables
Q1 = lmivar(2,[m n]);
X = lmivar(1,[n 1]);
Q2 = lmivar(2,[n q]);
L = lmivar(2,[n q]);
R = lmivar(1,[n 1]);
K = lmivar(2,[m n]);
eps1 = lmivar(1,[1 0]);
eps2 = lmivar(1,[1 0]);
eps3 = lmivar(1,[1 0]);
eps4 = lmivar(1,[1 0]);
eps5 = lmivar(1,[1 0]);
eps6 = lmivar(1,[1 0]);
eps7 = lmivar(1,[1 0]);
eps8 = lmivar(1,[1 0]);
eps9 = lmivar(1,[1 0]);
eps10 = lmivar(1,[1 0]);
eps11 = lmivar(1,[1 0]);
eps12 = lmivar(1,[1 0]);
eps13 = lmivar(1,[1 0]);
eps14 = lmivar(1,[1 0]);
eps15 = lmivar(1,[1 0]);
eps16 = lmivar(1,[1 0]);
eps17 = lmivar(1,[1 0]);
eps18 = lmivar(1,[1 0]);

% LMI terms
lmiterm([1 1 1 X],A,1,'s'); % LMI #1: A*X+X*A'
lmiterm([1 1 1 Q1],B,1,'s'); % LMI #1: B*Q1+Q1'*B'
lmiterm([1 1 1 eps1],1,1); % LMI #1: eps1*I
lmiterm([1 1 1 eps2],1,1); % LMI #1: eps2*I
lmiterm([1 1 1 eps3],1,1); % LMI #1: eps3*I
lmiterm([1 1 1 eps4],1,1); % LMI #1: eps4*I
lmiterm([1 1 1 eps5],1,1); % LMI #1: eps5*I
lmiterm([1 1 1 eps6],1,1); % LMI #1: eps6*I
lmiterm([1 1 1 eps7],1,1); % LMI #1: eps7*I
lmiterm([1 1 1 eps8],1,1); % LMI #1: eps8*I
lmiterm([1 1 1 eps9],1,1); % LMI #1: eps9*I
lmiterm([1 1 1 eps10],1,1); % LMI #1: eps10*I
lmiterm([1 1 1 eps11],1,1); % LMI #1: eps11*I
lmiterm([1 1 1 eps12],1,1); % LMI #1: eps12*I
lmiterm([1 1 1 eps13],1,1); % LMI #1: eps13*I
lmiterm([1 1 2 0],0); % LMI #1: 0
lmiterm([1 1 3 0],0); % LMI #1: 0
lmiterm([1 1 4 -X],au,1); % LMI #1: au*X'
lmiterm([1 1 5 -X],am*bu/(1-bm),1); % LMI #1: am*bu/(1-bm)*X'
lmiterm([1 1 6 -Q1],bu/(1-bm),1); % LMI #1: bu/(1-bm)*Q1'
lmiterm([1 1 7 -X],gm,1); % LMI #1: gm*X'

%%%% $$$ %%%% NONLINEAR TERM !!!
lmiterm([1 1 8 -X],gm,Ro); % LMI #1: am*bu/(1-bm)*X'
%%%% $$$ %%%% NONLINEAR TERM !!!

lmiterm([1 1 9 -X],(1+cm),C'); % LMI #1: (1+cm)*X'*C'
lmiterm([1 2 2 R],1,A,'s'); % LMI #1: R*A+A'*R
lmiterm([1 2 2 Q2],-1,C,'s'); % LMI #1: R*A+A'*R
lmiterm([1 2 2 eps14],1,1); % LMI #1: eps14*I
lmiterm([1 2 2 eps15],1,1); % LMI #1: eps15*I
lmiterm([1 2 2 eps16],1,1); % LMI #1: eps16*I
lmiterm([1 2 2 eps17],1,1); % LMI #1: eps17*I
lmiterm([1 2 2 eps18],1,1); % LMI #1: eps18*I
lmiterm([1 2 2 0],1); % LMI #1: I
lmiterm([1 2 3 0],0); % LMI #1: 0
lmiterm([1 2 4 0],0); % LMI #1: 0
lmiterm([1 2 5 0],0); % LMI #1: 0
lmiterm([1 2 3 0],0); % LMI #1: 0
lmiterm([1 2 6 0],0); % LMI #1: 0
lmiterm([1 2 7 0],0); % LMI #1: 0
lmiterm([1 2 8 0],0); % LMI #1: 0
lmiterm([1 2 9 0],0); % LMI #1: 0
lmiterm([1 2 10 0],am*bu/(1-bm)); % LMI #1: am*bu/(1-bm)*I
lmiterm([1 2 11 -K],bu/(1-bm),1); % LMI #1: bu/(1-bm)*K'
lmiterm([1 2 12 -L],bu/(1-bm)*C',Bp'); % LMI #1: bu/(1-bm)*C'*L'*Bp'
lmiterm([1 2 13 -L],cm*bu/(1-bm)*C',Bp'); % LMI #1: cm*bu/(1-bm)*C'*L'*Bp'
lmiterm([1 2 14 0],a*Bp'*B'); % LMI #1: a*Bp'*B'
lmiterm([1 2 15 -K],1,B'); % LMI #1: K'*B'
lmiterm([1 2 16 -L],C',Bp'*B'); % LMI #1: C'*L'*Bp'*B'
lmiterm([1 2 17 -L],cm*C',Bp'*B'); % LMI #1: cm*C'*L'*Bp'*B'
lmiterm([1 2 18 R],au,1); % LMI #1: au*R
lmiterm([1 2 19 R],a*Bp'*B',1); % LMI #1: a*Bp'*B'*R
lmiterm([1 2 20 -Q2],cm*C',1); % LMI #1: cm*C'*Q2'
lmiterm([1 3 3 0],-gama^2); % LMI #1: -gama^2*I
lmiterm([1 3 4 0],0); % LMI #1: 0
lmiterm([1 3 5 0],0); % LMI #1: 0
lmiterm([1 3 6 0],0); % LMI #1: 0
lmiterm([1 3 7 0],0); % LMI #1: 0
lmiterm([1 3 8 0],0); % LMI #1: 0
lmiterm([1 3 9 0],0); % LMI #1: 0
lmiterm([1 3 10 0],0); % LMI #1: 0
lmiterm([1 3 11 0],0); % LMI #1: 0
lmiterm([1 3 12 0],0); % LMI #1: 0
lmiterm([1 3 13 0],0); % LMI #1: 0
lmiterm([1 3 14 0],0); % LMI #1: 0
lmiterm([1 3 15 0],0); % LMI #1: 0
lmiterm([1 3 16 0],0); % LMI #1: 0
lmiterm([1 3 17 0],0); % LMI #1: 0
lmiterm([1 3 18 0],0); % LMI #1: 0
lmiterm([1 3 19 0],0); % LMI #1: 0
lmiterm([1 3 20 0],0); % LMI #1: 0
lmiterm([1 3 21 0],H'); % LMI #1: H'
lmiterm([1 3 22 R],H',1); % LMI #1: H'*R
lmiterm([1 4 4 eps1],-1,1); % LMI #1: -eps1*I
lmiterm([1 5 5 eps2],-1,1); % LMI #1: -eps2*I
lmiterm([1 6 6 eps4],-1,1); % LMI #1: -eps4*I
lmiterm([1 7 7 eps12],-1,1); % LMI #1: -eps12*I
lmiterm([1 8 8 eps17],-1,1); % LMI #1: -eps17*I
lmiterm([1 9 9 0],-1); % LMI #1: -I
lmiterm([1 10 10 eps3],-1,1); % LMI #1: -eps3*I
lmiterm([1 11 11 eps5],-1,1); % LMI #1: -eps5*I
lmiterm([1 12 12 eps6],-1,1); % LMI #1: -eps6*I
lmiterm([1 13 13 eps7],-1,1); % LMI #1: -eps7*I
lmiterm([1 14 14 eps8],-1,1); % LMI #1: -eps8*I
lmiterm([1 15 15 eps9],-1,1); % LMI #1: -eps9*I
lmiterm([1 16 16 eps10],-1,1); % LMI #1: -eps10*I
lmiterm([1 17 17 eps11],-1,1); % LMI #1: -eps11*I
lmiterm([1 18 18 eps14],-1,1); % LMI #1: -eps14*I
lmiterm([1 19 19 eps15],-1,1); % LMI #1: -eps15*I
lmiterm([1 20 20 eps16],-1,1); % LMI #1: -eps16*I
lmiterm([1 21 21 eps13],-1,1); % LMI #1: -eps13*I
lmiterm([1 22 22 eps18],-1,1); % LMI #1: -eps18*I

LMISYS3 = getlmis;

% cvector = [zeros(1,n^2) ];

% Finding Feasible Solution 
[~,xgen] = feasp(LMISYS3,[0 1000 -1 99 0]);

R = dec2mat(LMISYS3,xgen,1);
X = dec2mat(LMISYS3,xgen,2);
L = dec2mat(LMISYS3,xgen,3);
K = dec2mat(LMISYS3,xgen,4);
Q1 = dec2mat(LMISYS3,xgen,5);
Q2 = dec2mat(LMISYS3,xgen,6);

end
% % Finding Opimal Solution 
% target = 1e-10;
% [copt,xopt] = mincx(lmisys,cvector,[1e-2 1000 -1 100 0],xinit,target);

