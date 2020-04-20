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

%% Iteration Parameters
imax = 30;

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

%% Controller Parameters
gama = beta;

%% CVX Solver
for i = 1:imax
    
    cvx_begin SDP

        variable Q1(m,n) ;

        variable Q2(n,q) ;
        
        variable Q3(n,n) symmetric; 
        Q3 == semidefinite(n);

        variable R(n,n) symmetric; 
        R == semidefinite(n);

        variable L(n,q) ;

        variable K(m,n);

        variable X(n,n) symmetric; 
        X == semidefinite(n);

        variable eps1;
        eps1 == semidefinite(1);
        variable eps2;
        eps2 == semidefinite(1);
        variable eps3;
        eps3 == semidefinite(1);
        variable eps4;
        eps4 == semidefinite(1);
        variable eps5;
        eps5 == semidefinite(1);
        variable eps6;
        eps6 == semidefinite(1);
        variable eps7;
        eps7 == semidefinite(1);
        variable eps8;
        eps8 == semidefinite(1);
        variable eps9;
        eps9 == semidefinite(1);
        variable eps10;
        eps10 == semidefinite(1);
        variable eps11;
        eps11 == semidefinite(1);
        variable eps12;
        eps12 == semidefinite(1);
        variable eps13;
        eps13 == semidefinite(1);
        variable eps14;
        eps14 == semidefinite(1);
        variable eps15;
        eps15 == semidefinite(1);
        variable eps16;
        eps16 == semidefinite(1);
        variable eps17;
        eps17 == semidefinite(1);
        variable eps18;
        eps18 == semidefinite(1);
        

        minimize(norm(Q1-Kc*X,2)+norm(Q2-Ro*L,2)+norm(Q3-Ro*X,2));

        subject to

            Pii = [X*A'+A*X+Q1'*B'+B*Q1+(eps1+eps2+eps3+eps4+eps5+eps6+eps7+eps8+eps9+eps10+eps11+eps12+eps13)*eye(n) zeros(n,n) zeros(n,n_H) au*X' am*bu/(1-bm)*X' bu/(1-bm)*Q1' gm*X' gm*Q3' (1+cm)*X'*C' zeros(n,n+m+m+m+n+n+n+n+n+n+n) zeros(n,n+n);...
                    zeros(n,n) A'*R+R*A-Q2*C-C'*Q2'+(eps14+eps15+eps16+eps17+eps18+1)*eye(n) zeros(n,n_H) zeros(n,n+n+m+n+n+q) am*bu/(1-bm)*eye(n) bu/(1-bm)*K' bu/(1-bm)*C'*L'*Bp' cm*bu/(1-bm)*C'*L'*Bp' a*Bp'*B' K'*B' C'*L'*Bp'*B' cm*C'*L'*Bp'*B' au*R a*Bp'*B'*R cm*C'*Q2' zeros(n,n+n);...
                    zeros(n_H,n) zeros(n_H,n) -gama^2*eye(n_H) zeros(n_H,n+n+m+n+n+q) zeros(n_H,n+m+m+m+n+n+n+n+n+n+n) H' H'*R;...
                    au*X zeros(n,n+n_H) -eps1*eye(n) zeros(n,n+m+n+n+q) zeros(n,n+m+m+m+n+n+n+n+n+n+n) zeros(n,n+n);...
                    am*bu/(1-bm)*X zeros(n,n+n_H+n) -eps2*eye(n) zeros(n,m+n+n+q) zeros(n,n+m+m+m+n+n+n+n+n+n+n) zeros(n,n+n);...
                    bu/(1-bm)*Q1 zeros(m,n+n_H+n+n) -eps4*eye(m) zeros(m,n+n+q) zeros(m,n+m+m+m+n+n+n+n+n+n+n) zeros(m,n+n);...
                    gm*X zeros(n,n+n_H+n+n+m) -eps12*eye(n) zeros(n,n+q) zeros(n,n+m+m+m+n+n+n+n+n+n+n) zeros(n,n+n);...
                    gm*Q3 zeros(n,n+n_H+n+n+m+n) -eps17*eye(n) zeros(n,q) zeros(n,n+m+m+m+n+n+n+n+n+n+n) zeros(n,n+n);...
                    (1+cm)*C*X zeros(q,n+n_H+n+n+m+n+n) -eye(q) zeros(q,n+m+m+m+n+n+n+n+n+n+n) zeros(q,n+n);...
                    zeros(n) am*bu/(1-bm)*eye(n) zeros(n,n_H+n+n+m+n+n+q) -eps3*eye(n) zeros(n,m+m+m+n+n+n+n+n+n+n) zeros(n,n+n);...
                    zeros(m,n) bu/(1-bm)*K zeros(m,n_H+n+n+m+n+n+q+n) -eps5*eye(m) zeros(m,m+m+n+n+n+n+n+n+n) zeros(m,n+n);...
                    zeros(m,n) bu/(1-bm)*Bp*L*C zeros(m,n_H+n+n+m+n+n+q+n+m) -eps6*eye(m) zeros(m,m+n+n+n+n+n+n+n) zeros(m,n+n);...
                    zeros(m,n) bu*cm/(1-bm)*Bp*L*C zeros(m,n_H+n+n+m+n+n+q+n+m+m) -eps7*eye(m) zeros(m,n+n+n+n+n+n+n) zeros(m,n+n);...
                    zeros(n) a*B*Bp zeros(n,n_H+n+n+m+n+n+q+n+m+m+m) -eps8*eye(n) zeros(n,n+n+n+n+n+n) zeros(n,n+n);...
                    zeros(n) B*K zeros(n,n_H+n+n+m+n+n+q+n+m+m+m+n) -eps9*eye(n) zeros(n,n+n+n+n+n) zeros(n,n+n);...
                    zeros(n) B*Bp*L*C zeros(n,n_H+n+n+m+n+n+q+n+m+m+m+n+n) -eps10*eye(n) zeros(n,n+n+n+n) zeros(n,n+n);...
                    zeros(n) cm*B*Bp*L*C zeros(n,n_H+n+n+m+n+n+q+n+m+m+m+n+n+n) -eps11*eye(n) zeros(n,n+n+n) zeros(n,n+n);...
                    zeros(n) au*R zeros(n,n_H+n+n+m+n+n+q+n+m+m+m+n+n+n+n) -eps14*eye(n) zeros(n,n+n) zeros(n,n+n);...
                    zeros(n) a*R*B*Bp zeros(n,n_H+n+n+m+n+n+q+n+m+m+m+n+n+n+n+n) -eps15*eye(n) zeros(n,n) zeros(n,n+n);...
                    zeros(n) cm*Q2*C zeros(n,n_H+n+n+m+n+n+q+n+m+m+m+n+n+n+n+n+n) -eps16*eye(n) zeros(n,n+n);...
                    zeros(n,n+n) H zeros(n,n+n+m+n+n+q+n+m+m+m+n+n+n+n+n+n+n) -eps13*eye(n) zeros(n,n);...
                    zeros(n,n+n) R*H zeros(n,n+n+m+n+n+q+n+m+m+m+n+n+n+n+n+n+n+n) -eps18*eye(n)];
    %                 dum2 = n+n+n_H+n+n+m+n+n+q+n+m+m+m+n+n+n+n+n+n+n+n+n;


             0 >= Pii;
             
    cvx_end
    
    Kc = K;
    [eig(A+DeltaA)]
    [eig(A+DeltaA+(B+DeltaB)*Kc)]
    Ro = R;

end