clc
clear
A=[-209.443458108676,450.728478793718,0,0,0,0;-450.728478793718,-209.443458108676,0,0,0,0;0,0,-2.25680753735172,467.128469627320,0,0;0,0,-467.128469627320,-2.25680753735172,0,0;0,0,0,0,-0.460680331454473,87.5968670596031;0,0,0,0,-87.5968670596031,-0.460680331454473];

B=[-9.64474367821996,2.65660382188316;-7.76311006325836,-5.61613187293423;-0.117927290005266,0.0527206875641465;0.215603654718808,-0.128623447556356;0.0737096142170821,0.0391127484963023;-0.0970049785583168,-0.0786358112919118];

C=[1.88477251668799,-1.14284793707382,5.77230310808448,47.9525873535994,-5.11774457916239,-49.0010462763406];
H=-B(:,1);

dum1=size(A);
dum2=size(B);
dum3=size(H);
dum4=size(C);
n=dum1(1,1);
m=dum2(1,2);
m_H=dum3(1,2);
q=dum4(1,1);

DeltaA=eye(n,n)*0.01;

DeltaB=zeros(n,m);
DeltaB(1,1)=0;

Bp=inv(B'*B)*B';
Gama=eye(n)-B*Bp;


DeltaAm=Bp*DeltaA;
DeltaAu=Gama*DeltaB;
DeltaBm=Bp*DeltaA;
DeltaBu=Gama*DeltaB;
bm=norm(DeltaBm,2);
bu=norm(DeltaBu,2);
a=norm(DeltaA,2);
au=norm(DeltaAu,2);
am=norm(DeltaAm,2);
gm=0;
cm=0;
gama=.90;


function [LME, LMI, OBJ]=AFSMC(XLIST)
[X,Kh,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9]= XLIST(:)
LME=list(X-X')
LMI=list(-([X*A'+A*X+Kh'*B'+B*Kh+(eps1+eps2+eps3+eps4+eps5+eps6+eps7+eps8+eps9)*eye(n,n),zeros(n,m_H),(1+cm)*X'*C',au*X',gm*X'*Bp'*B',am*bu/(1-bm)*X',gm*bu/(1-bm)*X'*Bp',bu/(1-bm)*Kh',gm*X',zeros(n,n+m+n);zeros(m_H,n),-gama^2*eye(m_H,m_H),zeros(m_H,q+n+n+n+m+m+n),H'*Bp'*B',bu/(1-bm)*H'*Bp',H';(1+cm)*C*X,zeros(q,m_H),-eye(q,q),zeros(q,n+n+n+m+m+n+n+m+n);au*X,zeros(n,m_H+q),-eps1*eye(n,n),zeros(n,n+n+m+m+n+n+m+n);gm*B*Bp*X,zeros(n,m_H+q+n),-eps2*eye(n,n),zeros(n,n+m+m+n+n+m+n);am*bu/(1-bm)*X,zeros(n,m_H+q+n+n),-eps4*eye(n,n),zeros(n,m+m+n+n+m+n);gm*bu/(1-bm)*Bp*X,zeros(m,m_H+q+n+n+n),-eps5*eye(m,m),zeros(m,m+n+n+m+n);bu/(1-bm)*Kh,zeros(m,m_H+q+n+n+n+m),-eps7*eye(m,m),zeros(m,n+n+m+n);gm*X,zeros(n,m_H+q+n+n+n+m+m),-eps8*eye(n,n),zeros(n,n+m+n);zeros(n,n),B*Bp*H,zeros(n,q+n+n+n+m+m+n),-eps3*eye(n,n),zeros(n,m+n);zeros(m,n),bu/(1-bm)*Bp*H,zeros(m,q+n+n+n+m+m+n+n),-eps6*eye(m,m),zeros(m,n);zeros(n,n),H,zeros(n,q+n+n+n+m+m+n+n+m),-eps9*eye(n,n)]),X,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9)
OBJ=[]
endfunction
X0=ones(n,n);
Kh0=zeros(m,n);
eps1_0=1;
eps2_0=1;
eps3_0=1;
eps4_0=1;
eps5_0=1;
eps6_0=1;
eps7_0=1;
eps8_0=1;
eps9_0=1;

Init_guess=list(X0,Kh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0);

Ans_LMI=lmisolver(Init_guess, AFSMC);

X0=Ans_LMI(1);
Kh0=Ans_LMI(2);

K0=Kh0*X0^-1;
R0=eye(n,n);
Lh0=zeros(n,q);
eps1_0=1;
eps2_0=1;
eps3_0=1;
eps4_0=1;
eps5_0=1;
eps6_0=1;
eps7_0=1;
eps8_0=1;
eps9_0=1;
eps10_0=1;

gama0=1;

save('sys.dat',A,B,C,H,n,m,q,m_H,a,am,au,bm,bu,gm,cm,Bp);
save('init_vals1.dat',X0,Kh0,K0,R0,Lh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0,eps10_0,gama0);

for ii=1:2
    

clear

load('init_vals1.dat','X0','Kh0','K0','R0','Lh0','eps1_0','eps2_0','eps3_0','eps4_0','eps5_0','eps6_0','eps7_0','eps8_0','eps9_0','eps10_0','gama0');
load('sys.dat','A','B','C','H','n','m','q','m_H','a','am','au','bm','bu','gm','cm','Bp');
mprintf('%f',gama0)
    
function [LME, LMI, OBJ]=AFSMCobsv(XLIST)
[X,Kh,R,Lh,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10,gama]= XLIST(:)
LME=list(X-X',R-R')
LMI=list(-([X'*A'+A*X+Kh'*B'+B*Kh+(eps1+eps2+eps3+eps4+eps5+eps6+eps7)*eye(n,n),-B*K0,H,au*X',gm*X',am*bu/(1-bm)*X',bu/(1-bm)*Kh',gm*X',zeros(n,n+n+m+m+n+n);-K0'*B',A'*R+R*A+C'*Lh'+Lh*C+(eps8*am^2+eps9*au^2)*eye(n,n),R*H,zeros(n,n+n+n+m+n),a*Bp'*B',am*bu/(1-bm)*eye(n,n),bu/(1-bm)*K0',R*B,R,R;H',H'*R,-gama*eye(m_H,m_H),zeros(m_H,n+n+n+m+n+n+n+m+m+n+n);au*X,zeros(n,n+m_H),-eps2*eye(n,n),zeros(n,n+n+m+n+n+n+m+m+n+n);gm*X,zeros(n,n+m_H+n),-eps3*eye(n,n),zeros(n,n+m+n+n+n+m+m+n+n);am*bu/(1-bm)*X,zeros(n,n+m_H+n+n),-eps4*eye(n,n),zeros(n,m+n+n+n+m+m+n+n);bu/(1-bm)*Kh,zeros(m,n+m_H+n+n+n),-eps6*eye(m,m),zeros(m,n+n+n+m+m+n+n);gm*X,zeros(n,n+m_H+n+n+n+m),-eps10*eye(n,n),zeros(n,n+n+m+m+n+n);zeros(n,n),a*B*Bp,zeros(n,m_H+n+n+n+m+n),-eps1*eye(n,n),zeros(n,n+m+m+n+n);zeros(n,n),am*bu/(1-bm)*eye(n,n),zeros(n,m_H+n+n+n+m+n+n),-eps5*eye(n,n),zeros(n,m+m+n+n);zeros(m,n),bu/(1-bm)*K0,zeros(m,m_H+n+n+n+m+n+n+n),-eps7*eye(m,m),zeros(m,m+n+n);zeros(m,n),B'*R',zeros(m,m_H+n+n+n+m+n+n+n+m),-eps8*eye(m,m),zeros(m,n+n);zeros(n,n),R',zeros(n,m_H+n+n+n+m+n+n+n+m+m),-eps9*eye(n,n),zeros(n,n);zeros(n,n),R',zeros(n,m_H+n+n+n+m+n+n+n+m+m+n),-eps10*eye(n,n)]),X,R,eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8,eps9,eps10,gama)
OBJ=gama
endfunction

Init_guess=list(X0,Kh0,R0,Lh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0,eps10_0,gama0);

Ans_LMI=lmisolver(Init_guess,AFSMCobsv);

X0=Ans_LMI(1);
Kh0=Ans_LMI(2);
R0=Ans_LMI(3);
Lh0=Ans_LMI(4);
eps1_0=Ans_LMI(5);
eps2_0=Ans_LMI(6);
eps3_0=Ans_LMI(7);
eps4_0=Ans_LMI(8);
eps5_0=Ans_LMI(9);
eps6_0=Ans_LMI(10);
eps7_0=Ans_LMI(11);
eps8_0=Ans_LMI(12);
eps9_0=Ans_LMI(13);
eps10_0=Ans_LMI(14);
gama0=Ans_LMI(15);
Kind=K0-Kh0*X0^-1;
K0=Kh0*X0^-1;
L0=R0^-1*Lh0;
mprintf('%f',gama0)
save('init_vals1.dat',X0,Kh0,K0,R0,Lh0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,eps6_0,eps7_0,eps8_0,eps9_0,eps10_0,gama0);
end


