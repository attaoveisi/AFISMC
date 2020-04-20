clc
close all
clear all

m1 = 450*10^3;
m2 = 345*10^3;
m3 = 345*10^3;
m4 = 345*10^3;
c1 = 26.170*10^3;
c2 = 490*10^3;
c3 = 467*10^3;
c4 = 410*10^3;
k1 = 18.050*10^3;
k2 = 340000*10^3;
k3 = 326000*10^3;
k4 = 280000*10^3;

hb0 = 0;
Gama0 = 27;
pb0 = 0;
mu = 0;


kmax = 1000;
imax = 80;

deltaGama = 0.025;
deltapb = 0.0010;
deltahb = 0.00001;

hb = hb0;
Gama = Gama0;
pb = pb0;
EA = pb*eye(8,8);
C = eye(8,8);
D = zeros(8,1);
Dw = zeros(8,2);

Ms = blkdiag(m1,m2,m3,m4);
Cs = [c1+c2 -c2 0 0
    -c2 c2+c3 -c3 0
    0 -c3 c3+c4 -c4
    0 0 -c4 c4];
Ks = [k1+k2 -k2 0 0
    -k2 k2+k3 -k3 0
    0 -k3 k3+k4 -k4
    0 0 -k4 k4];

A = [zeros(4,4) eye(4,4)
    -inv(Ms)*Ks -inv(Ms)*Cs];
Bh = [0 0 0 0 -1/m1 1/m2 0 0]';
Bw = [0 0 0 0 k1/m1 0 0 0
    0 0 0 0 c1/m1 0 0 0]';

G = A;

for i=1:1
    
% Feasibility Problem


Gama = Gama-deltaGama;      %updating the Gama
hb = hb+deltahb;
pb = pb+deltapb;

 setlmis([]);
  X = lmivar(1,[8 1]);
  Qb = lmivar(1,[8 1]);
  Wb = lmivar(1,[8 1]);
  Zb = lmivar(1,[8 1]);
  T = lmivar(1,[8 1]);
  
  N1b = lmivar(2,[8 8]);
  N2b = lmivar(2,[8 8]);
  S1b = lmivar(2,[8 8]);
  S2b = lmivar(2,[8 8]);
  
  eps = lmivar(1,[1 0]);
  
  Lb = lmivar(2,[1 8]);
  
  R = lmivar(1,[8 1]);
  Rb = lmivar(1,[8 1]);
  Xb = lmivar(1,[8 1]);
  Tb = lmivar(1,[8 1]);
  
  % LMI 1
  % 1
  
  lmiterm([1 1 1 X],A,1,'s');                    
  lmiterm([1 1 1 N1b],1,1,'s');                      
  lmiterm([1 1 1 Qb],1,1);                        
  lmiterm([1 1 1 Wb],1,1);                        
  lmiterm([1 1 1 eps],1,G*G');
  
  lmiterm([1 1 2 Lb],Bh,1);
  lmiterm([1 1 2 N1b],-1,1);
  lmiterm([1 1 2 N2b],1,1,'s');
  lmiterm([1 1 2 N2b],-1,1);
  lmiterm([1 1 2 S1b],1,1);
  
  lmiterm([1 1 3 S1b],-1,1);

  lmiterm([1 1 4 0],Bw);
  
  lmiterm([1 1 5 0],0);
  
  lmiterm([1 1 6 N1b],hb,1);
  
  lmiterm([1 1 7 S1b],hb,1);
  
  lmiterm([1 1 8 eps],1,G*G');
  lmiterm([1 1 8 X],hb,A');
  
  lmiterm([1 1 9 0],0);
  
  lmiterm([1 1 10 X],1,C');
  lmiterm([1 1 10 Lb],D,1,'s');
  lmiterm([1 1 10 Lb],-D,1);
  
  lmiterm([1 1 11 X],1,EA');
  % 2
  
  lmiterm([1 2 2 Qb],-(1-mu),1);
  lmiterm([1 2 2 S2b],1,1,'s');
  lmiterm([1 2 2 N2b],-1,1,'s');
  
  lmiterm([1 2 3 S1b],-1,1);
  
  lmiterm([1 2 4 0],0);
  
  lmiterm([1 2 5 0],0);
  
  lmiterm([1 2 6 N2b],hb,1);
  
  lmiterm([1 2 7 S2b],hb,1);
  
  lmiterm([1 2 8 Lb],Bh,hb','s');
  lmiterm([1 2 8 Lb],-Bh,hb');
  
  lmiterm([1 2 9 0],0);
  
  lmiterm([1 2 10 0],0);
  
  lmiterm([1 2 11 0],0);
  % 3
  
  lmiterm([1 3 3 Wb],-1,1);
  
  lmiterm([1 3 4 0],0);
  
  lmiterm([1 3 5 0],0);
  
  lmiterm([1 3 6 0],0);
  
  lmiterm([1 3 7 0],0);
  
  lmiterm([1 3 8 0],0);
  
  lmiterm([1 3 9 0],0);
  
  lmiterm([1 3 10 0],0);
  
  lmiterm([1 3 11 0],0);
  % 4
  
  lmiterm([1 4 4 0],-Gama^2);
  
  lmiterm([1 4 5 0],0);
  
  lmiterm([1 4 6 0],0);
  
  lmiterm([1 4 7 0],0);
  
  lmiterm([1 4 8 0],hb*Bw');
  
  lmiterm([1 4 9 0],0);
  
  lmiterm([1 4 10 0],0);
  
  lmiterm([1 4 11 0],0);
  % 5s
  
  lmiterm([1 5 5 Zb],-hb,1);
  
  lmiterm([1 5 6 0],0);
  
  lmiterm([1 5 7 0],0);
  
  lmiterm([1 5 8 0],0);
  
  lmiterm([1 5 9 Zb],1,1);
  
  lmiterm([1 5 10 0],0);
  
  lmiterm([1 5 11 0],0);
  % 6
  
  lmiterm([1 6 6 Zb],-1,1);
  
  lmiterm([1 6 7 0],0);
  
  lmiterm([1 6 8 0],0);
  
  lmiterm([1 6 9 0],0);
  
  lmiterm([1 6 10 0],0);
  
  lmiterm([1 6 11 0],0);
  % 7
  
  lmiterm([1 7 7 Zb],-1,1);
  
  lmiterm([1 7 8 0],0);
  
  lmiterm([1 7 9 0],0);
  
  lmiterm([1 7 10 0],0);
  
  lmiterm([1 7 11 0],0);
  % 8
  
  lmiterm([1 8 8 eps],1,G*G');
  lmiterm([1 8 8 T],1,1);
  
  lmiterm([1 8 9 0],0);
  
  lmiterm([1 8 10 0],0);
  
  lmiterm([1 8 11 X],1,EA);
  % 9
  
  lmiterm([1 9 9 R],-1,1);
  
  lmiterm([1 9 10 0],0);
  
  lmiterm([1 9 11 0],0);
  % 10
  
  lmiterm([1 10 10 0],-1);
  
  lmiterm([1 10 11 0],0);
  % 11
  
  lmiterm([1 11 11 eps],-1,1);
  
  
  % LMI 2
  
  lmiterm([2 1 1 Rb],1,1);
  
  lmiterm([2 1 2 Xb],1,1);
  
  lmiterm([2 2 2 Tb],1,1);
  % LMI 3
  
  lmiterm([3 1 1 Rb],1,1);
  
  lmiterm([3 1 2 0],1);
  
  lmiterm([3 2 2 R],1,1);
  % LMI 4
  
  lmiterm([4 1 1 Xb],1,1);
  
  lmiterm([4 1 2 0],1);
  
  lmiterm([4 2 2 X],1,1);
  % LMI 5
  
  lmiterm([5 1 1 Tb],1,1);
  
  lmiterm([5 1 2 0],1);
  
  lmiterm([5 2 2 T],1,1);
  %LMI 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 
  
  lmiterm([6 1 1 R],1,1);
  lmiterm([7 1 1 T],1,1);
  lmiterm([8 1 1 X],1,1);
  lmiterm([9 1 1 Qb],1,1);
  lmiterm([10 1 1 Wb],1,1);
  lmiterm([11 1 1 Zb],1,1);
  lmiterm([12 1 1 Rb],1,1);
  lmiterm([13 1 1 Xb],1,1);
  lmiterm([14 1 1 Tb],1,1);
    
  LMISYS=getlmis;
  
  [t,x] = feasp(LMISYS);
  
  
  Xopt = dec2mat(LMISYS,x,1);
  Qbopt = dec2mat(LMISYS,x,2);
  Wbopt = dec2mat(LMISYS,x,3);
  Zbopt = dec2mat(LMISYS,x,4);
  Topt = dec2mat(LMISYS,x,5);
  N1bopt = dec2mat(LMISYS,x,6);
  N2bopt = dec2mat(LMISYS,x,7);
  S1bopt = dec2mat(LMISYS,x,8);
  S2bopt = dec2mat(LMISYS,x,9);
  epsopt = dec2mat(LMISYS,x,10);
  Lbopt = dec2mat(LMISYS,x,11);
  Ropt = dec2mat(LMISYS,x,12);
  Rbopt = dec2mat(LMISYS,x,12);
  Xbopt = dec2mat(LMISYS,x,13);
  Tbopt = dec2mat(LMISYS,x,14);
  

for k = 1:1;
    
   
cvx_begin SDP
    variable X(8,8) symmetric;
    X == semidefinite(8);
    variable Qb(8,8) symmetric;
    Qb == semidefinite(8);

    variable Wb(8,8) symmetric; 
    Wb == semidefinite(8);

    variable Zb(8,8) symmetric; 
    Zb == semidefinite(8);
    
    variable T(8,8) symmetric;
    T == semidefinite(8);
    
    variable N1b(8,8);
    variable N2b(8,8);
    variable S1b(8,8);
    variable S2b(8,8);
    variable eps;
    eps == semidefinite(1);
    
    variable Lb(1,8);
    
    variable R(8,8) symmetric;
    R == semidefinite(8);
    variable Rb(8,8) symmetric;
    Rb == semidefinite(8);
    variable Xb(8,8) symmetric;
    Xb == semidefinite(8);
    variable Tb(8,8) symmetric;
    Tb == semidefinite(8);
    
    minimize(trace(Rb*Ropt+Xb*Xopt+Tb*Topt+Rbopt*R+Xbopt*X+Tbopt*T));
    
    subject to
        Sib11 = A*X+X*A'+N1b+N1b'+Qb+Wb+eps*(G*G');
        Sib12 = Bh*Lb-N1b+N2b'+S1b;
        Sib18 = eps*(G*G')+hb*X*A';
        Sib22 = -(1-mu)*Qb+S2b+S2b'-N2b-N2b';
        Sib28 = hb*Lb'*Bh';
        Sib88 = eps*(G*G')-T;
        Sib111 = X*EA';
        Sib811 = X*EA';
        
        
        Sib = [Sib11 Sib12 -S1b Bw zeros(8,8) hb*N1b hb*S1b Sib18 zeros(8,8) X*C'+Lb'*D' Sib111
            Sib12' Sib22 -S1b zeros(8,2) zeros(8,8) hb*N2b hb*S2b Sib28 zeros(8,8) zeros(8,8) zeros(8,8) 
            -S1b' -S1b' -Wb zeros(8,2) zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8)
            Bw' zeros(8,2)' zeros(8,2)' -Gama^2*eye(2,2) zeros(2,8) zeros(2,8) zeros(2,8) hb*Bw' zeros(2,8) zeros(2,8) zeros(2,8)
            zeros(8,8) zeros(8,8) zeros(8,8) zeros(2,8)' -hb*Zb zeros(8,8) zeros(8,8) zeros(8,8) Zb zeros(8,8) zeros(8,8) 
            (hb*N1b)' (hb*N2b)' zeros(8,8) zeros(2,8)' zeros(8,8) -Zb zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8)
            (hb*S1b)' (hb*S2b)' zeros(8,8) zeros(2,8)' zeros(8,8) zeros(8,8) -Zb zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8) 
            Sib18' Sib28' zeros(8,8) (hb*Bw')' zeros(8,8) zeros(8,8) zeros(8,8) Sib88 zeros(8,8) zeros(8,8) Sib811
            zeros(8,8) zeros(8,8) zeros(8,8) zeros(2,8)' Zb' zeros(8,8) zeros(8,8) zeros(8,8) -R zeros(8,8) zeros(8,8)
            (X*C'+Lb'*D')' zeros(8,8) zeros(8,8) zeros(2,8)' zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8) zeros(8,8) -eye(8,8) zeros(8,8)
            (Sib111)' zeros(8,8) zeros(8,8) zeros(2,8)' zeros(8,8) zeros(8,8) zeros(8,8) (Sib811)' zeros(8,8) zeros(8,8) -eps*eye(8,8)];
        
        
        LMI1 = [Rb Xb
            Xb Tb];
        LMI2 = [Rb eye(8,8)
            eye(8,8) R];
        LMI3 = [Xb eye(8,8)
            eye(8,8) X];
        LMI4 = [Tb eye(8,8)
            eye(8,8) T];
        
        0>=Sib;
        LMI1>=0;
        LMI2>=0;
        LMI3>=0;
        LMI4>=0;
               
        
cvx_end


Xbopt = Xb;
Rbopt = Rb;
Tbopt = Tb;
Xopt = X;
Ropt = R;
Topt = T;

index = Xopt*inv(Topt)*Xopt-Ropt;
[rrr,ppp] = chol(index+index');


if strcmp('Solved', cvx_status)==1;
    
    if 0 == ppp
        
        display('feasible solution found !')
        break
    
    end
    
else 
    continue
end


    
end
k
i
if k == kmax;
    break
end
    
end


K_controller = Lb*inv(X);
% K_controller = 10^11*[1.0982 -0.4195 -0.1287 0.0215 -0.1086 -0.1091 -0.1308 -0.1413];
 

uncontrolled_sys = ss(A,Bw,C,Dw);
controlled_sys = ss(A+Bh*K_controller,Bw,C-D*K_controller,Dw);

figure
bodemag(uncontrolled_sys);
hold on
bodemag(controlled_sys);
% title('Bode magnitude comparison between the open loop and close loop system','fontsize',20);
% legend('Uncontrolled','Controlled(4 Digits)','fontsize',14);
% xlabel('Frequency (rad/s)','fontsize',20,'fontweight','b')
% ylabel('Magnitude (db)','fontsize',20,'fontweight','b')

Gama
pb
hb

load elcentro_UP;
Disturbance = elcentro_UP;

Input = Disturbance(:,2);
Time = Disturbance(:,1);

[Output1 time1 states1] = lsim(uncontrolled_sys,[Input,Input],Time);
[Output2 time2 states2] = lsim(controlled_sys,[Input,Input],Time);

figure
plot(time1,Output1(:,1),'--r','LineWidth',1.5)
hold on
plot(time2,Output2(:,1),'-k','LineWidth',1.5)
title('Comparison of Uncontrolled and Controlled(state feedback) first floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,3),'--r','LineWidth',1.5)
hold on
plot(time2,Output2(:,3),'-k','LineWidth',1.5)
title('Comparison of Uncontrolled and Controlled(state feedback) second floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time2,Output2(:,5),'--r','LineWidth',1.5)
hold on
plot(time1,Output1(:,5),'-k','LineWidth',1.5)
title('Comparison of Uncontrolled and Controlled(state feedback) 3rd floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,7),'--r','LineWidth',1.5)
hold on
plot(time2,Output2(:,7),'-k','LineWidth',1.5)
title('Comparison of Uncontrolled and Controlled(state feedback) 4th floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

% acceleration

figure
plot(time1,Output1(:,2),'--r','LineWidth',2)
hold on
plot(time2,Output2(:,2),'-k','LineWidth',2)
title('Comparison of Uncontrolled and Controlled(state feedback) first floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,4),'--r','LineWidth',2)
hold on
plot(time2,Output2(:,4),'-k','LineWidth',2)
title('Comparison of Uncontrolled and Controlled(state feedback) second floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,6),'-k','LineWidth',2)
hold on
plot(time2,Output2(:,6),'--r','LineWidth',2)
title('Comparison of Uncontrolled and Controlled(state feedback) 2rd floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,8),'--r','LineWidth',2)
hold on
plot(time2,Output2(:,8),'-k','LineWidth',2)
title('Comparison of Uncontrolled and Controlled(state feedback) 4th floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

Control_input = K_controller*states2';

figure
plot(time2,Control_input,'-k','LineWidth',2)
title('control effort','fontsize',20);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of Control input','fontsize',20,'fontweight','b')

% uncertainty$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% uncertainty$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% uncertainty$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% uncertainty$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% uncertainty$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% uncertainty$$$$$$$$$$$$$$$$$$$$$$$$$$
% uncertainty$$$$$$$$$$$$

pb = 0.04;

f1 = ureal('f1',pb,'Percentage',200);
f2 = ureal('f2',pb,'Percentage',200);
f3 = ureal('f3',pb,'Percentage',200);
f4 = ureal('f4',pb,'Percentage',200);

f1a = 1/max(f1.Range)*pb;
f2a = 1/max(f2.Range)*pb;
f3a = 1/max(f3.Range)*pb;
f4a = 1/max(f4.Range)*pb;

f11 = f1*f1a;
f21 = f2*f2a;
f31 = f3*f1a;
f41 = f4*f2a;

k1 = 18.050*10^3;
k2 = 340000*10^3;
k3 = 326000*10^3;
k4 = 280000*10^3;

k1 = k1*(1+f11);
k2 = k2*(1+f21);
k3 = k3*(1+f31);
k4 = k4*(1+f41);

DeltaA = pb*[zeros(4,8);-(k1+k2)*f1/m1 k2*f2/m1 0 0 -(c1+c2)*f1/m1 c2*f2/m1 0 0;
    k2*f1/m2 -(k2+k3)*f2/m2 k3*f3/m2 0 c2*f1/m2 -(c2+c3)*f2/m2 c3*f3/m2 0;
    0 k3*f2/m3 -(k3+k4)*f3/m3 k4*f4/m4 0 c3*f2/m3 -(c3+c4)*f3/m3 c4*f4/m4;
    0 0 k4*f3/m4 -k4*f4/m4 0 0 c4*f3/m4 -c4*f4/m4];
    



uncontrolled_sys = ss(A+DeltaA,Bw,C,Dw);
controlled_sys = ss(A+DeltaA+Bh*K_controller,Bw,C-D*K_controller,Dw);

figure
bodemag(uncontrolled_sys);
hold on
bodemag(controlled_sys);


[Output1 time1 states1] = lsim(uncontrolled_sys,[Input,Input],Time);
[Output2 time2 states2] = lsim(controlled_sys,[Input,Input],Time);

figure
plot(time1,Output1(:,1),'--r','LineWidth',1.5)
hold on
plot(time2,Output2(:,1),'-k','LineWidth',1.5)
title('Comparison of Uncontrolled and Controlled(state feedback) first floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,3),'--r','LineWidth',1.5)
hold on
plot(time2,Output2(:,3),'-k','LineWidth',1.5)
title('Comparison of Uncontrolled and Controlled(state feedback) second floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time2,Output2(:,5),'--r','LineWidth',1.5)
hold on
plot(time1,Output1(:,5),'-k','LineWidth',1.5)
title('Comparison of Uncontrolled and Controlled(state feedback) 3rd floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,7),'--r','LineWidth',1.5)
hold on
plot(time2,Output2(:,7),'-k','LineWidth',1.5)
title('Comparison of Uncontrolled and Controlled(state feedback) 4th floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

% acceleration

figure
plot(time1,Output1(:,2),'--r','LineWidth',2)
hold on
plot(time2,Output2(:,2),'-k','LineWidth',2)
title('Comparison of Uncontrolled and Controlled(state feedback) first floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,4),'--r','LineWidth',2)
hold on
plot(time2,Output2(:,4),'-k','LineWidth',2)
title('Comparison of Uncontrolled and Controlled(state feedback) second floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,6),'-k','LineWidth',2)
hold on
plot(time2,Output2(:,6),'--r','LineWidth',2)
title('Comparison of Uncontrolled and Controlled(state feedback) 2rd floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

figure
plot(time1,Output1(:,8),'--r','LineWidth',2)
hold on
plot(time2,Output2(:,8),'-k','LineWidth',2)
title('Comparison of Uncontrolled and Controlled(state feedback) 4th floor','fontsize',20);
legend('Uncontrolled','Controlled','fontsize',14);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of vibration (m)','fontsize',20,'fontweight','b')

Control_input = K_controller*states2';

figure
plot(time2,Control_input,'-k','LineWidth',2)
title('control effort','fontsize',20);
xlabel('Time(s)','fontsize',20,'fontweight','b')
ylabel('Amplitude of Control input','fontsize',20,'fontweight','b')
% Robuststab
[stabmarg,destabunc,report,info] = robuststab(controlled_sys)

% Mussv

[M,Delta,BlkStruct] = lftdata(controlled_sys);
szDelta = size(Delta);
M11 = M(1:szDelta(2),1:szDelta(1));
omega=1:1:1000;
M11_g = frd(M11,omega);
mubnds = mussv(M11_g,BlkStruct,'s');
semilogx(mubnds)

% Wcnorm

[maxgain,maxgainunc,info] = wcgain(controlled_sys);
maxgain
bodemag(controlled_sys.Nominal,'b-',usubs(controlled_sys,maxgainunc),'r--',{.01,1000})
legend('Nominal','Worst case')
