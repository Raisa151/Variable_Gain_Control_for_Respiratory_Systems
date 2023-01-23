clc
clear all
close all

%%
R_lung = 5/1000;
C_lung = 20;
R_leak = 60/1000;
R_hose = 4.5/1000;
wn = 2*pi*30;

A_h = -(1/R_hose + 1/R_leak)/(R_lung*C_lung*(1/R_lung+1/R_hose+1/R_leak));
B_h = (1/R_hose)/(R_lung*C_lung*(1/R_lung+1/R_hose+1/R_leak));
C_h = [( 1/R_lung/(1/R_lung+1/R_hose+1/R_leak)),  -(1/R_hose+1/R_leak)/(R_lung*(1/R_lung+1/R_hose+1/R_leak))]';
D_h = [ (1/R_hose/(1/R_lung+1/R_hose+1/R_leak)) , (1/R_hose / (R_lung*(1/R_lung+1/R_hose+1/R_leak)))]';

fprintf("Value of A_h B_h C_h & D_h  is : %f\n %f\n ", A_h,B_h);
fprintf("Value of C_h is: \n")
disp(C_h)
fprintf("Value of D_h is: \n")
disp(D_h)

%%
%patient hose H_s
s = tf("s");
H_s = C_h*(inv(s*eye(1)-A_h))*B_h + D_h;
fprintf("Patient hose transfer function is: \n")
H_s

%Blower B_s
zeta = 1;
B_s = (wn^2)/(s^2 + 2*zeta*wn*s + wn^2);

fprintf("Blower transfer function is: \n")
B_s

fprintf("Overall plant transfer function is : \n")
B_s*H_s


t = 0:0.001:16;
pset_t = 5*(t>=0)+(75*((t-1).*(t-1>=0)))-(75*((t-1.2).*(t-1.2>=0)))-15*(t-5>=0);%slope of the ramp 75mbar/sec
pset_t_secondcycle =(75*(((t-8)-1).*((t-8)-1>=0)))-(75*(((t-8)-1.2).*((t-8)-1.2>=0)))-15*((t-8)-5>=0);
u = pset_t+pset_t_secondcycle;
figure(1)
plot(t,[pset_t+pset_t_secondcycle])
xlabel("Time[s]")
ylabel("Paw[mbar]")
grid on

ki =[0,0.4,10];

%Paw plot for different ki
for i = 1:length(ki)
    
    C_s(i) = ki(i)/s;
    closed_loopsys_Paw(i) = (1+C_s(i))*((B_s*H_s(1))/(1+B_s*H_s(1)*C_s(i)));
    fprintf("Closed loop Paw transfer function is : \n")
    closed_loopsys_Paw(i)

    [yPaw,taw,x] = lsim(closed_loopsys_Paw(i),u,t);
    hold on
    plot(taw,yPaw)
end
grid on
hold off

%%
t = 0:0.001:16;
pset_t = 5*(t>=0)+(75*((t-1).*(t-1>=0)))-(75*((t-1.2).*(t-1.2>=0)))-15*(t-5>=0);%slope of the ramp 75mbar/sec
pset_t_secondcycle =(75*(((t-8)-1).*((t-8)-1>=0)))-(75*(((t-8)-1.2).*((t-8)-1.2>=0)))-15*((t-8)-5>=0);
u = pset_t+pset_t_secondcycle;

ki =[0,0.4,10];
for i = 1:length(ki)
    
    C_s(i) = ki(i)/s;
    H_s0(i)= (1+C_s(i))*((B_s*H_s(1))/(1+(B_s*H_s(1)*C_s(i))));
    closed_loopsys_Qpat(i) = (B_s*H_s(2))*(1+C_s(i) - (H_s0(i)*C_s(i)));
       
    fprintf("Closed loop Qpat transfer function is : \n")
    closed_loopsys_Qpat(i)

    [ypat,tpat,x] = lsim(closed_loopsys_Qpat(i),u,t);
    figure(2)
    plot(tpat,ypat*60)
    xlabel("Time[s]")
    ylabel("Qpat[mL/min]")
    hold on
end
grid on
hold off

%%
figure(3)
bode(B_s*H_s(1))
%%
% root_locus
integ_controller = 1/s;
rlocus(B_s*H_s(1)*integ_controller);
sgrid(zeta,wn)
axis([-250 20 -400 400])