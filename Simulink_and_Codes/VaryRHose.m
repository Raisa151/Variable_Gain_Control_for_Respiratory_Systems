clc
clear all
close all

%%
R_lung = 5/1000;
C_lung = 20;
R_leak = 60/1000;
R_hose = [4.5/500 4.5/1000 4.5/1500];
wn = 2*pi*30;

for i = 1:length(R_hose)
    A_h(i) = -(1/R_hose(i) + 1/R_leak)/(R_lung*C_lung*(1/R_lung+1/R_hose(i)+1/R_leak));
    B_h(i) = (1/R_hose(i))/(R_lung*C_lung*(1/R_lung+1/R_hose(i)+1/R_leak));
    C_h(:,i) = [( 1/R_lung/(1/R_lung+1/R_hose(i)+1/R_leak)),  -(1/R_hose(i)+1/R_leak)/(R_lung*(1/R_lung+1/R_hose(i)+1/R_leak))]';
    D_h(:,i) = [ (1/R_hose(i)/(1/R_lung+1/R_hose(i)+1/R_leak)) , (1/R_hose(i) / (R_lung*(1/R_lung+1/R_hose(i)+1/R_leak)))]';
end
%%
s = tf("s");
%Blower B_s
zeta = 1;
B_s = (wn^2)/(s^2 + 2*zeta*wn*s + wn^2);

fprintf("Blower transfer function is: \n")
B_s

for i = 1:3
    H_s(: ,i) = C_h(:,i)*(inv(s*eye(1)-A_h(i)))*B_h(i) + D_h(:,i);
    fprintf("Patient hose transfer function for R_hose %f is: \n",R_hose(i))
    H_s(:,i)
    
    fprintf("Overall plant transfer function for R_hose %f is  is : \n",R_hose(i))
    B_s*H_s(:,i)
end

%%
t = 0:0.001:16;
pset_t = 5*(t>=0)+(75*((t-1).*(t-1>=0)))-(75*((t-1.2).*(t-1.2>=0)))-15*(t-5>=0);%slope of the ramp 75mbar/sec
pset_t_secondcycle =(75*(((t-8)-1).*((t-8)-1>=0)))-(75*(((t-8)-1.2).*((t-8)-1.2>=0)))-15*((t-8)-5>=0);
u = pset_t+pset_t_secondcycle;

figure(1)
grid on
plot(t,[pset_t+pset_t_secondcycle])
xlabel("Time[s]")
ylabel("Paw[mbar]")
title("Effect of Varying RHose on Paw")
for i = 1:length(R_hose)
    
    C_s = 0.4/s;
    closed_loopsys_Paw(i) = (1+C_s)*((B_s*H_s(1,i))/(1+B_s*H_s(1,i)*C_s));
    fprintf("Closed loop Paw transfer function is : \n")
    closed_loopsys_Paw(i)

    [yPaw,taw,x] = lsim(closed_loopsys_Paw(i),u,t);
    hold on
    plot(taw,yPaw)
    
end
grid on
hold off
 
%%
for i = 1:length(R_hose)
    
    C_s = 0.4/s;
    H_s0(i)= (1+C_s)*((B_s*H_s(1,i))/(1+(B_s*H_s(1,i)*C_s)));
    closed_loopsys_Qpat(i) = (B_s*H_s(2,i))*(1+C_s - (H_s0(i)*C_s));
       
    fprintf("Closed loop Qpat transfer function is : \n")
    closed_loopsys_Qpat(i)

    [ypat,tpat,x] = lsim(closed_loopsys_Qpat(i),u,t);
    figure(2)
    plot(tpat,ypat*60)
    xlabel("Time[s]")
    ylabel("Qpat[mL/min]")
    title("Effect of Varying RHose on Qpat")
    hold on
end
grid on
hold off