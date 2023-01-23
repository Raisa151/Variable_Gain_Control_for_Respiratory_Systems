t = 0:0.001:10;
unitstep = (t-2)>=0;
ramp = (t-2).*unitstep;
%plot(t,unitstep)
% hold on
plot(t,ramp(t-2))