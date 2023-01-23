clc
clear all
close all

s = tf('s')
k = -4.985
g = (s+40)/ (s+10)

h = 1/(s+20)
r = 1/s
forward= (c+r)*g
close = feedback(forward,h)

stepinfo(G)
step(G,0:0.003:10)
