clc
clear all
close all
respiratoryFinal     % model is your model name  
print('-sname','-dbitmap','new_name')
im = imread('new_name.bmp')
% You can then save it as a jpg file
imwrite(im,'respiratotyFinal.jpg')