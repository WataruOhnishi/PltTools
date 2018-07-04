% Nyquist plot example
clear; close all;
s = tf('s');

w = 100*2*pi;
data1.sys = w^2/(s^2 + 2*0.1*w*s + w^2);
data1.color = 'b'; 
data1.name = 'system1';

pfig = pltNyquist({data1});
expfig('plot/nyquist1','-pdf','-png');
