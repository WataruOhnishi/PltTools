% Bode plot example
clear; close all;
s = tf('s');

data1.sys = 10*2*pi/(s+10*2*pi);
data1.color = 'b'; 
data1.name = 'system1';

w = 100*2*pi;
data2.sys = w^2/(s^2 + 2*0.2*w*s + w^2);
data2.color = 'r'; 
data2.name = 'system2';

pfig = pltBode({data1,data2});
pfig.LegendLoc = 'southwest';
expfig('plot/bode1','-pdf','-png');
