% pole-zero plot example
clear; close all;
s = tf('s');

w = 100*2*pi;
sys = 10*2*pi/(s+10*2*pi)/s * (s^2 + 2*0.7*w*s + w^2)/(s^2 + 2*0.2*w*s + w^2)*(s-100)/(s+100);
pfig = pltPzmap(sys);
expfig('plot/pzmap1','-pdf','-png');
pfig = pltPzmap(c2d(sys,100e-6,'zoh'));
expfig('plot/pzmap2','-pdf','-png');
xlim([0.9,1.1]); ylim([-0.1,0.1]);
expfig('plot/pzmap2_zoom','-pdf','-png');
