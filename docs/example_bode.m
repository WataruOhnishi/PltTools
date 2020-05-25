% Bode plot example
clear; close all;
s = tf('s');
%%
sys1 = 10*2*pi/(s+10*2*pi);

w = 100*2*pi;
sys2 = w^2/(s^2 + 2*0.2*w*s + w^2);

pfig = pltBode({sys1,sys2});
pfig.LegendLoc = 'southwest';

%%
clear data1 data2 option
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

%%
clear data1 data2 option
data1.sys = 10*2*pi/(s+10*2*pi);
data1.color = 'b'; 
data1.name = 'system1';

w = 100*2*pi;
data2.sys = frd(w^2/(s^2 + 2*0.2*w*s + w^2),logspace(0,2,20),'FrequencyUnit','Hz');
data2.color = 'r'; 
data2.style = ''; 
data2.name = 'system2';

option.fmin = 1; option.fmax = 100; option.pmin = -180; option.pmax = 180; 
pfig = pltBode({data1,data2},option);
pfig.LegendLoc = 'southwest';
expfig('plot/bode2','-pdf','-png');
