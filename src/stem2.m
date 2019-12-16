function stem2(x,y,linespec)
% stem plot for pubfig
linespec = ['''',linespec,''''];
for k = 1:length(x)
    plot([x(k),x(k)],[0,y(k)],eval(linespec));
end
