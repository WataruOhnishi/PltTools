function pfig = pltPzmap(sys)
%pltPzmap - plot pole zero map
% pfig = pltPzmap(sys)
% Author    : Wataru Ohnishi, University of Tokyo, 2014
%%%%%


h = figure;
pzmap(sys);
xRange = h.CurrentAxes.XLim;
yRange = h.CurrentAxes.YLim;
close(h);

pole_sys = pole(sys);
zero_sys = zero(sys);

if sys.Ts == 0 % continuous
    hfig = figure;
    xlabel('Real axis');
    ylabel('Imaginary axis');
    hold on;
    plot([max(abs(xRange)),-max(abs(xRange))]*10,[0 0],'k-','linewidth',0.5);
    plot([0 0],[max(abs(yRange)),-max(abs(yRange))]*10,'k-','linewidth',0.5);
    plot(real(pole_sys) ,imag(pole_sys),'kx');
    plot(real(zero_sys), imag(zero_sys),'ko');
    box on;
    xlim(xRange); ylim(yRange);
    hold off
    
    pfig = pubfig(hfig);
    pfig.Grid = 'off';
    pfig.MinorTick = 'off';
    pfig.Dimension = [12,10];
    pfig.LineWidth = [1.5,1.5,1,1];
    pfig.MarkerSize = [10 10 6 6];
    
else % discrete
    uc_point = 1e3; % unitcircleのplot点数
    uc_x = zeros(1,uc_point);
    uc_y = zeros(1,uc_point);
    for k = 1:1:uc_point
        uc_x(k) = cos(2*(pi+2*2*pi/uc_point)*k/uc_point+pi/2);
        uc_y(k) = sin(2*(pi+2*2*pi/uc_point)*k/uc_point+pi/2);
    end
    
    hfig = figure;
    xlabel('Real axis');
    ylabel('Imaginary axis');
    hold on;
    plot(uc_x,uc_y,'k--','linewidth',1);
    plot([max(abs(xRange)),-max(abs(xRange))]*10,[0 0],'k-','linewidth',0.5);
    plot([0 0],[max(abs(yRange)),-max(abs(yRange))]*10,'k-','linewidth',0.5);
    plot(real(pole_sys) ,imag(pole_sys),'kx');
    plot(real(zero_sys), imag(zero_sys),'ko');
    box on;
    xlim(xRange); ylim(yRange);
    hold off
    pfig = pubfig(hfig);
    pfig.Grid = 'off';
    pfig.MinorTick = 'off';
    pfig.Dimension = [12,10];
    pfig.LineWidth = [1.5,1.5,1,1,1];
    pfig.MarkerSize = [10 10 6 6 6];
    
end

