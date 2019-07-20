function [pfig, hfig] = pltNyquist(in,option)
%pltNyquist - plot Nyquist diagram
%
% pfig = pltNyquist(in,option)
%   in              : {data1, data2, ..., data5}
% Required
%   data1.sys       : tf, ss, zpk, frd
% Optional
%   data1.color     : line color, ('r','g','b','k','m','c','g2','b2','b3') or RGB space
%   data1.name = 'test';% legend name
%   data1.style = '-';  % '-', '--', '-.', '.'
%   data1.marker = '*'; % '*', '.'
%   option.xmin = -1.5; % min x axis
%   option.xmax = 1;    % max x axis
%   option.ymin = -1.5; % min y axis
%   option.ymax = 1;    % max y axis
%   option.fmin = 1;    % freq min
%   option.fmax = 1e3;  % freq max
%   option.title = '';  % title
%   option.gmdb = 6;    % gain margin for circle condition
%   option.pmdeg = 30;  % phase margin for circle condition
%   option.Smax_dB = 6;  % max of sensitivity function

% Advanced options
%   data1.P         : plant model
%   data1.Cfb       : FB controller model
%   data1.Csh       : shaping filters
%   option.plot     : 'sys' : plot data{k}.sys
%                     'Gop' : open loop (PCfbCsh)
% Author    : Wataru Ohnishi, University of Tokyo, 2017


N = length(in);

if nargin < 2
    option = struct;
end

data = cell(1,N); % to accept pltNyquist(tf(1)) or pltNyquist({tf(1)})
for k = 1:1:N
    if ~iscell(in), in = {in}; end
    try data{k}.sys = in{k}.sys; data{k} = in{k}; catch, data{k}.sys = in{k}; end
    
    % fdel
    if isfield(option,'fmin')
        [~,kmin] = min(data{k}.sys.freq - option.fmin);
        freqdel = data{k}.sys.freq(1:kmin);
        fprintf('freqs from %.1f Hz to %.1fHz deleted\n',data{k}.sys.freq(1),data{k}.sys.freq(kmin));
        data{k}.sys = fdel(data{k}.sys,freqdel);
    end
    if isfield(option,'fmax')
        kmax = knnsearch(data{k}.sys.freq,option.fmax);
        freqdel = data{k}.sys.freq(kmax:end);
        fprintf('freqs from %.1f Hz to %.1fHz deleted\n',data{k}.sys.freq(kmax),data{k}.sys.freq(end));
        data{k}.sys = fdel(data{k}.sys,freqdel);
    end
end


try option.plot; catch, option.plot = 'sys'; end
for k = 1:1:N
    if strcmp(option.plot,'Gop'), data{k}.sys = data{k}.P*data{k}.Cfb*data{k}.Csh; end
end

try option.xmin; catch, option.xmin = -1.5; option.xtick = 0.5; option.xmax = 1; end
try option.ymin; catch, option.ymin = -1.5; option.ytick = 0.5; option.ymax = 1; end

try option.fmin; catch, option.fmin = 0.1; option.fmax = 1e4; end
freq = logspace(log10(option.fmin),log10(option.fmax),3000);

uc_point = 1e3;
uc_x = zeros(1,uc_point);
uc_y = zeros(1,uc_point);
for k = 1:1:uc_point
    uc_x(k) = cos(2*(pi+2*2*pi/uc_point)*k/uc_point+pi/2);
    uc_y(k) = sin(2*(pi+2*2*pi/uc_point)*k/uc_point+pi/2);
end

colorlist = {'b','r','k','m','g','c','g2','b2','b3'};
for k = 1:1:N
    try data{k}.color; catch, data{k}.color = colorlist{mod(k,9)}; end
    try data{k}.name; catch, data{k}.name = num2str(k); end
end

for k = 1:1:N
    try data{k}.sys.ResponseData;
    catch
        % if data is not frd, disp isstable
        Gcl = feedback(data{k}.sys,1);
        try
            if isstable(Gcl) == 0
                fprintf('data{%d} is unstable\n',k);
            end
        catch
        end
        data{k}.sys = frd(data{k}.sys,freq,'Hz');
    end
    
    try data{k}.style; catch, data{k}.style = '-'; end
    
    % color options
    data{k}.color = str2rgb(data{k}.color);
    
    % nyquist data
    [data{k}.re,data{k}.im,data{k}.w] = nyquist(data{k}.sys);
end


hfig = figure;
for k = 1:1:N
    h = plot(squeeze(data{k}.re),squeeze(data{k}.im)); hold on;
    set(h,'Color',data{k}.color);
    try set(h,'linestyle',data{k}.style); catch, end
    try set(h,'Marker',data{k}.marker); catch, end
end

xlabel('Real axis'); ylabel('Imaginary axis');
set(gca,'xtick',option.xmin:option.xtick:option.xmax);
set(gca,'ytick',option.ymin:option.ytick:option.ymax);
xlim([option.xmin,option.xmax]);
ylim([option.ymin,option.ymax]);

plot(uc_x,uc_y,'k--');
plot([0,0],[min(-2,option.ymin),max(2,option.ymax)],'k-'); plot([min(-2,option.xmin),max(2,option.xmax)],[0,0],'k-'); % real and imag axes
hm = plot(-1,0,'kx');
set(hm,'MarkerSize',12);
try title(option.title); catch, end

if isfield(option,'gmdb')
    [sigma, rm] = stabCircle(option.gmdb,option.pmdeg);
    [xstab,ystab] = circle(-sigma,0,rm);
    plot(xstab,ystab,'k:');
end
if isfield(option,'Smax_dB')
    [xSmax,ySmax] = circle(-1,0,1/(db2mag(option.Smax_dB)));
    plot(xSmax,ySmax,'k--');
end



multiLegend(data);

pfig = pubfig(hfig);pfig.Grid = 'off';
pfig.LegendLoc = 'southeast';
