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
%   option.Smax = tf(2);  % max of sensitivity function
%   option.multiFRDcolor = true;
%
% Author    : Wataru Ohnishi, University of Tokyo, 2017


if nargin < 2
    option = struct;
end

data = iscellinput(in);
option = setoptions(data,option);
freq = option.freq;
N = length(data);
data = shapeData(data,freq,option);


hfig = figure;
pltUc_Axis(option)
if isfield(option,'title'), title(option.title); end
pltGmPm(option);
pltSmax(data,option);
p = plot(-1,0,'ko','MarkerFaceColor','k');
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

colormap('lines');
main_pltNyquist(data,N,option)
multiLegend(data);

xlabel('Real axis'); ylabel('Imaginary axis');
set(gca,'xtick',option.xmin:option.xtick:option.xmax);
set(gca,'ytick',option.ymin:option.ytick:option.ymax);
xlim([option.xmin,option.xmax]);
ylim([option.ymin,option.ymax]);
daspect([1,1,1]);

pfig = pubfig(hfig); pfig.Grid = 'off';
pfig.LegendLoc = 'southeast';
pfig.MarkerSize = 7;
end


function data = iscellinput(in)
if iscell(in), N = length(in); else, N = 1; end
data = cell(1,N); % to accept pltBode(tf(1)) or pltBode({tf(1)})
for k = 1:1:N
    if iscell(in) == 0
        in = {in};
    end
    if isfield(in{k},'sys')
        data{k}.sys = in{k}.sys;
        data{k} = in{k};
    else
        data{k}.sys = in{k};
    end
end
end


function option = setoptions(data,option)
if ~isfield(option,'xmin'), option.xmin = -1.5; option.xtick = 0.5; option.xmax = 1; end
if ~isfield(option,'ymin'), option.ymin = -1.5; option.ytick = 0.5; option.ymax = 1; end
if ~isfield(option,'multiFRDcolor'), option.multiFRDcolor = false; end    
if ~isfield(option,'fmin')
    isfrd = false;
    freqtemp = [];
    for k = 1:length(data)
        if isa(data{k}.sys,'frd')
            freqtemp = [freqtemp;data{k}.sys.freq];
            isfrd = true;
        end
    end
    if isfrd
        option.fmin = min(freqtemp); option.fmax = max(freqtemp);
    else % all parametric model
        option.fmin = 0.1; option.fmax = 1e4;
    end
end
option.freq = logspace(log10(option.fmin),log10(option.fmax),1000);
end


function data = shapeData(data,freq,option)
colorlist = {'b','r','k','m','g','c','g2','b2','b3'};
for k = 1:1:length(data)
    % convert to FRD
    if isnumeric(data{k}.sys)
        data{k}.sys = tf(data{k}.sys);
    end
    %     if ~isfield(data{k}.sys,'ResponseData')
    if ~isa(data{k}.sys,'frd')
        [mag,phase,w] = bode(data{k}.sys,freq*2*pi);
        data{k}.sys = frd(squeeze(mag).*exp(1j*deg2rad(squeeze(phase))),w/2/pi,'FrequencyUnit','Hz');
    end
    
    if (strcmp(data{k}.sys.FrequencyUnit,'rad/s') == 1)
        data{k}.sys.Frequency = data{k}.sys.Frequency/2/pi;
        data{k}.sys.FrequencyUnit = 'Hz';
    end
    
    if (strcmp(data{k}.sys.FrequencyUnit,'rad/TimeUnit') && strcmp(data{k}.sys.TimeUnit,'seconds'))
        data{k}.sys.Frequency = data{k}.sys.Frequency/2/pi;
        data{k}.sys.FrequencyUnit = 'Hz';
    end
    
    if ~isfield(data{k},'style'), data{k}.style = '-'; end
    if ~isfield(data{k},'color'), data{k}.color = colorlist{mod(k,9)}; end
    if ~isfield(data{k},'name'), data{k}.name = num2str(k); end
    
    % color options
    data{k}.color = str2rgb(data{k}.color);
    
    % pltNyquist
    if ~isa(data{k}.sys,'frd')
        Gcl = feedback(data{k}.sys,1);
        if isstable(Gcl) == 0
            fprintf('data{%d} is unstable\n',k);
        end
    end
    
    % fdel
    if isfield(option,'fmin')
        [~,kmin] = min(abs(data{k}.sys.freq - option.fmin));
        freqdel = data{k}.sys.freq(1:kmin-1);
        if data{k}.sys.freq(1) ~= data{k}.sys.freq(kmin)
            fprintf('freqs from %.1f Hz to %.1fHz deleted\n',data{k}.sys.freq(1),data{k}.sys.freq(kmin));
        end
        data{k}.sys = fdel(data{k}.sys,freqdel);
    end
    if isfield(option,'fmax')
        [~,kmax] = min(abs(data{k}.sys.freq - option.fmax));
        freqdel = data{k}.sys.freq(kmax+1:end);
        if data{k}.sys.freq(kmax) ~= data{k}.sys.freq(end)
            fprintf('freqs from %.1f Hz to %.1fHz deleted\n',data{k}.sys.freq(kmax),data{k}.sys.freq(end));
        end
        data{k}.sys = fdel(data{k}.sys,freqdel);
    end
    
end
end


function pltUc_Axis(option)
[uc_x,uc_y] = unitcircle;
p1 = plot(uc_x,uc_y,'k:'); hold on
p2 = plot([0,0],[min(-2,option.ymin),max(2,option.ymax)],'k-');
p3 = plot([min(-2,option.xmin),max(2,option.xmax)],[0,0],'k-'); % real and imag axes
set(get(get(p1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end


function [uc_x,uc_y] = unitcircle(uc_point)
if nargin < 1, uc_point = 1e3; end
uc_x = zeros(1,uc_point);
uc_y = zeros(1,uc_point);
for k = 1:1:uc_point
    uc_x(k) = cos(2*(pi+2*2*pi/uc_point)*k/uc_point+pi/2);
    uc_y(k) = sin(2*(pi+2*2*pi/uc_point)*k/uc_point+pi/2);
end
end


function main_pltNyquist(data,N,option)
for k = 1:1:N
    if length(data{k}.sys) > 1
        data{k}.re = [];
        data{k}.im = [];
        for kk = 1:length(data{k}.sys)
            [re,im,~] = nyquist(data{k}.sys(:,:,kk));
            
            h = plot(squeeze(re),squeeze(im)); hold on;
            if ~option.multiFRDcolor, set(h,'Color',data{k}.color); end
            try set(h,'linestyle',data{k}.style); catch, end
            try set(h,'Marker',data{k}.marker); catch, end
            if kk > 1, set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); end
        end
    else
        [data{k}.re,data{k}.im,data{k}.w] = nyquist(data{k}.sys);
        h = plot(squeeze(data{k}.re),squeeze(data{k}.im)); hold on;
        set(h,'Color',data{k}.color);
        try set(h,'linestyle',data{k}.style); catch, end
        try set(h,'Marker',data{k}.marker); catch, end
        
    end
end

end


function pltGmPm(option)
if isfield(option,'gmdb')
    [sigma, rm] = stabCircle(option.gmdb,option.pmdeg);
    [xstab,ystab] = circle(-sigma,0,rm);
    p = plot(xstab,ystab,'k:');
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

if isfield(option,'gmpmdot')
    if option.gmpmdot == true
        p = plot(-sigma,0,'ko','MarkerFaceColor','k');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
end


function pltSmax(data,option)
if isfield(option,'Smax')
    % find unique Smax
    if ~isfield(option.Smax,'ResponseData')
        option.Smax = frd(option.Smax,data{1}.sys.freq,'FrequencyUnit',data{1}.sys.FrequencyUnit);
    end
    Smax_gain = squeeze(option.Smax.ResponseData);
    Smax_array = Smax_gain(1);
    idx = find(diff(Smax_gain));
    Smax_array = [Smax_array;Smax_array(idx+1);];
    for k = 1:length(Smax_array)
        [xSmax,ySmax] = circle(-1,0,1/Smax_array(k));
        p = plot(xSmax,ySmax,'k-.');
        set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
end
