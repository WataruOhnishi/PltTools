function [pfig, hfig] = pltBodemag(in,option)
%pltBode - plot Bode diagram
%
% pfig = pltBode(in,option)
%   in              : {data1, data2, ..., data5}
% Required
%   data1.sys       : tf, ss, zpk, frd
% Optional
%   data1.color = 'b'; %line color, ('r','g','b','k','m','c','g2','b2','b3') or RGB space
%   data1.name = 'test'; %legend name
%   data1.style = '-'; %'-', '--', '-.', '.'
%   data1.marker = '*'; %'*', '.'
%   option.fmin = 1; % min freq
%   option.fmax = 1000; % max freq;
%   option.gmin = -40; % gain min
%   option.gtick = 20; % gain tick
%   option.gmax = 40; % gain min
%   option.foption = 'lin'; % freq plot options 'log' or 'lin'
%   option.Legpos   : legend in gain or phase plot 'gain' or 'phase'
%   option.title = 'Bode plot'; % title
%   cohflag         : coherence plot
% Author    : Wataru Ohnishi, University of Tokyo, 2020
%%%%%


N = length(in);

if nargin < 2
    option = struct;
end

data = cell(1,N); % to accept pltBode(tf(1)) or pltBode({tf(1)})
for k = 1:1:N
    if iscell(in) == 0
        in = {in};
    end
    try
        data{k}.sys = in{k}.sys;
        data{k} = in{k};
    catch
        data{k}.sys = in{k};
    end
end

option.plot = 'sys';

if ~isfield(option,'fmin'), option.fmin = 1; option.fmax = 1000; end
if ~isfield(option,'pmin'), option.pmin = -360; option.pmax = 0; end
if ~isfield(option,'foption'), option.foption = 'log'; end
if ~isfield(option,'Legpos'), option.Legpos = 'gain'; end
if ~isfield(option,'LegendLoc'), option.LegendLoc = 'best'; end
if ~isfield(option,'freq'), freq = logspace(log10(option.fmin),log10(option.fmax),1000); else, freq = option.freq; end

if isfield(option,'m') && isfield(option,'wc')
    s = tf('s');
    data_con.sys = (s/option.wc)^option.m;
    data_con.name = 'constraint';
    data_con.color = 'k';
    data_con.style = '--';
    data = [data,data_con];
end
if isfield(option,'Smax')
    data_con2.sys = option.Smax;
    data_con2.name = 'constraint';
    data_con2.color = 'k';
    data_con2.style = ':';
    data = [data,data_con2];
end
N = length(data);

colorlist = {'b','r','k','m','g','c','g2','b2','b3'};
for k = 1:1:N
    % convert to FRD
    if isnumeric(data{k}.sys)
        data{k}.sys = tf(data{k}.sys);
    end
    if ~isfield(data{k}.sys,'ResponseData')
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
    
    try data{k}.style; catch, data{k}.style = '-'; end
    try data{k}.color; catch, data{k}.color = colorlist{mod(k,9)}; end
    try data{k}.name; catch, data{k}.name = num2str(k); end
    
    % color options
    data{k}.color = str2rgb(data{k}.color);
end

if ~isfield(option,'gmin')
    h = figure;
    for k = 1:N
        bode(data{k}.sys); hold on
    end
    option.gmin = h.Children(3).YLim(1);
    option.gmax = h.Children(3).YLim(end);
    option.gtick = h.Children(3).YTick(2)-h.Children(3).YTick(1);
    close(h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hfig = figure;
for k = 1:1:N
    if strcmp(option.foption, 'log')
        if strcmp(data{k}.style,'')
            h = scatter(data{k}.sys.frequency,mag2db(abs(squeeze(data{k}.sys.ResponseData))),'filled','o'); hold on
            set(h,'MarkerEdgeColor',data{k}.color);
            set(h,'MarkerFaceColor',data{k}.color);
            set(gca,'xscale','log');
        else
            h = semilogx(data{k}.sys.frequency,mag2db(abs(squeeze(data{k}.sys.ResponseData)))); hold on;
            try set(h,'Color',data{k}.color); catch, end
            try set(h,'linestyle',data{k}.style); catch, end
            try set(h,'Marker',data{k}.marker); catch, end
        end
    else
        h = plot(data{k}.sys.frequency,mag2db(abs(squeeze(data{k}.sys.ResponseData)))); hold on;
        try set(h,'Color',data{k}.color); catch, end
        try set(h,'linestyle',data{k}.style); catch, end
        try set(h,'Marker',data{k}.marker); catch, end
    end
end
axis([option.fmin, option.fmax, option.gmin, option.gmax]);
set(gca,'ytick',option.gmin:option.gtick:option.gmax);
ylabel('Magnitude [dB]');
if strcmp(option.Legpos,'gain'), multiLegend(data); end
grid on; hold on;
try title(option.title); catch, end
xlabel('Frequency [Hz]');

pfig = pubfig(hfig);
pfig.LegendLoc = option.LegendLoc;
