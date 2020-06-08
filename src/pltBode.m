function [pfig, hfig, ax] = pltBode(in,option)
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
%   option.pmin = -180; % phase min
%   option.ptick = 90; % phase tick
%   option.pmax = 180; % phase max
%   option.gmin = -40; % gain min
%   option.gtick = 20; % gain tick
%   option.gmax = 40; % gain min
%   option.foption = 'lin'; % freq plot options 'log' or 'lin'
%   option.Legpos   : legend in gain or phase plot 'gain' or 'phase'
%   option.title = 'Bode plot'; % title
%   option.noPhae = false;
%   option.multiFRDcolor = true; % change color for multi frd
% Author    : Wataru Ohnishi, University of Tokyo, 2017
%%%%%

if nargin < 2
    option = struct;
end

data = iscellinput(in);
option = setoptions(data,option);
freq = option.freq;
[option,data] = setconstraints(data,option);

N = length(data);
data = shapeData(data,freq);
option = findGainRange(data,option);

% plot figure
ax = nan;
hfig = figure;
if ~option.noPhase
    ax(1) = subplot(2,1,1);
    pltGain(data,option);
    
    ax(2) = subplot(2,1,2);
    pltPhase(data,option);
    linkaxes(ax,'x');
else
    pltGain(data,option);
end


if exist('pubfig','file')
    pfig = pubfig(hfig);
    pfig.LegendLoc = option.LegendLoc;
else
    pfig = [];
end
if option.noPhase, ax = []; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if ~isfield(option,'fmin'), option.fmin = 1; option.fmax = 1000; end
if ~isfield(option,'pmin'), option.pmin = -360; option.pmax = 0; end
if ~isfield(option,'ptick'), option.ptick = 90; end
if ~isfield(option,'foption'), option.foption = 'log'; end
if ~isfield(option,'Legpos'), option.Legpos = 'gain'; end
if ~isfield(option,'LegendLoc'), option.LegendLoc = 'best'; end
if ~isfield(option,'phasePlot'), option.phasePlot = 1:length(data); end
if ~isfield(option,'freq'), option.freq = logspace(log10(option.fmin),log10(option.fmax),1000); end
if ~isfield(option,'noPhase'), option.noPhase = false; end
if ~isfield(option,'multiFRDcolor'), option.multiFRDcolor = false; end
end

function [option,data] = setconstraints(data,option)
if isfield(option,'m')
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
end

function data = shapeData(data,freq)
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
end
end

function option = findGainRange(data,option)
if ~isfield(option,'gmin')
    h = figure;
    for k = 1:length(data)
        bode(data{k}.sys); hold on
    end
    xlim([option.fmin,option.fmax]);
    option.gmin = h.Children(3).YLim(1);
    option.gmax = h.Children(3).YLim(end);
    if abs(option.gmin) > 1, option.gmin = round(option.gmin); end
    if abs(option.gmax) > 1, option.gmax = round(option.gmax); end
    option.gtick = h.Children(3).YTick(2)-h.Children(3).YTick(1);
    close(h);
end
end

function pltGain(data,option)
gain = @(data) mag2db(abs(squeeze(data.ResponseData)));
for k = 1:1:length(data)
    for kk = 1:length(data{k}.sys)
        if strcmp(option.foption, 'log')
            if strcmp(data{k}.style,'')
                h = scatter(data{k}.sys.frequency,gain(data{k}.sys(:,:,kk)),'filled','o'); hold on
                if ~option.multiFRDcolor
                    set(h,'MarkerEdgeColor',data{k}.color);
                    set(h,'MarkerFaceColor',data{k}.color);
                end
                set(gca,'xscale','log');
            else
                h = semilogx(data{k}.sys.frequency,gain(data{k}.sys(:,:,kk))); hold on;
                if ~option.multiFRDcolor, try set(h,'Color',data{k}.color); catch, end; end
                try set(h,'linestyle',data{k}.style); catch, end
                try set(h,'Marker',data{k}.marker); catch, end
            end
        else
            h = plot(data{k}.sys.frequency,gain(data{k}.sys(:,:,kk))); hold on;
        end
        if kk > 1, set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); end
    end
end
axis([option.fmin, option.fmax, option.gmin, option.gmax]);
set(gca,'ytick',option.gmin:option.gtick:option.gmax);
ylabel('Magnitude [dB]');
if strcmp(option.Legpos,'gain'), multiLegend(data); end
grid on; hold on;
if isfield(option,'title'), title(option.title); end
end

function pltPhase(data,option)
for k = option.phasePlot
    for kk = 1:length(data{k}.sys)
        phasedeg = rad2deg(angle(squeeze(data{k}.sys(:,:,kk).ResponseData)));
        for kk = 1:1:length(phasedeg)
            while phasedeg(kk) > option.pmax
                phasedeg(kk) = phasedeg(kk) - 360;
            end
            while phasedeg(kk) < option.pmin
                phasedeg(kk) = phasedeg(kk) + 360;
            end
        end
        
        if strcmp(option.foption, 'log')
            if strcmp(data{k}.style,'')
                h = scatter(data{k}.sys.frequency,phasedeg,'filled','o'); hold on
                if ~option.multiFRDcolor
                    set(h,'MarkerEdgeColor',data{k}.color);
                    set(h,'MarkerFaceColor',data{k}.color);
                end
                set(gca,'xscale','log');
            else
                h = semilogx(data{k}.sys.frequency,phasedeg); hold on;
                if ~option.multiFRDcolor, try set(h,'Color',data{k}.color); catch, end; end
                set(h,'linestyle',data{k}.style);
                if isfield(data{k},'marker'), set(h,'Marker',data{k}.marker); end
            end
        else
            h = plot(data{k}.sys.frequency,phasedeg); hold on;
        end
        if kk > 1, set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); end
    end
end
axis([option.fmin, option.fmax, option.pmin, option.pmax]);
set(gca,'ytick',option.pmin:option.ptick:option.pmax);
xlabel('Frequency [Hz]');
ylabel('Phase [deg]');
if strcmp(option.Legpos,'phase'), multiLegend(data); end
grid on;
end
