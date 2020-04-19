function [pfig, hfig, ax] = pltBode(in,option,cohflag)
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
%   cohflag         : coherence plot
% Advanced options
%   data1.P         : plant model
%   data1.Cfb       : FB controller model
%   data1.Csh       : shaping filters
%   option.plot     : 'P'   : plant
%                     'PCsh': plant shaping (PCsh)
%                     'S'   : sensitivity (1/(1+PCfbCsh))
%                     'T'   : complementary sensitivity (PCfbCsh/(1+PCfbCsh))
%                     'Gop' : open loop (PCfbCsh)
% Author    : Wataru Ohnishi, University of Tokyo, 2017
%%%%%


N = length(in);

if nargin < 2
    option = struct;
end

if nargin < 3
    cohflag = 0;
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

try option.plot; catch, option.plot = 'sys'; end
for k = 1:1:N
    if strcmp(option.plot,'P'), data{k}.sys = data{k}.P; end
    if strcmp(option.plot,'PCsh'), data{k}.sys = data{k}.P*data{k}.Csh; end
    if strcmp(option.plot,'S'), data{k}.sys = feedback(1,data{k}.P*data{k}.Cfb*data{k}.Csh); end
    if strcmp(option.plot,'T'), data{k}.sys = feedback(data{k}.P*data{k}.Cfb*data{k}.Csh,1); end
    if strcmp(option.plot,'Gop'), data{k}.sys = data{k}.P*data{k}.Cfb*data{k}.Csh; end
end

if ~isfield(option,'fmin'), option.fmin = 1; option.fmax = 1000; end
if ~isfield(option,'pmin'), option.pmin = -360; option.pmax = 0; end
if ~isfield(option,'ptick'), option.ptick = 90; end
if ~isfield(option,'foption'), option.foption = 'log'; end
if ~isfield(option,'Legpos'), option.Legpos = 'gain'; end
if ~isfield(option,'LegendLoc'), option.LegendLoc = 'best'; end
if ~isfield(option,'phasePlot'), option.phasePlot = 1:length(in); end
if ~isfield(option,'freq'), freq = logspace(log10(option.fmin),log10(option.fmax),1000); else, freq = option.freq; end

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

ax = nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hfig = figure;
if ~cohflag, ax(1) = subplot(2,1,1);
elseif cohflag, ax(1) = subplot(3,1,1); end
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
    end
end
axis([option.fmin, option.fmax, option.gmin, option.gmax]);
set(gca,'ytick',option.gmin:option.gtick:option.gmax);
ylabel('Magnitude [dB]');
if strcmp(option.Legpos,'gain'), multiLegend(data); end
grid on; hold on;
try title(option.title); catch, end


if ~cohflag, ax(2) = subplot(2,1,2);
elseif cohflag, ax(2) = subplot(3,1,2); end
for k = option.phasePlot
    phasedeg = rad2deg(angle(squeeze(data{k}.sys.ResponseData)));
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
            set(h,'MarkerEdgeColor',data{k}.color);
            set(h,'MarkerFaceColor',data{k}.color);
            set(gca,'xscale','log');
        else
            h = semilogx(data{k}.sys.frequency,phasedeg); hold on;
            try set(h,'Color',data{k}.color); catch, end
            try set(h,'linestyle',data{k}.style); catch, end
            try set(h,'Marker',data{k}.marker); catch, end
        end
    else
        h = plot(data{k}.sys.frequency,phasedeg); hold on;
    end
    
end
axis([option.fmin, option.fmax, option.pmin, option.pmax]);
set(gca,'ytick',option.pmin:option.ptick:option.pmax);
if ~cohflag, xlabel('Frequency [Hz]'); end
ylabel('Phase [deg]');
if strcmp(option.Legpos,'phase'), multiLegend(data); end
grid on; hold on;



if cohflag, ax(3) = subplot(3,1,3);
    for k = 1:1:N
        if ~isempty(data{k}.sys.UserData)
            if strcmp(option.foption, 'log')
                if strcmp(data{k}.style,'')
                    hc = scatter(data{k}.sys.frequency,data{k}.sys.UserData,'filled','o'); hold on
                    set(hc,'MarkerEdgeColor',data{k}.color);
                    set(hc,'MarkerFaceColor',data{k}.color);
                    set(gca,'xscale','log');
                else
                    hc = semilogx(data{k}.sys.frequency,data{k}.sys.UserData); hold on;
                end
            else
                hc = plot(data{k}.sys.frequency,data{k}.sys.UserData); hold on;
            try set(hc,'Color',data{k}.color); catch, end
            try set(hc,'linestyle',data{k}.style); catch, end
            try set(hc,'Marker',data{k}.marker); catch, end
            end
        else
            % for bug for scatter
            hc = plot(data{k}.sys.freq(1)/10,2);
        end
    end
    axis([option.fmin, option.fmax, 0, 1]);
    set(gca,'ytick',0:0.2:1);
    xlabel('Frequency [Hz]');
    ylabel('Coherence [-]');
    grid on; hold on;
end

% if ~cohflag, linkaxes([ax1,ax2],'x');
% elseif cohflag, linkaxes([ax1,ax2,ax3],'x'); end
linkaxes(ax,'x');

pfig = pubfig(hfig);
pfig.LegendLoc = option.LegendLoc;
if cohflag, pfig.Dimension = [19,16]; end

