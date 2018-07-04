function out = str2rgb(str)
%str2rgb - Output rgb space from string
%
% out = str2rgb(str)
% color list
%      'b': blue 
%      'g': green 
%      'r': red
%      'k': black 
%      'm': magenta 
%      'c': cyan 
%     'g2': deep blue green 
%     'b2': sky blue
%     'b3': deepskyblue
% Author  : Wataru Ohnishi, University of Tokyo, 2017

% colorlist = {'b','g','r','m','k','c','g2','b2','b3'};
if strcmp(str, 'b')
    out = [0,0,1];
elseif strcmp(str, 'g')
    out = [0,0.5,0];
elseif strcmp(str, 'r')
    out = [1,0,0];
elseif strcmp(str, 'k')
    out = [0,0,0];
elseif strcmp(str, 'm')
    out = [1,0,1];
elseif strcmp(str, 'c')
    out = [0,1,1];
elseif strcmp(str, 'g2')
    out = [0,187,193]/255;
    % http://www.colordic.org/picker/
    % http://ironodata.info/rgb.php?color=00BBC1
elseif strcmp(str, 'b2')
    out = [135,206,235]/255;
    % sky blue
    % http://members2.jcom.home.ne.jp/3264mvsj/lecture/color.html
elseif strcmp(str, 'b3')
    out = [0,191,255]/255;
    % deepskyblue
    % http://members2.jcom.home.ne.jp/3264mvsj/lecture/color.html
elseif length(str) == 3
    out = str;
else
    error('str = ''r'',''g'',''b'',''k'',''m'', ''c'', ''g2'', ''b2'', ''b3'' or RGB space');
end

end
