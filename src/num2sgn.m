function s = num2sgn(val, n)
% num2sgn - convert number to string consdering significant digits
%  val    - input value (number)
%  n      - number of significant digits (>0)
% Author    : Wataru Ohnishi, University of Tokyo, 2017
%%%%%

    s = num2str(val,n);

    if isempty(strfind(s,'e')) == 0 % %e style
        if abs(val) > 1
           s = str2double(s);
           s = sprintf('%d',s);
        else
            n2 = sprintf('%d',abs(ceil(log10(val)))+n);
            s = str2double(s);
            strn = strcat('s = sprintf(''%.', num2str(n2), 'f'',s);');
            eval(strn);
        end
    end
end
