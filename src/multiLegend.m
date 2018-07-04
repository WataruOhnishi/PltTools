function multiLegend(data,pre,post)
if nargin < 2
    pre = [];
end
if nargin < 3
    post = [];
end

N = length(data);
slegends = cell(1,N);
for kk = 1:N
    slegends{kk} = data{kk}.name;
end

if ~isempty(pre)
    slegends = {pre,slegends{1:end}};
end
if ~isempty(post)
    slegends = {slegends{1:end},post};
end

legend(slegends);

end       
        