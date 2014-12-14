function p = mutatepath2(p,G)

% Select two random nodes on path
% n = sort(ceil(rand(1,2)*length(p)));

% Select one node on path, and the next one l steps ahead
l = min(length(p)-1,2);
n = ceil(rand*(length(p)-l));
n = [n min(length(p)-1,n+l)];

% Find a new path over the interval
[~,pfill] = graphshortestpath(G,p(n(1)),p(n(2)));

% Insert into path
p = [p(1:n(1)-1) pfill p(n(2)+1:end)];

% Remove any loops introduced
p = removeLoops(p);
% succ = zeros(1,size(G,1));
% for i=1:length(p)-1
% 	succ(p(i)) = p(i+1);
% end
% for i=2:length(p)
% 	p(i) = succ(p(i-1));
% 	if p(i)==p(end), break; end
% end
% p = p(1:i);
