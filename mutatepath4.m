function p = mutatepath4(p,G,Gt)

if length(p)<3, return; end

% Select one node on path, and the next one 2 steps ahead
n = ceil(rand*(length(p)-2))+1;

% Switch node to a common neighbour
neigh = intersect(find(Gt(:,p(n-1))),find(G(:,p(n+1))));
if isempty(neigh), return; end
new = neigh(ceil(rand*length(neigh)));

% Insert into path
p = [p(1:n-1) new p(n+1:end)];


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

