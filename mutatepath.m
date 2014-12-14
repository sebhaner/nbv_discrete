function p = mutatepath(p1,Gt)

% Select two random nodes on path
n = sort(randperm(length(p1),2));

% Find a new path over the interval
% [~,pred] = graphminspantree(sprand(Gt),p(n(1)));
% pfill = graphpred2path(pred,p(n(2)));
% pfill = lerw_mex(Gt,p(n(1)),p(n(2)));
pfill = dfs_rec_mex(Gt,p1(n(1)),p1(n(2)));
% [~,pfill] = graphshortestpath(Gt',p(n(1)),p(n(2)));

% Insert into path
p = [p1(1:n(1)-1) pfill p1(n(2)+1:end)];

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

if isempty(p), keyboard; end