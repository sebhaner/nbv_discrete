function p = refinePath(p,G,Gt,costfun)

% Select vertices on the path in random order, switch to neighbors,
% select change which reduces cost most in every step
cost = costfun(p);
ind = randperm(length(p)-2)+1;
for n=ind
	neigh = setdiff(intersect(find(Gt(:,p(n-1))),find(G(:,p(n+1)))),p);
	for i=1:length(neigh)
		np = [p(1:n-1) neigh(i) p(n+1:end)];
		nc = costfun(np);
		if nc<cost, cost = nc; p = np; end
	end	
end

% p = removeLoops(p);