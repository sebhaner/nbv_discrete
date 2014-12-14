function path = dfs_rec(Gt,source,dest)

n = size(Gt,1);
% Gt = G';
pred = NaN(n,1);
pred(source) = 0;

rec(source);
pred(pred==-1) = NaN;
path = graphpred2path(pred',dest);

	function r = rec(node)
		if node==dest, r = true; return; end
		neigh = find(Gt(:,node) & isnan(pred));
		ind = randperm(length(neigh));
		for i=ind
			pred(neigh(i)) = node;
			if rec(neigh(i)), r = true; return; end
			pred(neigh(i)) = -1; % Do not revisit, takes too long
		end
		r = false;
	end

end
			