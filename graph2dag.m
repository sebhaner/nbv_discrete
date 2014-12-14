function [Gdag edgeind] = graph2dag(G,node,mode)


d = graphshortestpath(G,node);

[i j v] = find(G);

keep = false(nnz(G),1);

for k=1:nnz(G)
	switch mode
		case 'toward'
			if d(j(k)) < d(i(k)) % Keep edges strictly reducing the distance to target node
				keep(k) = true;
			end
		case 'from'
			if d(j(k)) > d(i(k)) % Keep edges strictly increasing the distance from start node
				keep(k) = true;
			end
		otherwise
			error('"mode" should be "from" or "toward"');
	end
end

Gdag = sparse(i(keep),j(keep),v(keep),size(G,1),size(G,2));

edgeind = 1:nnz(G);
edgeind = edgeind(keep);