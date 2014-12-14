% Given a graph G and binary edge indicators x,
% compute the vertices visited on the path from
% source to dest in the correct order.
function path = walkPath(G,x,source,dest,quiet)


[i j] = find(G);

P = sparse(i,j,x,size(G,1),size(G,2));

[~,pred] = graphtraverse(P,source,'method','dfs');
path = graphpred2path(pred,dest);


if nargin<5, quiet = 0; end
if ~quiet && (numel(path)~=nnz(x)+1 || path(end)~=dest)
	warning('Input is not a simple path from source to dest')
end