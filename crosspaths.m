function p = crosspaths(p1,p2)

nodes = max([p1 p2])+2;
n = length(p1)+length(p2)+1;
indi = zeros(n,1);
indj = zeros(n,1);

q = 1;
for i=1:length(p1)-1
	indi(q) = p1(i);
	indj(q) = p1(i+1);
	q = q+1;
end
for i=1:length(p2)-1
	indi(q) = p2(i);
	indj(q) = p2(i+1);
	q = q+1;
end

indi(q) = nodes-1; indj(q) = p1(1);
indi(q+1) = nodes-1; indj(q+1) = p2(1);
indi(q+2) = p1(end); indj(q+2) = nodes;
indi(q+3) = p2(end); indj(q+3) = nodes;

Gt = sparse(indj,indi,ones(q+3,1),nodes,nodes);

% [~,pred] = graphminspantree(sprand(G+G'),nodes-1);
% p = graphpred2path(pred,nodes);
% p = lerw_mex(Gt,nodes-1,nodes);
p = dfs_rec_mex(Gt,nodes-1,nodes);

p = p(2:end-1);

if isempty(p), keyboard; end