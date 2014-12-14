% Generate the edge incidence matrix for a graph G.
% Edges are in column-major order as they appear in G.
% A(i,j) =  1 if edge j enters node i
% A(i,j) = -1 if edge j exits node i
function A = edgeincidence(G)

[i,j,~] = find(G);

indi = zeros(numel(i)*2,1);
indj = zeros(numel(i)*2,1);
mval = zeros(numel(i)*2,1);

q = 1;
for k=1:numel(i)
	
% 	A(i(k),k) = A(i(k),k) - 1;
% 	A(j(k),k) = A(j(k),k) + 1;
	
	indi(q) = i(k);
	indj(q) = k;
	mval(q) = -1;
	
	indi(q+1) = j(k);
	indj(q+1) = k;
	mval(q+1) = 1;
	q = q+2;

end

A = sparse(indi,indj,mval,size(G,1),nnz(G));