function ind = getsymmetricindex(G)

% If G is a symmetric graph and the nodes are numbered as by find(G),
% which node index corresponds to x_ji given the index of x_ij?
% Returns ind s.t. ind(ind_ij) = ind_ji

[row,col,~] = find(G);

noncompactedind = (col-1)*size(G,1) + row;

ind = nan(nnz(G),1);
for n=1:nnz(G)
	
	nci = (row(n)-1)*size(G,1) + col(n);
	i = find(noncompactedind==nci);
	
	if ~isempty(i) % If G is not symmetric, keep the NaN
		ind(n) = i;
	end
	
end
	