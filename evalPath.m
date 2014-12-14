% type 1: formulation 1, trace
% type 2: formulation 2, trace
% type 3: formulation 1, maximum eigenvalue
% type 4: formulation 2, maximum eigenvalue
function [val len geo] = evalPath(path,G,info,lambda,type,penalty_weight)

if nargin<6, penalty_weight = 1e6; end

len = 0;
for k=1:length(path)-1
	v = G(path(k),path(k+1));
	if v==0, error('Not a valid path, no edge (%u,%u)',path(k),path(k+1)); end
	len = len + v;
end


geo = 0;
for p=1:size(info{1},3)
	I = info{1}(:,:,p); % I0
	for kk=path(:)', I = I + info{kk+1}(:,:,p); end
	if type<=2
		geo = geo + trace(inv(I));
	else
		geo = max(geo,1/min(eig(I)));
	end
end

if type==1 || type==3
	val = len + geo/lambda;
else
	penalty = penalty_weight*max(0,len-lambda)^2;
	val = geo + penalty;
end


