function [p ind edgeind] = trimProblem(p)

if isempty(p.costfun.L)
% 	error('Only applicable to constrained length formulation; specify L');
	return;
end

% Any node n s.t. d(source,n) + d(n,dest) > L can be removed
ds = graphshortestpath(p.G,p.source);
dt = graphshortestpath(p.G',p.dest);

rem = sum([ds; dt])>p.costfun.L;

ind = 1:size(p.G,1);
ind(rem) = [];

source_orig = p.source;
dest_orig = p.dest;

p.source = find(ind==p.source);
p.dest = find(ind==p.dest);

[i j] = find(p.G);
H = sparse(i,j,1:nnz(p.G),size(p.G,1),size(p.G,2));
H(rem,:) = [];
H(:,rem) = [];
[~,~,edgeind] = find(H);

p.G(rem,:) = [];
p.G(:,rem) = [];

p.x(rem) = [];
p.y(rem) = [];
p.z(rem) = [];

p.P(rem) = [];
p.Pshow(rem) = [];

rem = find(rem);
p.info(rem+1) = []; % info{1} is I_0, need an offset

% Re-index old solution paths to fit new graph

if isfield(p,'gasol') && ~isempty(p.gasol.path)
	if ~isempty(intersect(rem,p.gasol.path)) || p.gasol.path(1)~=source_orig || p.gasol.path(end)~=dest_orig
		fprintf('Old GA solution not valid\n');
		p.gasol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.gasol.path = revind(p.gasol.path);
		assert(~any(isnan(p.gasol.path)));
	end
end

if isfield(p,'sdpsol') && ~isempty(p.sdpsol.path)
	if ~isempty(intersect(rem,p.sdpsol.path)) || p.sdpsol.path(1)~=source_orig || p.sdpsol.path(end)~=dest_orig
		fprintf('Old SDP path not valid\n');
		p.sdpsol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.sdpsol.path = revind(p.sdpsol.path);
		assert(~any(isnan(p.sdpsol.path)));
	end
end

if isfield(p,'lpsol') && ~isempty(p.lpsol.path)
	if ~isempty(intersect(rem,p.lpsol.path)) || p.lpsol.path(1)~=source_orig || p.lpsol.path(end)~=dest_orig
		fprintf('Old LP path not valid\n');
		p.lpsol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.lpsol.path = revind(p.lpsol.path);
		assert(~any(isnan(p.lpsol.path)));
	end
end

if isfield(p,'rsol') && ~isempty(p.rsol.path)
	if ~isempty(intersect(rem,p.rsol.path)) || p.rsol.path(1)~=source_orig || p.rsol.path(end)~=dest_orig
		fprintf('Old reduced GA path not valid\n');
		p.rsol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.rsol.path = revind(p.rsol.path);
		assert(~any(isnan(p.rsol.path)));
	end
end

if isfield(p,'rhcsol') && ~isempty(p.rhcsol.path)
	if ~isempty(intersect(rem,p.rhcsol.path)) || p.rhcsol.path(1)~=source_orig || p.rhcsol.path(end)~=dest_orig
		fprintf('Old RHC path not valid\n');
		p.rhcsol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.rhcsol.path = revind(p.rhcsol.path);
		assert(~any(isnan(p.rhcsol.path)));
	end
end