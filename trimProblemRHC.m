function [p ind edgeind] = trimProblemRHC(p,steps)

d = graphshortestpath(spones(p.G),p.source); % Get distance as number of steps
rem =  setdiff(find(d > steps),p.dest); % Keep destination node

% Redefine info at edge nodes as sum of info along shortest path
% to destination and add corresponding arcs to destination.
% All arcs from edge nodes back into the active set are also removed.

nop = size(p.info{1},3);
edgenodes = find(d==steps);
for n=edgenodes(:)'
	[dist,path] = graphshortestpath(p.G,n,p.dest);
	for pp=path(2:end-1)
		for i=1:nop
			p.info{n+1}(:,:,i) = p.info{n+1}(:,:,i) + p.info{pp+1}(:,:,i);
		end
	end
	p.G(n,:) = 0;
	p.G(n,p.dest) = dist;
end

ind = 1:size(p.G,1);
ind(rem) = [];

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

p.info(rem+1) = []; % info{1} is I_0, need an offset


% Re-index old solution paths to fit new graph

if isfield(p,'gasol') && ~isempty(p.gasol.path)
	if ~isempty(intersect(rem,p.gasol.path)) || p.gasol.path(1)~=p.source || p.gasol.path(end)~=p.dest
% 		fprintf('Old GA solution not valid\n');
		p.gasol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.gasol.path = revind(p.gasol.path);
		assert(~any(isnan(p.gasol.path)));
	end
end

if isfield(p,'sdpsol') && ~isempty(p.sdpsol.path)
	if ~isempty(intersect(rem,p.sdpsol.path)) || p.sdpsol.path(1)~=p.source || p.sdpsol.path(end)~=p.dest
% 		fprintf('Old SDP path not valid\n');
		p.sdpsol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.sdpsol.path = revind(p.sdpsol.path);
		assert(~any(isnan(p.sdpsol.path)));
	end
end

if isfield(p,'lpsol') && ~isempty(p.lpsol.path)
	if ~isempty(intersect(rem,p.lpsol.path)) || p.lpsol.path(1)~=p.source || p.lpsol.path(end)~=p.dest
% 		fprintf('Old LP path not valid\n');
		p.lpsol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.lpsol.path = revind(p.lpsol.path);
		assert(~any(isnan(p.lpsol.path)));
	end
end

if isfield(p,'rsol') && ~isempty(p.rsol.path)
	if ~isempty(intersect(rem,p.rsol.path)) || p.rsol.path(1)~=p.source || p.rsol.path(end)~=p.dest
% 		fprintf('Old reduced GA path not valid\n');
		p.rsol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.rsol.path = revind(p.rsol.path);
		assert(~any(isnan(p.rsol.path)));
	end
end

if isfield(p,'rhcsol') && ~isempty(p.rhcsol.path)
	if ~isempty(intersect(rem,p.rhcsol.path)) || p.rhcsol.path(1)~=p.source || p.rhcsol.path(end)~=p.dest
% 		fprintf('Old RHC path not valid\n');
		p.rhcsol.path = [];
	else
		revind = nan(size(p.G,1),1);
		revind(ind) = 1:length(ind);
		p.rhcsol.path = revind(p.rhcsol.path);
		assert(~any(isnan(p.rhcsol.path)));
	end
end