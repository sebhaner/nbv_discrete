function p = runReducedGA(p,varargin)

initialize = 0;
maxiter = 100000;
T = 200;
nodes = [];

for i=1:2:length(varargin)
	switch varargin{i}
		case 'maxiter', maxiter = varargin{i+1};
		case 'init', initialize = varargin{i+1};
		case 'nodes', nodes = varargin{i+1};
		case 'T', T = varargin{i+1};
		otherwise, error('Parameter parse error');
	end
end


if ~isfield(p,'sdpsol') || isempty(p.sdpsol.alpha), error('Run "solvePlanningSDP" first'); end
if ~isfield(p,'rsol')
	p.rsol = [];
end
if ~isfield(p.rsol,'beta') || isempty(p.rsol.beta) || initialize || (~isempty(nodes) && nodes~=p.rsol.nodes)
	fprintf('Rounding SDP solution...\n');
	if ~isempty(nodes)
		p.rsol.nodes = nodes;
	end
	if ~isfield(p.rsol,'nodes') || isempty(p.rsol.nodes)
		error('Specify number of nodes');
	end
	p.rsol.beta = roundSDPsol(p);
end
if ~p.rsol.beta(p.source) || ~p.rsol.beta(p.dest)
	error('Start and destination nodes not in alpha');
end

%%%%%%
alpha = p.rsol.beta;

n = nnz(alpha);
nodeind = find(alpha);
source = find(nodeind==p.source);
dest = find(nodeind==p.dest);
% Construct a reduced graph containing only the selected nodes
G = p.G(alpha,alpha);
% Set distances to shortest path dist. in full graph
sp = cell(n,n);
for i=1:n
	[d s] = graphshortestpath(p.G,nodeind(i),nodeind);
	d(isinf(d)) = 0;
	G(i,:) = d; % Not actually used right now
	sp(i,:) = s;
end


[type param] = parseProblem(p);


npop1 = 50;
npop2 = 20;
pcross = 0.1;
pmutate = 0.3;
pmutate2 = 0.2;
pmutate3 = 0.2;
pmutate4 = 0.4;
penalty_weight = 10000;
iterations = maxiter;


% Initialize population
Gt = G';
npop = npop1;
pop = cell(npop,1);
for i=1:npop
% 	pop{i} = lerw_mex(Gt,source,dest);
	pop{i} = dfs_rec_mex(Gt,source,dest);
end


% sfigure(2)
tbest = 0;
plot3(p.x(:),p.y(:),p.z(:),'k.','markers',1); hold on;
plot3(p.x(alpha),p.y(alpha),p.z(alpha),'k.','markers',5);
plot3(p.X(1,:),p.X(2,:),p.X(3,:),'r.');
h = plot3(0,0,0,'ms');
h2 = plot3(0,0,0,'b');
m = median([p.x(~isnan(p.x))' p.y(~isnan(p.y))' p.z(~isnan(p.z))'],1);
camtarget(m); campos(m+[0 -50 0]); camup([0 0 -1]);
axis equal;
axis off;
hold off;

cost_history = [];

fprintf('cost          path length   geometric cost\n');
fprintf('------------------------------------------\n');
last_change = 1;
for iter=1:iterations
	
	% Terminate if stagnated
	if iter-last_change > T, break; end
	
	% Local refinement of best path
	if mod(iter,5)==0
		pop{1} = refinePath(pop{1},G,Gt,@(pp)evalReducedPath(pp,p,param,type,penalty_weight,sp));
	end
	if mod(iter,5)==0
		n = ceil(rand*length(pop));
		pop{n} = refinePath(pop{n},G,Gt,@(pp)evalReducedPath(pp,p,param,type,penalty_weight,sp));
	end
	
	% Evaluate population
	fitness = 1./cellfun(@(pp)evalReducedPath(pp,p,param,type,penalty_weight,sp),pop);
	
	% Plot best
	[best ind] = max(fitness);
	bp = pop{ind};
	if best>tbest
		[~,~,~,tp] = evalReducedPath(bp,p,param,type,penalty_weight,sp);
		set(h2,'XData',p.x(tp),'YData',p.y(tp),'ZData',p.z(tp));
		set(h,'XData',p.x(nodeind(bp)),'YData',p.y(nodeind(bp)),'ZData',p.z(nodeind(bp)));  drawnow
		tbest = best;
		last_change = iter;
	end
	if mod(iter,10)==0
		[val len geo] = evalReducedPath(bp,p,param,type,penalty_weight,sp);
		fprintf('%-13.3f %-13.3f %-13.3f\n',val,len,geo);
	end
	% Select new population via roulette wheel
	npop = round(npop1/(iter*0.1) + npop2);
	newpop = cell(npop,1);
	cs = cumsum(fitness);
	for i=1:npop
		r1 = rand; r2 = rand;
		p1 = pop{find(cs>r1*cs(end),1)};
		p2 = pop{find(cs>r2*cs(end),1)};
		if rand<pcross, p1 = crosspaths(p1,p2); end
		if rand<pmutate, p1 = mutatepath(p1,Gt); end
		if rand<pmutate2, p1 = mutatepath2(p1,G); end
		if rand<pmutate3, p1 = mutatepath3(p1,G); end
		if rand<pmutate4, p1 = mutatepath4(p1,G,Gt); end
		newpop{i} = p1;
	end
	pop = newpop;
	% Elitism
	pop{1} = bp;
	
	pause(0.001);

	cost_history(end+1) = 1/best; %#ok
end

[val,~,~,tp] = evalReducedPath(bp,p,param,type,penalty_weight,sp);
p.rsol.obj = val;
p.rsol.path = tp;
p.rsol.cost_history = cost_history;


function [val len geo tp] = evalReducedPath(path,p,param,type,penalty_weight,sp)

tp = [];
for i=1:length(path)-1
	tp = [tp sp{path(i),path(i+1)}(1:end-1)];
end
tp = [tp sp{path(end-1),path(end)}(end)];

[val len geo] = evalPath_mex(tp,p.G,p.info,param,type,penalty_weight);
if length(unique(tp))<length(tp), val = inf; end

