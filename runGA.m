function p = runGA(p,varargin)

maxiter = 50000;
initialize = 1;
T = 200; % stop if more than T iterations without improvement

for i=1:2:length(varargin)
	switch varargin{i}
		case 'maxiter', maxiter = varargin{i+1};
		case 'init', initialize = varargin{i+1};
		case 'T', T = varargin{i+1};
		otherwise, error('Parameter parse error');
	end
end

[type param] = parseProblem(p);

if type==2 || type==4
	d = graphshortestpath(p.G,p.source);
	if param < d(p.dest)
		fprintf('L is too small (%.2f/%.2f), problem is infeasible. Aborting\n',param,d(p.dest));
		return;
	end
	[pt ind] = trimProblem(p); % remove nodes which cannot be part of the solution
else
	pt = p;
end


[path obj cost_history] = runGA_2(pt,type,param,maxiter,initialize,T);

if type==2 || type==4
	% Restore trimmed problem solution to original indexing
	path_orig = ind(path);
else
	path_orig = path;
end

p.gasol.path = path_orig;
p.gasol.obj = obj;
p.gasol.cost_history = cost_history;



function [bp obj cost_history] = runGA_2(p,type,param,iterations,initialize,T)


npop1 = 100;
npop2 = 60;
pcross = 0.7;
pmutate = 0.2;
pmutate2 = 0.1;
pmutate3 = 0.1;
pmutate4 = 0.3;
penalty_weight = 10000;


% Initialize population
Gt = p.G';
npop = npop1;
pop = cell(npop,1);
for i=1:npop
% 	pop{i} = lerw_mex(Gt,source,dest);
	pop{i} = dfs_rec_mex(Gt,p.source,p.dest);
end

if initialize
	% Seed with previous SDP solution if available
	if isfield(p,'sdpsol') && isfield(p.sdpsol,'path') && ~isempty(p.sdpsol.path)
		pop{1} = p.sdpsol.path(:)';
	end
	% Seed with previous LP solution if available
	if isfield(p,'lpsol') && isfield(p.lpsol,'path') && ~isempty(p.lpsol.path)
		pop{2} = p.lpsol.path(:)';
	end
	% Seed best previous GA solution if available
	if isfield(p,'gasol') && isfield(p.gasol,'path')  && ~isempty(p.gasol.path)
		pop{3} = p.gasol.path(:)';
	end
	% Seed best reduced GA solution if available
	if isfield(p,'rsol') && isfield(p.rsol,'path')  && ~isempty(p.rsol.path)
		pop{4} = p.rsol.path(:)';
	end
end

% sfigure(2)
tbest = 0;
plot3(p.x(:),p.y(:),p.z(:),'k.','markers',1); hold on;
plot3(p.X(1,:),p.X(2,:),p.X(3,:),'r.');
h = plot3(0,0,0,'b');
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
		pop{1} = refinePath(pop{1},p.G,Gt,@(pp)evalPath_mex(pp,p.G,p.info,param,type,penalty_weight));
	end
	if mod(iter,1)==0
		n = ceil(rand*length(pop));
		pop{n} = refinePath(pop{n},p.G,Gt,@(pp)evalPath_mex(pp,p.G,p.info,param,type,penalty_weight));
	end
	
	% Evaluate population
	fitness = 1./cellfun(@(pp)evalPath_mex(pp,p.G,p.info,param,type,penalty_weight),pop);
	
	% Plot best
	[best ind] = max(fitness);
	bp = pop{ind};
	if best>tbest
		set(h,'XData',p.x(bp),'YData',p.y(bp),'ZData',p.z(bp));  drawnow
		tbest = best;
		last_change = iter;
	end
	if mod(iter,10)==0
		[val len geo] = evalPath_mex(bp,p.G,p.info,param,type,penalty_weight);
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
		if rand<pmutate2, p1 = mutatepath2(p1,p.G); end
		if rand<pmutate3, p1 = mutatepath3(p1,p.G); end
		if rand<pmutate4, p1 = mutatepath4(p1,p.G,Gt); end
		newpop{i} = p1;
	end
	pop = newpop;
	% Elitism
	pop{1} = bp;
	
	pause(0.001);

	cost_history(end+1) = 1/best; %#ok
end

obj = evalPath_mex(bp,p.G,p.info,param,type,penalty_weight);
