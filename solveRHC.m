function p_orig = solveRHC(p,varargin)

% Receding horizon control setup:
% Limit each stage to the nodes reachable within
% k steps. Edges to nodes outside this set are replaced
% by edges to the destination node with weights equal
% to the shortest distance there in the original graph.
% The corresponding I at the edge is set to the accumulated
% information expected along that shortest path.
% This will not work very well with a free start supernode
% since then every node is within one step.

steps = 1;
solver = 'ga';

for i=1:2:length(varargin)
	switch varargin{i}
		case 'steps', steps = varargin{i+1};
		case 'solver', solver = varargin{i+1};
		otherwise, error('Parameter parse error');
	end
end

type = parseProblem(p);
if type~=2 && type~=4
	error('L must be specified for RHC');
end

tpath = p.source;
nop = size(p.info{1},3);
L = p.costfun.L;
p_orig = p;

for iter=1:50
	
	[pt ind] = trimProblemRHC(p,steps);
	pt.costfun.L = pt.costfun.L*0.5;
	switch solver
		case 'bnb'
			pt = solvePlanningSDP(pt,'bnb',1);
			path = ind(pt.sdpsol.path(1:end-1));
		case 'ga'
			pt = runGA(pt,'maxiter',30,'T',inf);
			path = ind(pt.gasol.path(1:end-1));
		otherwise
			error('"solver" must be "bnb" or "ga"');
	end

	% Add accumulated information to I_0 (the path traveled so far)
	for n=path(1:end-1)
		for pp=1:nop
			p.info{1}(:,:,pp) = p.info{1}(:,:,pp) + p.info{n+1}(:,:,pp);
		end
	end
	
	tpath = [tpath path(2:end)]; %#ok
	p.source = path(end);
	
	if ~isempty(L)
		[~,len] = evalPath(tpath,p.G,p.info,0,1); % just get the path length
		p.costfun.L = L - len;
	end
	
	if length(path)==1
		tpath(end+1) = p.dest; %#ok
		break;
	end
		
	sfigure(1);
	p.rhcsol.path = tpath;
	plotProblem(p,'rhcpath',1,'edges',0);

end

p_orig.rhcsol.path = tpath;
[val len] = evalPath(tpath,p_orig.G,p_orig.info,L,type,0);
p_orig.rhcsol.obj = val;
p_orig.rhcsol.length = len;

sfigure(1);
plotProblem(p_orig,'rhcpath',1,'edges',0);