function p = solvePlanningSDP(p,varargin)

relaxation_only = 1;
verbosity = 1;

for i=1:2:length(varargin)
	switch varargin{i}
		case 'bnb', relaxation_only = ~varargin{i+1};
		case 'verbosity', verbosity = varargin{i+1};
		otherwise, error('Parameter parse error');
	end
end

[type param] = parseProblem(p);

if type==2 || type==4
	d = graphshortestpath(p.G,p.source);
	if param < d(p.dest)
		fprintf('L is too small, problem is infeasible. Aborting\n');
		return;
	end
	[pt ind edgeind] = trimProblem(p); % remove nodes which cannot be part of the solution
else
	pt = p;
end

tic
[path xx alpha lb] = sdpbb_multi(pt.G,pt.source,pt.dest,pt.info,param,type,relaxation_only,pt.x,pt.y,pt.z,verbosity);
if verbosity, fprintf('SDP solution took %.2f s\n',toc); end

if type==2 || type==4
	% Restore trimmed problem solution to original indexing
	path_orig = ind(path);
	alpha_orig = zeros(size(p.G,1),1);
	alpha_orig(ind) = alpha;
	xx_orig = zeros(nnz(p.G),1);
	xx_orig(edgeind) = xx;
else
	path_orig = path;
	xx_orig = xx;
	alpha_orig = alpha;
end

if ~(isempty(path_orig) && isfield(p,'sdpsol') && isfield(p.sdpsol,'path')) % do not overwrite old bnb solution
	p.sdpsol.path = path_orig; % empty if branch-and-bound not run
end
p.sdpsol.x = xx_orig;
p.sdpsol.alpha = alpha_orig;
p.sdpsol.obj = lb;
if ~isempty(path_orig)
	p.sdpsol.pathobj = evalPath(path_orig,p.G,p.info,param,type,0);
else
	p.sdpsol.pathobj = [];
end