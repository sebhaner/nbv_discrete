function out = bb3(obj,C,x,alpha,t,UB)%,t,info,G,source,dest,sol,plotfun)

sol=[];
opts = sdpsettings;
% opts.usex0 = true; % Not applicable to most interior-point solvers
opts.verbose = 0;

tic
[~,~,~,problem] = export(C,obj,opts);
fprintf('\n Yalmip parser: %.2f s\n',toc);
Q = problem.Q;
f = problem.f;
c = problem.c;
nb = length(x);
problem.lb(1:nb) = 0;
problem.ub(1:nb) = 1;
problem.options.mosek.MSK_IPAR_PRESOLVE_USE = 'MSK_PRESOLVE_MODE_OFF'; % bug in Mosek 7.0.0.90 presolver

order = 'depth';
int_tol = 1e-4;
gap_tol = 0.01;
prune_tol = 0.00001;



root = struct('p',problem,'lb',-inf,'depth',0);
qs = cell(1,1); % Set of nodes
qs{1} = root;

ncount = 0;

% UB = inf;
UBpath = inf;
gap = inf;
pathgap = inf;
min_fathomed = inf;


fprintf('\n Node Open  Lower        Upper        Gap(%%)  Upper path   Gap(%%) Time(s)  Status\n');
fprintf('--------------------------------------------------------------------------------------\n');
timer = tic;
while ~isempty(qs) && (gap>gap_tol || isempty(sol))
	
	ncount = ncount + 1;
	
	% Pop and expand next node
	[node qs] = selectNextNode(qs,order);
	
	% Solve node problem
	output = callmosek(node.p);
	
	lb = computecost(f,c,Q,output.Primal,node.p);
	
	xopt = output.Primal(1:nb);
	
	if output.problem==12 || output.problem==1 || output.problem==-4
		msg = 'Infeasible problem';
		lb = inf;
	else
	
		% Is solution integer?
		dx = abs(round(xopt) - xopt);
		if all(dx < int_tol)
			if lb < UB % New incumbent?
				UB = lb;
				sol = output.Primal;
				
				if checkfeasiblefast(problem,sol,1e-4)~=1, keyboard; end
% 				plotfun(round(xopt)); drawnow
				% Prune nodes
				bad = cellfun(@(s)s.lb,qs) >= UB*(1 - prune_tol);
				qs(bad) = [];
			end
			min_fathomed = min(min_fathomed,lb);
			msg = sprintf('Integer solution, lb: %f',lb);
			order = 'best';
		else
			msg = sprintf('Solved, lb: %f',lb);
			% Find an upper bound, update UB if better
			% First, try to round it
			ispath = false;
			xx = output.Primal;
			xx(1:nb) = round(xx(1:nb));
			feasible = checkfeasiblefast(problem,xx,1e-4);
			if ~feasible % if not feasible, try something else
				% ------ Try to generate feasible solution close to relaxation --------
% 				[x_f alfa_f t_f] = findFeasiblePath(xopt,G,source,dest,info,length(t)/3,@(xx)computecost(f,c,Q,xx,problem));
% 				xx = [x_f; alfa_f; t_f];
% 				if checkfeasiblefast(problem,xx,1e-4)
% 					ispath = true;
% 					feasible = true;
% 				end
			end
			if feasible
% 				disp(xx')
				
				obj_val = computecost(f,c,Q,xx,problem);

				if ispath && obj_val < UBpath
					UBpath = obj_val;
% 					plotfun(xx(1:nb));
				end

				if  obj_val < UB % New incumbent?
					UB = obj_val;

					sol = xx;
					msg = [msg ', new UB'];

					% Prune nodes
					bad = cellfun(@(s)s.lb,qs) >= UB*(1 - prune_tol);
					qs(bad) = [];
				end
			end
			

			% Do we need to continue this branch?
			if lb < UB
	
				% Branch on next fractional variable
				nonint = find(dx>int_tol & (node.p.lb(1:nb)==0 & node.p.ub(1:nb)==1));
				[~,ind] = sort(dx(nonint),'descend');
% 				if isempty(nonint), nonint = 1; ind=1; end
				var = nonint(ind(1));
				
				p0 = node.p; p0.x0 = output.Primal; p0.ub(var) = 0;
				p1 = node.p; p1.x0 = output.Primal; p1.lb(var) = 1;
				node0 = struct('p',p0,'lb',lb,'depth',node.depth+1);
				node1 = struct('p',p1,'lb',lb,'depth',node.depth+1);
				% Put 'most promising' first in queue
				if xopt(var)>0.5, tn = node0; node0 = node1; node1 = tn; end
				qs = [qs node0 node1]; %#ok
				msg = [msg ', lb<UB, branching'];
			end
			
		end
		
	end
	
	LB = min([cellfun(@(s)s.lb,qs) lb min_fathomed]);
% 	fprintf('minopen: %f lb: %f min_fathomed: %f\n',min(cellfun(@(s)s.lb,qs)), lb, min_fathomed);
	gap = (UB-LB)/(eps+abs(LB));
	pathgap = (UBpath-LB)/(eps+abs(LB));
	
% 	if LB>UB, keyboard; end
	
	fprintf(' %-4u %-5u %-11.5e  %-11.5e  %-7.2f %-12.5e %-6.2f %-8.1f %s\n',...
		ncount,length(qs),LB,UB,gap*100,UBpath,pathgap*100,toc(timer),msg);
	
end

out = [];
if ~isempty(sol)
	out.x = round(sol(1:nb));
	out.alpha = sol(nb+1:nb+length(alpha));
	out.t = sol(nb+length(alpha)+1:nb+length(alpha)+length(t));
	out.bestUB = UB;
else
	warning('sol empty');
end

fprintf(' Final cost: %.6f\n',UB);

end

function [node qs] = selectNextNode(qs,order)
	
	switch order
		case 'depth'
			[~,ind] = max(cellfun(@(s)s.depth,qs));
			
		case 'breadth'
			[~,ind] = min(cellfun(@(s)s.depth,qs));
			
		case 'best'
			[~,ind] = min(cellfun(@(s)s.lb,qs));
		otherwise
			error('Unknown order')
	end
	
	node = qs{ind};
	qs(ind) = [];
end
