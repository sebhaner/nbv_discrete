function p = solveLinearApprox(p)

% Do a "linear approximation" by maximizing the trace
% instead of minimizing the trace-of-inverse:
%
% min. c'*x - trace(sum(I_i*alpha_i))/lambda
% s.t. A_G*x == b
%	  0 <= x <= 1
%
% This is equivalent to the shortest path problem on
% the graph with modified weights, which may be negative.
% If there are no negative cycles it can be solved using
% e.g. the Bellman-Ford algorithm. With the LP formulation,
% we can still extract a feasible path even with negative
% cycles present. Unfortunately, with negative cycles the
% problem is still NP-hard if we enforce the no-loop
% constraints.
%
% Alternate formulation is
% min. -trace(sum(I_i*alpha_i))
% s.t. A_G*x == b
%		c'*x <= L
%		x in {0,1}^n

[type,param] = parseProblem(p);


[x alpha obj] = solveLP(p.G,p.source,p.dest,p.info,param,type,@(xx)plotfun(xx,p.G,p.x,p.y,p.z));


p.lpsol.x = x;
p.lpsol.alpha = alpha;
p.lpsol.lpobj = obj;
p.lpsol.path = walkPath(p.G,round(x),p.source,p.dest);
p.lpsol.obj = evalPath(p.lpsol.path,p.G,p.info,param,type,0);


% plotEdges(p.lpsol.x,p.G,p.x,p.y,p.z);



function [x alpha obj] = solveLP(G,source,dest,info,param,type,plotfun)

verbosity = 2;

n = size(G,1);
[row,col,c] = find(G);
ne = nnz(G);
nop = size(info{1},3);

A = edgeincidence(G);
b = zeros(n,1);
b(source) = -1;
b(dest) = 1;

A2 = A; A2(A2==-1) = 0;  % incident edges
A3 = A; A3(A3==1) = 0; A3 = -A3; % exiting edges
ind_ns = setdiff(1:n,source);
A(dest,:) = []; % one constraint is redundant and can be removed
b(dest) = [];

tt = tic;

% Precompute the trace of I_i
tr = zeros(n,1);
for i=1:n
	for j=1:nop
		tr(i) = tr(i) + trace(info{i+1}(:,:,j));
	end
end

% Init YALMIP variables
yalmip('clear');
x = binvar(ne,1);
alpha = sdpvar(n,1);
% Objective function
if type==1 || type==3
	f = sum(c.*x) - sum(alpha.*tr)/param;
	Clen = [];
else
	f = -sum(alpha.*tr);
	Clen = sum(c'*x) <= param;
end

% Linear graph constraints
C0 = A*x == b; % Sum incident/exiting edges = 0, except -1 for source node, 1 for terminal
C1 = A2(ind_ns,:)*x == alpha(ind_ns); % alfa_i = sum of incident edges
C2 = A3(source,:)*x == alpha(source); % alfa_source = sum of exiting edges
C3 = alpha <= 1; % Pass a node at most once

C = C0 + C1 + C2 + C3 + Clen + (0 <= x <= 1); %#ok

if verbosity, fprintf('YALMIP parser: %.2f s\n',toc(tt)); end

% Iterate while adding no-loop constraints
Ccycle = [];
for iter=1:30
	solvesdp(C+Ccycle,f,sdpsettings('verbose',0));
	dx = round(double(x));
	path = walkPath(G,dx,source,dest,1);
	[val len] = evalPath(path,G,info,param,type,0);
	fprintf('Iteration %u: %.2f (%.2f)\n',iter,val,len);
	
	
% 	if nargin>5, plotfun(dx); drawnow; end
	
	% Find cycles in solution
	H = sparse(row,col,dx,size(G,1),size(G,2)); % solution graph
	[~,CC] = graphconncomp(H,'Weak',true);
	
	cset = setdiff(CC(A2*dx>0),CC(source)); % remove unconnected nodes and source conn. comp.
	
	if isempty(cset), break; end
	for cc=cset
		Csum = 0;
		for j=find(CC==cc) % accumulate edges in cycle
			Csum = Csum + x(edgelinearindex(G,j,find(H(j,:))));
		end
		Ccycle = [Ccycle, Csum <= nnz(CC==cc)-1]; %#ok
	end
	
	
end

x = double(x);
alpha = double(alpha);
obj = double(f);



function plotfun(xx,G,x,y,z)
[r c] = find(G);
plot3(x(:),y(:),z(:),'k.','markersize',1);
hold on;
for i=find(xx(:)')
	line([x(r(i)) x(c(i))],[y(r(i)) y(c(i))],[z(r(i)) z(c(i))]);
end
hold off;

