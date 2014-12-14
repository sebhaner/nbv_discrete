function [path x alpha lb] = sdpbb_multi(G,source,dest,info,param,type,relaxation_only,px,py,pz,verbosity)

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
A(dest,:) = []; % one constraint is redundant and could be removed
b(dest) = [];

tt = tic;
solver = 0; % 0: YALMIP bnb solver, 1 own.

% Init YALMIP variables
yalmip('clear');
fprintf('Parsing problem...\n');
if solver==0 && ~relaxation_only
	x = binvar(ne,1);
else
	x = sdpvar(ne,1);
end
alpha = sdpvar(n,1);
if type<=2 % trace
	t = sdpvar(3*nop,1);
else
	t = sdpvar(1,1); % maximum eigenvalue
end

% Semidefinite constraints
Ctr = [];
for p=1:nop
	I = info{1}(:,:,p);
	for i=1:n, I = I + info{i+1}(:,:,p)*alpha(i); end
	if type<=2
		Ctr = Ctr + ([I [1 0 0]'; [1 0 0] t(3*p-2)] >= 0); % tr(inv(sum(alfa_i*info_i))) <= t
		Ctr = Ctr + ([I [0 1 0]'; [0 1 0] t(3*p-1)] >= 0); % implies sum(t) >= the trace
		Ctr = Ctr + ([I [0 0 1]'; [0 0 1] t(3*p)  ] >= 0);
	else
		Ctr = Ctr + (I + t*eye(3) >= 0); %  min(eig(I)) + t >= 0 <=> max(eig(inv(I))) <= -1/t (we may assume t<=0 since I psd)
	end
end

% Objective function
Cp = [];
switch type
	case 1  % Trade-off formulation, trace
		f = sum(c.*x) + sum(t)/param;
	
	case 2  % Alternate formulation, trace
		f = sum(t);
		Cp = sum(c.*x) <= param;
	
	case 3 % Trade-off formulation, max e.v.
		u = sdpvar(1,1);
		f = sum(c.*x) + u/param;
		Cp = [-t 1; 1 u] >= 0; % -1/t <= u
	
	case 4  % Alternate formulation, max e.v.
		u = sdpvar(1,1);
		f = u; % Could just use t but nice to get the actual objective value
		Cp = ([-t 1; 1 u] >= 0) + (sum(c.*x) <= param);
end

% Linear graph constraints
C0 = A*x == b; % Sum incident/exiting edges = 0, except -1 for source node, 1 for terminal
C1 = A2(ind_ns,:)*x == alpha(ind_ns); % alfa_i = sum of incident edges
C2 = A3(source,:)*x == alpha(source); % alfa_source = sum of exiting edges
C3 = alpha <= 1; % Pass a node at most once

C = C0 + C1 + C2 + C3 + Ctr + Cp;

% Attempt better relaxation (but much slower)
% if relaxation_only
% 	y = 2*[x; alpha]-1;
% 	H = sdpvar(ne+n+1);
% 	Ch = (H>=0) + (H(1,:)==[1 y']) + (diag(H)==1);
% 	C = C + Ch;
% end

% Add linear 2-cycle constraints
symind = getsymmetricindex(G);
C7 = []; for i=1:ne, si = symind(i); if ~isnan(si), C7 = C7 + (x(i)+x(si) <= 1); end; end; C = C + C7;

if verbosity, fprintf('YALMIP parser: %.2f s\n',toc(tt)); end

if relaxation_only
	solvesdp(C+(0<=x<=1),f,sdpsettings('verbose',verbosity));
	path = [];
	x = double(x);
	alpha = double(alpha);
	lb = double(f);
	if verbosity, fprintf('Objective value: %g\n',lb); end
% 	if type>=3, disp([double(u) double(t) double(abs(-1/t - u))]); end
	return;
end


opts = sdpsettings('solver','bnb');
opts.bnb.gaptol = 0.01;
opts.bnb.presolve = 0;
opts.bnb.method = 'depthbest';
opts.bnb.branchrule = 'max';
opts.bnb.maxiter = 1000;

UB = inf;

Ccycle = [];

% solvesdp(C + Ccycle + [0<=x<=1],f,sdpsettings('dualize',1,'removeequalities',1))

for iter=1:100
	
	if solver==0
		solvesdp(C + Ccycle,f,opts);
	else
		out = bb3(f,C + Ccycle,x,alpha,t,UB);%t,info,G,source,dest,UB,sol,@(xx)plotSol(xx,px,py,pz,row,col,'b'));
		assign(x,out.x);
		assign(alpha,out.alpha);
		assign(t,out.t);
		assign(f,out.bestUB);
	end
	
	% Plot solution path
	plotSol(x,px,py,pz,row,col,'g');
	
	% Find cycles in solution
	H = sparse(row,col,double(x),size(G,1),size(G,2)); % solution graph
	[~,CC] = graphconncomp(H,'Weak',true);
	
	cset = setdiff(CC(round(double(alpha))>0),CC(source)); % do not include path or solitary nodes
	
	if isempty(cset), break; end
	for cc=cset
		Csum = 0; Csumr = 0;
		for j=find(CC==cc) % accumulate edges in cycle
			Csum = Csum + x(edgelinearindex(G,j,find(H(j,:))));
			Csumr = Csumr + x(edgelinearindex(G,find(H(j,:)),j)); % include reverse cycle as well
		end
		Ccycle = [Ccycle, Csum <= nnz(CC==cc)-1, Csumr <= nnz(CC==cc)-1]; %#ok
	end
	
end

x = double(x);
alpha = double(alpha);
lb = double(f);

path = walkPath(G,round(x),source,dest);

% if type>=3, disp([double(u) double(t) double(abs(-1/t - u))]); end
end

function plotSol(x,px,py,pz,row,col,color)
% 	sfigure(3);
	path = reshape([row(double(x)~=0)'; col(double(x)~=0)'],[],1);
	for kk=1:2:length(path)-1
		plot3(px(path(kk:kk+1)),py(path(kk:kk+1)),pz(path(kk:kk+1)),[color '+-'])
		hold on
	end
	hold off
	axis equal
	drawnow

end
