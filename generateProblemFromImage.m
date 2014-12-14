function prob = generateProblemFromImage(im,varargin)

structurenoise = 0;
gridnoise = 0;
poi = [];
rot = [];
conn_d = 1.0;
maxanglediff = 180;
supernodes = [0 0];
fov = 90/180*pi;
zmax = 20;
X = [];

for i=1:2:length(varargin)
	switch varargin{i}
		case 'structurenoise', structurenoise = varargin{i+1};
		case 'gridnoise',gridnoise = varargin{i+1};
		case 'poi', poi = varargin{i+1};
		case 'rotations', rot = varargin{i+1};
		case 'neigh_rad', conn_d = varargin{i+1};
		case 'maxangle', maxanglediff = varargin{i+1};
		case 'free_ends', supernodes = varargin{i+1};
		case 'fov', fov = varargin{i+1};
		case 'zmax', zmax = varargin{i+1};
		case 'X', X = varargin{i+1};
		otherwise, error('Parameter parse error');
	end
end

if ~xor(isempty(rot),isempty(poi))
	error('Either "rotations" or "poi" must be specified');
end

if isempty(X), X = map2cloud(im,structurenoise); end

tic
if ~isempty(rot)
	if size(rot,1)~=3, error('"rotations" must be 3-by-n'); end
	[G x y z P source dest Pshow] = generategraph(im,conn_d,gridnoise,rot,maxanglediff,'rot',supernodes);
else
	if size(poi,1)~=3, error('"poi" must be 3-by-n'); end
	[G x y z P source dest Pshow] = generategraph(im,conn_d,gridnoise,poi,maxanglediff,'lookat',supernodes);
end
fprintf('Generate graph: %.2f s\n',toc)

% Generate information matrices
info = cell(size(G,1),1);
tic
info{1} = repmat(eye(3,3),[1,1,size(X,2)])*1e-2; % Initial information I0
for i=1:size(G,1)
	if ~isempty(P{i})
		info{i+1} = compGeomCost_raw(P{i},X,fov/180*pi,zmax);
	else
		info{i+1} = zeros(3,3,size(X,2));
	end
end
fprintf('Generate information matrices: %.2f s\n',toc)

prob.X = X;
prob.G = G;
prob.x = x;
prob.y = y;
prob.z = z;
prob.P = P;
prob.Pshow = Pshow;
prob.source = source;
prob.dest = dest;
prob.info = info;
prob.costfun.F = 'trace';
prob.costfun.L = [];
prob.costfun.lambda = 10;


