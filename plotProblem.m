function plotProblem(p,varargin)

plotsdpsol = 0;
plotgapath = 0;
plotsdppath = 0;
plotlppath = 0;
plotrhcpath = 0;
plotrpath = 0;
plotrnodes = 0;
plotedges = 1;
plotX = 1;

for i=1:2:length(varargin)
	switch varargin{i}
		case 'sdpsol', plotsdpsol = varargin{i+1};
		case 'gapath', plotgapath = varargin{i+1};
		case 'sdppath', plotsdppath = varargin{i+1};
		case 'lppath', plotlppath = varargin{i+1};
		case 'rhcpath', plotrhcpath = varargin{i+1};
		case 'rpath', plotrpath = varargin{i+1};
		case 'rnodes', plotrnodes = varargin{i+1};
		case 'edges', plotedges = varargin{i+1};
		case 'points', plotX = varargin{i+1};
		otherwise, error('Parameter parse error');
	end
end


bp = [];
if plotgapath==1 && isfield(p,'gasol') && isfield(p.gasol,'path') && ~isempty(p.gasol.path)
	bp = p.gasol.path;
end
if plotsdppath==1 && isfield(p,'sdpsol') && isfield(p.sdpsol,'path') && ~isempty(p.sdpsol.path)
	bp = p.sdpsol.path;
end
if plotlppath==1 && isfield(p,'lpsol') && isfield(p.lpsol,'path') && ~isempty(p.lpsol.path)
	bp = p.lpsol.path;
end
if plotrhcpath==1 && isfield(p,'rhcsol') && isfield(p.rhcsol,'path') && ~isempty(p.rhcsol.path)
	bp = p.rhcsol.path;
end
if plotrpath==1 && isfield(p,'rsol') && isfield(p.rsol,'path') && ~isempty(p.rsol.path)
	bp = p.rsol.path;
end

if ~plotX, p.X = []; cla; end
hcam = [];
cscale = 0.15; if size(p.Pshow,2)>1, cscale = 0.1; end
for i=1:size(p.Pshow,2)
	P = p.Pshow(:,i);
	if i==1, pX = p.X; else pX = []; end
	[~,h] = plotCams(P(setdiff(1:length(P),[bp(:)' p.source p.dest])),pX,'cmodel',1,'wf',0,...
		'ccolor',[128 128 255]/255,'cscale',cscale,'msize',5,'pcolor','r');
	plotCams(P([p.source p.dest]),[],'cmodel',1,'wf',0,'ccolor',[255 160 0]/255,'cscale',cscale);
	hcam = [hcam; h(:)]; %#ok
end


hold on
color = [193 221 198]/255;
if plotgapath==2 && isfield(p,'gasol') && isfield(p.gasol,'path') && ~isempty(p.gasol.path)
	plot3(p.x(p.gasol.path),p.y(p.gasol.path),p.z(p.gasol.path),'linewidth',5,'color',color);
end
if plotsdppath==2 && isfield(p,'sdpsol') && isfield(p.sdpsol,'path') && ~isempty(p.sdpsol.path)
	plot3(p.x(p.sdpsol.path),p.y(p.sdpsol.path),p.z(p.sdpsol.path),'linewidth',5,'color',color);
end
if plotlppath==2 && isfield(p,'lpsol') && isfield(p.lpsol,'path') && ~isempty(p.lpsol.path)
	plot3(p.x(p.lpsol.path),p.y(p.lpsol.path),p.z(p.lpsol.path),'linewidth',5,'color',color);
end
if plotrhcpath==2 && isfield(p,'rhcsol') && isfield(p.rhcsol,'path') && ~isempty(p.rhcsol.path)
	plot3(p.x(p.rhcsol.path),p.y(p.rhcsol.path),p.z(p.rhcsol.path),'linewidth',5,'color',color);
end
if plotrpath==2 && isfield(p,'rsol') && isfield(p.rsol,'path') && ~isempty(p.rsol.path)
	plot3(p.x(p.rsol.path),p.y(p.rsol.path),p.z(p.rsol.path),'linewidth',5,'color',color);
end


if ~isempty(bp)
	for i=1:size(p.Pshow,2)
		P = p.Pshow(:,i);
		plotCams(P(bp(2:end-1)),[],'cmodel',1,'wf',0,'ccolor',[0 255 0]/255,'cscale',cscale);
		% plotCams(p.Pshow(bp([1 end])),[],'cmodel',1,'wf',0,'ccolor',[255 160 0]/255,'cscale',0.15);
	end
end


set(hcam,'FaceAlpha',0.2,'EdgeAlpha',0.1);
set(gcf,'renderer','opengl');
plot3(p.x(bp),p.y(bp),p.z(bp),'linewidth',2,'color',[0 128 128]/255)
hold off

axis([-100 100 -10 10 -100 100])
if plotX
	camtarget([(max([p.x p.X(1,:)])+min([p.x p.X(1,:)]))/2 ...
			(max([p.y p.X(2,:)])+min([p.y p.X(2,:)]))/2 ... 
			(max([p.z p.X(3,:)])+min([p.z p.X(3,:)]))/2]);
else
	camtarget([(max(p.x)+min(p.x))/2 ...
			(max(p.y)+min(p.y))/2 ... 
			(max(p.z)+min(p.z))/2]);
end

campos(camtarget + [0 -50 -2])
camup([0 0 -1]);
camva(15);

axis equal

if plotsdpsol && isfield(p,'sdpsol') && ~isempty(p.sdpsol.x)
	plotEdges(p.sdpsol.x,p.G,p.x,p.y,p.z);
elseif plotedges
	plotEdges(ones(nnz(p.G),1)*0.1,p.G,p.x,p.y,p.z);
end

if plotrnodes && isfield(p,'rsol') && ~isempty(p.rsol.beta)
	hold on
	plot3(p.x(p.rsol.beta),p.y(p.rsol.beta),p.z(p.rsol.beta),'ks','markers',6);
	hold off
end


% All this silliness to ensure only the thickest edge is plotted between
% two particular spatial locations, as there may be many nodes in one place
function plotEdges(xx,G,x,y,z) 

[r c] = find(G);
xx = max(0,min(1,xx));

map = containers.Map;

for i=1:numel(xx)
	sn = r(i);
	en = c(i);
	key = hashedge(x(sn),x(en),z(sn),z(en));

	if ~map.isKey(key)
		map(key) = [xx(i) sn en];
	else
		rec = map(key);
		rec(1) = max(rec(1),xx(i));
		map(key) = rec;
	end
end

for k=map.keys
	rec = map(cell2mat(k));
	sn = rec(2);
	en = rec(3);
	xxt = rec(1);
	
% 	ind = find(excpath==sn); if ~isempty(ind) && excpath
	if xxt>0.01
		line([x(sn) x(en)],[y(sn) y(en)],[z(sn) z(en)], ...
			'color',max(0,min(1,[.8 .8 .8]*(1-xxt)+[0.8 0.5 0.8]*xxt)),'linewidth',xxt*5+1);
	end
end

function h = hashedge(x1,x2,y1,y2)

h = sortrows([x1 y1; x2 y2]);
if any(isnan(h))
	h = '';
else
	h = char(h(:)');
end
	
	
	
	
	