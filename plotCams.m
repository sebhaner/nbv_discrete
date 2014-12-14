function varargout = plotCams(P,X,varargin)

cmodel = 1;
pstyle = '.';
pcolor = 'b';
ccolor = 'k';
cecolor = [1 1 1]*0.0;
bgcolor = 'w';
cmap = [];
cscale = 'auto';
msize = 1;
wf = 1;
shownum = 0;

i = 1;
while i<=length(varargin)
	switch varargin{i}
		case 'cmodel', cmodel = varargin{i+1}; i=i+2;
		case 'pstyle', pstyle = varargin{i+1}; i=i+2;
		case 'ccolor', ccolor = varargin{i+1}; i=i+2;
		case 'cecolor', cecolor = varargin{i+1}; i=i+2;
		case 'cscale', cscale = varargin{i+1}; i=i+2;
		case 'wf', wf = varargin{i+1}; i=i+2;
		case 'shownum', shownum = varargin{i+1}; i=i+2;
		case 'pcolor', pcolor = varargin{i+1}; i=i+2;
		case 'bgcolor', bgcolor = varargin{i+1}; i=i+2;
		case 'colormap', cmap = varargin{i+1}; i=i+2;
		case 'msize', msize = varargin{i+1}; i=i+2;
		otherwise, error('Parameter parse error');
	end
end

figexist = strcmp(get(gcf,'Tag'),'3dview');
ax = gca;
cuv = get(ax,'CameraUpVector');
cpos = get(ax,'CameraPosition');
ctrg = get(ax,'CameraTarget');
va = get(ax,'CameraViewAngle');
pmode = get(ax,'Projection');
hp = [];

if(~isempty(X))
	if ischar(pcolor)
		hp = plot3(X(1,:),X(2,:),X(3,:),[pstyle pcolor],'MarkerSize',msize);
	elseif numel(pcolor)==3
		hp = plot3(X(1,:),X(2,:),X(3,:),pstyle,'Color',pcolor,'MarkerSize',msize);
	elseif isempty(cmap)
		hp = plot3c(X(1:3,:),pcolor,32,'MarkerSize',msize,'Marker',pstyle);
	else
		hp = plot3cm(X(1:3,:),pcolor,cmap,'MarkerSize',msize,'Marker',pstyle);
	end
	M = median(X(1:3,:),2);
% 	if strcmpi(cscale,'auto')
		% make sure outliers don't impact
% 		small = sum(X.^2)<prctile(sum(X.^2),99);
% 		cscale = sqrt(max(eig(cov(X(:,small)')))) * 0.1;
% 	end
	% 	V = [median(abs(X(1,:)-M(1))); median(abs(X(2,:)-M(2))); median(abs(X(3,:)-M(3)))];
end
hc = [];
if ~isempty(P)
	ind = find(~cellfun(@isempty,P))';
	if ~isempty(ind)
		if strcmpi(cscale,'auto')
			t = zeros(3,length(ind));
			for i=1:length(ind)
				t(:,i) = -P{ind(i)}(:,1:3)'*P{ind(i)}(:,4);
			end
			t = mean(sqrt(sum((t(:,2:end)-t(:,1:end-1)).^2)));
			if t==0, t = 0.1; end
			cscale = t*0.3;
		end
		
		if length(ccolor)>3 && ~ischar(ccolor)
			[cind_t cmap] = getColorIndex(ccolor(ind),32);
			cind = zeros(1,length(P)); cind(ind) = cind_t;
			if cmodel<4
				for c=ind
					hc(end+1) = plotCamera(P{c},wf,cmodel,c,cmap(cind(c),:),cscale,shownum,cecolor); %#ok
				end
			else
				cc = zeros(3,length(ind));
				for c=ind
					cc(:,c) = -P{c}(:,1:3)'*P{c}(:,4);
				end
				hold on; plot3cm(cc,cind,cmap,'Marker','*'); hold off;
			end
			
		else
			
			if cmodel<4
				for c=ind
					hc(end+1) = plotCamera(P{c},wf,cmodel,c,ccolor,cscale,shownum,cecolor); %#ok
				end
			else
				cc = zeros(3,length(ind));
				for c=ind
					cc(:,c) = -P{c}(:,1:3)'*P{c}(:,4);
				end
				hold on; plot3(cc(1,:),cc(2,:),cc(3,:),['*' ccolor]); hold off;
			end
			
		end
	end
end

if figexist
	set(ax,'CameraViewAngleMode','manual','CameraViewAngle',va);
	set(ax,'CameraUpVector',cuv,'CameraPosition',cpos,'CameraTarget',ctrg,'Projection',pmode);
	set(ax,'XColor',[1 0 0],'YColor',[0 1 0],'ZColor',[0 0 1]);
else
	if ~isempty(X) && all(isfinite(M)), ct = M; else ct = [0.5; 0.5; 0.5]; end
	set(ax,'CameraTarget',ct(1:3));
	set(gcf,'Tag','3dview','Toolbar','none');
	cameratoolbar(gcf,'Show');
	cameratoolbar(gcf,'SetMode','orbit');
	cameratoolbar(gcf,'SetCoordSys','y');
end
set(gcf,'Color',bgcolor);
% set(gcf,'Renderer','zbuffer');
set(ax,'Color',bgcolor);
axis equal;
axis off;

if nargout>0, varargout{1} = hp; varargout{2} = hc; end

% if ~isempty(X), axis([M(1)-V(1) M(1)+V(1) M(2)-V(2) M(2)+V(2) M(3)-V(3) M(3)+V(3)]); end
% axis(axis + [-1 1 -1 1 -1 1]*10);

