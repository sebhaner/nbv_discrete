function hc = plotCamera(P,fast,model,num,fc,scale,showtext,ec)
if nargin<8, ec = 'k'; end
if nargin<7, showtext = 0; end
if nargin<6, scale = 0.05; end
if nargin<5, fc = []; end
if nargin<4, num = []; end
if nargin<3, model = 0; end
if nargin<2, fast = 0; end

if fast, ec = fc; fc = 'none'; end

if model==0
	vertices = [-1 -1 1; 1 -1 1; 1 1 1; -1 1 1; -1 1 -1; 1 1 -1; 1 -1 -1;
				-1 -1 -1; -1 1 -3; 1 1 -3; 1 -1 -3; -1 -1 -3; 0 0 0; -0.5 1 1;
				0.5 1 1; 0.5 1 -1; -0.5 1 -1; -0.5 2 0; 0.5 2 0]*scale;

	faces = [1 2 3 4; 2 3 6 7; 3 4 5 6; 1 2 7 8;
			 1 4 5 8; 5 6 7 8; 9 10 13 13; 10 11 13 13;
			 11 12 13 13; 12 9 13 13; 14 15 19 18; 16 17 18 19;
			 15 16 19 19; 14 17 18 18];
elseif model==1
	vertices = [-1 1 -3; 1 1 -3; 1 -1 -3; -1 -1 -3; 0 0 0; 0 1.5 0]*scale;
	faces = [1 2 3 4; 1 2 5 5; 3 4 5 5; 1 4 5 5; 2 3 5 5; 5 5 6 6];
% 	faces = [1 2 3; 1 3 4; 1 2 5; 3 4 5; 1 4 5; 2 3 5; 5 5 6]; % Asymptote only handles triangles (now fixed in exporter)
elseif model==2
	vertices = [-1 1 -3; 1 1 -3; 1 -1 -3; -1 -1 -3; 0 0 0; 0 1.5 0; 0 0 -3]*scale;
	faces = [1 2 3 4; 5 5 6 6; 7 7 5 5];
else
	vertices = [0 0 0; 0 1.5 0; 0 0 -3]*scale;
	faces = [1 1 2 2; 1 1 3 3];
end

flip = [1 0 0; 0 -1 0; 0 0 -1];
Pt = [P(:,1:3)' -P(:,1:3)'*P(:,4)];
tvert = Pt*[flip*vertices'; ones(1,size(vertices,1))];

color = [repmat([0 0.2 0.8],8,1); repmat([0.4 1 0.6],5,1); repmat([1 0.6 0],6,1)];
color = color(1:size(vertices,1),:);

if isempty(fc)
	hc = patch('Vertices',tvert(1:3,:)','Faces',faces,'FaceVertexCData',color);
	shading faceted;
else
	hc = patch('Vertices',tvert(1:3,:)','Faces',faces,'FaceColor',fc,'EdgeColor',ec);
end

if showtext
	tp = Pt*[0 -2.5*scale 0 1]';
	text(tp(1),tp(2),tp(3),num2str(num),'Background',[1 1 1],'FontSize',8);
end
