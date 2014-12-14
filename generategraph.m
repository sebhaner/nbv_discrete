function [G x y z P start finish Pshow] = generategraph(im,conn_d,noise,rotations,maxanglediff,mode,supernode)

imbw = im2bw(im,0.001);
[i j] = find(imbw);
nb = numel(i);
x = j(:)';
y = zeros(size(x));
z = i(:)';

x = x + randn(size(x))*noise;
y = y + randn(size(y))*noise;
z = z + randn(size(z))*noise;

nrot = size(rotations,2);
x = repmat(x,1,nrot);
y = repmat(y,1,nrot);
z = repmat(z,1,nrot);

n = numel(x);

ny = y;
P = cell(n,1);
Pshow = cell(n,1);
for r=1:nrot
	for c=1:nb
		ci = (r-1)*nb+c;
		if strcmpi(mode,'rot')
			P{ci} = rodr(rotations(:,r))*[eye(3) -[x(ci); y(ci); z(ci)]];
			ny(ci) = y(ci) + (r-1)*0.2*0; % Add vertical space for plotting purposes
			Pshow{ci} = rodr(rotations(:,r))*[eye(3) -[x(ci); ny(ci); z(ci)]];
		elseif strcmpi(mode,'lookat')
			P{ci} = lookat_m([x(ci); y(ci); z(ci)],rotations(:,r))*[eye(3) -[x(ci); y(ci); z(ci)]];
			ny(ci) = y(ci) + (r-1)*0.2*0; % Add vertical space for plotting purposes
			Pshow{ci} = lookat_m([x(ci); y(ci); z(ci)],rotations(:,r))*[eye(3) -[x(ci); ny(ci); z(ci)]];
		end
	end
end

indi = zeros(8*n,1);
indj = zeros(8*n,1);
ival = zeros(8*n,1);
q = 0;
conn_d2 = conn_d^2;
for i=1:n-1
	for j=i+1:n
		d2 = sum([x(i)-x(j) y(i)-y(j) z(i)-z(j)].^2);
		if d2>0 && d2<=conn_d2
			anglediff = norm(irodr(P{i}(:,1:3)'*P{j}(:,1:3)));
			if anglediff/pi*180 <= maxanglediff+1e-4
				% Add edge between node i and j with weight d
				q = q+1;
				indi(q) = i;
				indj(q) = j;
				ival(q) = sqrt(d2);
				% Add reverse edge
				q = q+1;
				indi(q) = j;
				indj(q) = i;
				ival(q) = sqrt(d2);
			end
		end
% 		if d2==0 && norm(irodr(P{i}(:,1:3)'*P{j}(:,1:3)))/pi*180 < 48
% 			q = q+1;
% 			indi(q) = i;
% 			indj(q) = j;
% 			ival(q) = 1.9;
% 			q = q+1;
% 			indi(q) = j;
% 			indj(q) = i;
% 			ival(q) = 1.9;
% 		end
	end
end

G = sparse(indi(1:q),indj(1:q),ival(1:q),n,n);

y = ny;

% Try to find start and finish as green and red pixels
green = im(:,:,2)>im(:,:,1) & im(:,:,2)>im(:,:,3);
red   = im(:,:,1)>im(:,:,2) & im(:,:,1)>im(:,:,3);
ind = find(imbw);
start = find(green(ind),1);
finish = find(red(ind),1);

% Add supernodes
if nargin>=6
	w = eps/10;
	if supernode(1) % create free starting point
		G = [G sparse(size(G,1),1); ones(1,size(G,2))*w 0];
		start = size(G,1);
		P{size(G,1)} = []; Pshow{size(G,1)} = [];
% 		x(end+1) = mean(x); z(end+1) = mean(z); y(end+1) = 5;
		x(end+1) = NaN; z(end+1) = NaN; y(end+1) = NaN;
	end
	if supernode(2) % create free end point
		if supernode(1)
			G = [G [ones(size(G,1)-1,1)*w; 0]; sparse(1,size(G,2)) 0];
		else
			G = [G ones(size(G,1),1)*w; sparse(1,size(G,2)) 0];
		end
		finish = size(G,1);
		P{size(G,1)} = []; Pshow{size(G,1)} = [];
% 		x(end+1) = mean(x); z(end+1) = mean(z); y(end+1) = -5;
		x(end+1) = NaN; z(end+1) = NaN; y(end+1) = NaN;
	end
end

				