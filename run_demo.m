% Discrete optimal view path planning
%
% This file sets up and solves some of the problems in the paper
%	Haner and Heyden, "Discrete optimal view path planning", VISAPP 2015
% 
% It requires YALMIP (http://users.isy.liu.se/johanl/yalmip/) and while
% e.g. SeDuMi will work (sort of) a faster convex solver such as MOSEK
% (http://www.mosek.com) is highly recommended.

%% Basic:
X = bsxfun(@plus,randn(3,10)*1,[-15 -2 8]');
p = generateProblemFromImage(ones(15,15,3,'uint8')*255,...
	'X',X,'poi',mean(X,2),'neigh_rad',1.5,'fov',90,'zmax',inf);
p.source = 218;
p.dest = 8;
p.costfun.F = 'maxeig';
p.costfun.L = 30;
p.costfun.lambda = [];

plotProblem(p,'edges',0); camva(25)
disp('Press any key to solve planning SDP...'); pause

p = solvePlanningSDP(p);
plotProblem(p,'sdpsol',1); camva(25)
disp('Press any key to run reduced GA...'); pause

p = runReducedGA(p,'nodes',40);
plotProblem(p,'rpath',1,'edges',0); camva(25)
disp('Press any key to run GA...'); pause

p = runGA(p);
plotProblem(p,'gapath',1,'edges',0); camva(25)

fprintf('Objective value gap: %f (%.2f%%)\n',p.gasol.obj-p.sdpsol.obj,...
	(p.gasol.obj-p.sdpsol.obj)/p.sdpsol.obj*100);



%% Basic 2:
X = bsxfun(@plus,randn(3,10)*2,[-5 -2 8]');
p = generateProblemFromImage(ones(15,25,3,'uint8')*255,...
	'X',X,'poi',mean(X,2),'neigh_rad',1.5,'fov',90,'zmax',inf);
p.source = 368;
p.dest = 8;
p.costfun.F = 'maxeig';
p.costfun.L = 30;
p.costfun.lambda = [];

plotProblem(p,'edges',0); camva(25)
disp('Press any key to solve planning SDP...'); pause

p = solvePlanningSDP(p);
plotProblem(p,'sdpsol',1); camva(25)
disp('Press any key to run reduced GA...'); pause

p = runReducedGA(p,'nodes',40);
plotProblem(p,'rpath',1,'edges',0); camva(25)
disp('Press any key to run GA...'); pause

p = runGA(p);
plotProblem(p,'gapath',1,'edges',0); camva(25)

fprintf('Objective value gap: %f (%.2f%%)\n',p.gasol.obj-p.sdpsol.obj,...
	(p.gasol.obj-p.sdpsol.obj)/p.sdpsol.obj*100);



%% DAG
X = bsxfun(@plus,randn(3,10)*2,[-5 -2 8]');
p = generateProblemFromImage(ones(15,15,3,'uint8')*255,...
	'X',X,'poi',mean(X,2),'neigh_rad',1.5,'fov',90,'zmax',inf);
p.source = 218;
p.dest = 8;
p.costfun.F = 'maxeig';
p.costfun.L = 40;
p.costfun.lambda = [];

p.G = graph2dag(p.G,p.dest,'toward');

plotProblem(p,'edges',0); camva(25)
disp('Press any key to solve LP...'); pause

p = solveLinearApprox(p);
plotProblem(p,'lppath',1,'edges',0); camva(25)
disp('Press any key to solve planning SDP optimally with branch-and-bound...'); pause

p = solvePlanningSDP(p,'bnb',1);
plotProblem(p,'sdppath',1,'edges',0); camva(25)

fprintf('Objective value gap: %f (%.2f%%)\n',p.lpsol.obj-p.sdpsol.obj,...
	(p.lpsol.obj-p.sdpsol.obj)/p.sdpsol.obj*100);



%% Obstacles:
p = generateProblemFromImage(imread('obstacle.png'),...
	'poi',[10.5 -1 10.5; 4.5 -1 4.5]',...
	'neigh_rad',1.5,...
	'fov',90,...
	'zmax',7);

p.costfun.F = 'trace';
p.costfun.L = 22;
p.costfun.lambda = [];

plotProblem(p);
disp('Press any key to solve planning SDP...'); pause

p = solvePlanningSDP(p);
plotProblem(p,'sdpsol',1);
disp('Press any key to solve linear approximation...'); pause

p = solveLinearApprox(p);
plotProblem(p,'lppath',1,'edges',0);
disp('Press any key to run GA seeded with LP solution...'); pause

p = runGA(p);
plotProblem(p,'gapath',1,'edges',0);

fprintf('Objective value gap: %f (%.2f%%)\n',p.gasol.obj-p.sdpsol.obj,...
	(p.gasol.obj-p.sdpsol.obj)/p.sdpsol.obj*100);



%% Explore a room
p = generateProblemFromImage(imread('room7x7.png'),...
	'rotations',[0 pi 0; 0 pi*3/2 0; 0 pi/2 0; 0 0 0]',...
	'maxangle',180,...
	'neigh_rad',1.5,...
	'fov',90,...
	'zmax',5);

p.costfun.F = 'maxeig';
p.costfun.L = 20;
p.costfun.lambda = [];

plotProblem(p);
disp('Press any key to solve planning SDP...'); pause

p = solvePlanningSDP(p);
plotProblem(p,'sdpsol',1);
disp('Press any key to run GA...'); pause

p = runGA(p);
plotProblem(p,'gapath',1,'edges',0);

fprintf('Objective value gap: %f (%.2f%%)\n',p.gasol.obj-p.sdpsol.obj,...
	(p.gasol.obj-p.sdpsol.obj)/p.sdpsol.obj*100);


%% Mixed
p = generateOmnidirectional(imread('mixed.png'),...
	'maxangle',90,'neigh_rad',1.5,'fov',90,'zmax',15);

p.costfun.F = 'trace';
p.costfun.L = 43;
p.costfun.lambda = [];

plotProblem(p)
disp('Press any key to solve planning SDP...'); pause

p = solvePlanningSDP(p);
plotProblem(p,'sdpsol',1);
disp('Press any key to run reduced GA...'); pause

p = runReducedGA(p,'nodes',30);
plotProblem(p,'rpath',1,'edges',0);
disp('Press any key to run GA...'); pause

p = runGA(p);
plotProblem(p,'gapath',1,'edges',0);

fprintf('Objective value gap: %f (%.2f%%)\n',p.gasol.obj-p.sdpsol.obj,...
	(p.gasol.obj-p.sdpsol.obj)/p.sdpsol.obj*100);



%% Monument problem
load('monument_aligned');

[x y] = meshgrid(linspace(-8,8,20),linspace(-8,8,20));
d = sqrt(x.^2+y.^2);
im = ones(size(x));
im(d<2 | d>8) = 0;
im = repmat(im,[1 1 3])*255;

Xt = bsxfun(@plus,X,[10.5 0 10.5]');


p = generateProblemFromImage(im,...
	'X',Xt(:,center),'poi',mean(Xt,2),'neigh_rad',1.5,'fov',60,'zmax',inf,'free_ends',[0 1]);
p.X = Xt;
p.source = 70;
% p.dest = 71; % uncomment to fix end point destination
p.costfun.F = 'maxeig';
p.costfun.L = 30;
p.costfun.lambda = [];

plotProblem(p,'edges',0)
disp('Press any key to solve planning SDP...'); pause

p = solvePlanningSDP(p);
plotProblem(p,'sdpsol',1);
disp('Press any key to run GA algorithm on reduced graph...'); pause

p = runReducedGA(p,'nodes',40);
disp('Press any key to refine solution with GA algorithm on full graph...'); pause

p = runGA(p);
plotProblem(p,'gapath',1,'edges',0);

fprintf('Objective value gap: %f (%.2f%%)\n',p.gasol.obj-p.sdpsol.obj,...
	(p.gasol.obj-p.sdpsol.obj)/p.sdpsol.obj*100);


%% Castle problem
load('orebro_aligned');

[x y] = meshgrid(linspace(-20,20,50),linspace(-20,20,50));
d = sqrt(x.^2+y.^2);
im = ones(size(x));
im(d<5 | d>12) = 0;
im = repmat(im,[1 1 3])*255;

Xt = bsxfun(@plus,X*0.4,[27 0 25.5]');


p = generateProblemFromImage(im,...
	'X',Xt(:,center),'poi',mean(Xt,2),'neigh_rad',1.5,'fov',40,'zmax',7,'free_ends',[0 0]);
p.X = Xt;
p.source = 70;
p.dest = 421;
% p.dest = 69;
p.costfun.F = 'trace';
p.costfun.L = 50;
p.costfun.lambda = [];

plotProblem(p,'edges',0); camva(30)
disp('Press any key to solve planning SDP...'); pause

p = solvePlanningSDP(p);
plotProblem(p,'sdpsol',1); camva(30)
disp('Press any key to run GA algorithm on reduced graph...'); pause

p = runReducedGA(p,'nodes',40);
disp('Press any key to refine solution with GA algorithm on full graph...'); pause

p = runGA(p);
plotProblem(p,'gapath',1,'edges',0); camva(30)

fprintf('Objective value gap: %f (%.2f%%)\n',p.gasol.obj-p.sdpsol.obj,...
	(p.gasol.obj-p.sdpsol.obj)/p.sdpsol.obj*100);

