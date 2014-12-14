function beta = roundSDPsol(p)

alpha = max(0,min(1,p.sdpsol.alpha));
type = parseProblem(p);
pm = inf;


% for i=1:10000
% 	a = alpha > rand(size(alpha));
% 	if nnz(a)==p.rsol.nodes
% % 		s = sum(alpha.*a);
% 		c = evalGeoCost_mex(find(a),p.info,type);
% 		if c<pm, pm = c; am = a; end
% 	end
% end

for i=1:1000
	a = alpha.*rand(size(alpha));
	a([p.source p.dest]) = 2; % make sure they are in the top
	[~,ind] = sort(a,'descend');
	a = ind(1:p.rsol.nodes);
	c = evalGeoCost_mex(a,p.info,type);
	if c<pm, pm = c; am = a; end
end
beta = false(size(alpha));
beta(am) = true;


fprintf('Cost: %.4f\n',pm);



% 	function geo = geocost(path)
% 		
% 		geo = 0;
% 		for pp=1:size(p.info{1},3)
% 			I = p.info{1}(:,:,pp); % I0
% 			for kk=path(:)', I = I + p.info{kk+1}(:,:,pp); end
% 			if type<=2
% 				geo = geo + trace(inv(I));
% 			else
% 				geo = max(geo,1/min(eig(I)));
% 			end
% 		end
% 		
% 	end

end