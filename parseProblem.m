function [type param] = parseProblem(p)

if ~xor(isempty(p.costfun.L),isempty(p.costfun.lambda))
	error('Either L or lambda must be empty');
end

if isempty(p.costfun.L)
	type = 1;
	param = p.costfun.lambda;
else
	type = 2;
	param = p.costfun.L;
end
if strcmp(p.costfun.F,'maxeig')
	type = type+2;
elseif ~strcmp(p.costfun.F,'trace')
	error('costfun.F must be either "trace" or "maxeig"');
end
