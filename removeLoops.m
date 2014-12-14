function p = removeLoops(p)

succ = zeros(1,max(p));
for i=1:length(p)-1
	succ(p(i)) = p(i+1);
end
for i=2:length(p)
	p(i) = succ(p(i-1));
	if p(i)==p(end), break; end
end
p = p(1:i);
