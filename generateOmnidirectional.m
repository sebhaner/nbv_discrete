function p1 = generateOmnidirectional(im,varargin)


p1 = generateProblemFromImage(im,'rotations',[0 0 0]',varargin{:});
p2 = generateProblemFromImage(im,'rotations',[0 pi/2 0]',varargin{:});
p3 = generateProblemFromImage(im,'rotations',[0 pi 0]',varargin{:});
p4 = generateProblemFromImage(im,'rotations',[0 3*pi/2 0]',varargin{:});

% Add information matrices
for i=2:size(p1.info)
	for j=1:size(p1.info{1},3)
		p1.info{i}(:,:,j) = p1.info{i}(:,:,j) + p2.info{i}(:,:,j)+ p3.info{i}(:,:,j) + p4.info{i}(:,:,j);
	end
end

p1.Pshow = [p1.Pshow p2.Pshow p3.Pshow p4.Pshow];