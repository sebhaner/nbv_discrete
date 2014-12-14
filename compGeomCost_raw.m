function I = compGeomCost_raw(P,X,fov,zmax)

nop = size(X,2);
w = irodr(P(:,1:3));
I = zeros(3,3,nop);
lim = tan(fov/2);
Rinv = inv(eye(2)*100); % measurement covariance

for p=1:nop

	weight = zcheck(P,X(:,p),lim,zmax);
	
	if weight>0
		J = jacexpPt(w,P(:,4),X(:,p));
		I(:,:,p) = J'*Rinv*J*weight;
	else
		I(:,:,p) = zeros(3);
	end

end


function w = zcheck(P,X,lim,zmax)
x = P*[X; 1];
z = x(3); y = x(2)/z; x = x(1)/z;

w = z>0 && x<lim && x>-lim && y<lim && y>-lim && z<zmax;
