% diffusion - alternating direction implicit (ADI) scheme
%    - need tridiagonal.m
%
function var=diffusion_adi(var,nx,ny,ka,dx,dy,dt)
	x=2:nx+1; y=2:ny+1;
	var1=var; var2=var;
	ax=eye(nx+2);
	for i=x
		ax(i,i-1:i+1)=[-ka/(dx^2) 2*ka/(dx^2)+2/dt -ka/(dx^2)];
	end
	for j=y
		qx=var(j,1:nx+2);
		qx(x)=var(j,x)*(2/dt)+ka*(var(j,x+1)-2*var(j,x)+var(j,x-1))/(dx^2);
		var1(j,1:nx+2)=tridiagonal(ax,qx,nx+2);
	end
	var1(1,:)=var(1,:); var1(ny+2,:)=var(ny+2,:);
	ay=eye(ny+2);
	for j=y
		ay(j,j-1:j+1)=[-ka/(dy^2) 2*ka/(dy^2)+2/dt -ka/(dy^2)];
	end
	for i=x
		qy=var1(1:nx+2,i);
		qy(y)=var1(y,i)*(2/dt)+ka*(var1(y+1,i)-2*var1(y,i)+var1(y-1,i))/(dy^2);
		var2(1:ny+2,i)=tridiagonal(ay,qy,ny+2);
	end
	var2(:,1)=var1(:,1); var2(:,nx+2)=var1(:,nx+2);
	var=var2;
end
