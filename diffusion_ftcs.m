% diffusion - forward in time centered in space (FTCS) scheme
function var=diffusion_ftcs(var,nx,ny,ka,dx,dy,dt)
	x=2:nx+1; y=2:ny+1;
	var(y,x) = var(y,x) + ...
	dt*ka*( (var(y+1,x)+var(y-1,x)-2*var(y,x))/(dy^2) ...
	       +(var(y,x+1)+var(y,x-1)-2*var(y,x))/(dx^2) );
end