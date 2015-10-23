%  S H A L L O W   W A T E R   M O D E L
%         by Yue Ying   2015
%
%  State variables:  u, v, h
%  Dimensions: nx,ny; the actual state is stored in (2:nx+1, 2:ny+1)
%    with 1 layer of fictitious variables to help calculate gradients
%    along the boundary.
%
%  Model equations:
%    dh/dt = -u*dh/dx -v*dh/dy -h*(du/dx+dv/dy) +ka*(d2/dx2+d2/dy2)h
%    du/dt = -u*du/dx -v*du/dy -g*dh/dx +f*v +ka*(d2/dx2+d2/dy2)u
%    dv/dt = -u*dv/dx -v*dv/dy -g*dh/dy -f*u +ka*(d2/dx2+d2/dy2)v
%
%  Rossby deformation radius l_R = sqrt(g*h_base)/f
%    if initial h perturbation has scale larger than l_R, the initial (u,v)=(ug,vg)
%    will remain in balance with the height field. Otherwise, u,v,h will undergo 
%    geostrophic adjustment.
%

clear all; close all;

% grid and numeric parameters
nx=37;
ny=37;
dt=100;
dx=100000; dy=dx;
nt=1000;

% physics parameters
g=9.8;
f=0.0001;
ka=0.0001*dx^2;

% initial condition
x=2:nx+1; y=2:ny+1;
h_base=8000;
p_amp=100;
[j,i]=meshgrid(1:ny+2,1:nx+2);
h(1:ny+2,1:nx+2)=h_base+p_amp*peaks(nx+2)/10;  %perturb h
%h(1:ny+2,1:nx+2)=h_base+p_amp*exp(-((j'-20).^2+(i'-20).^2)/(5^2));
u(1:ny+2,1:nx+2)=0; %no wind
v(1:ny+2,1:nx+2)=0;
%u(y,x)=-(g/f)*(h(y+1,x)-h(y-1,x))/(2*dy); %geostrophic wind (ug,vg)
%v(y,x)= (g/f)*(h(y,x+1)-h(y,x-1))/(2*dx);

% model integration in time
for t=1:nt
	%Euler forward
	%u(y,x)=u(y,x)+dt*rhs_u(u,v,h,nx,ny,dx,dy,g,f);
	%v(y,x)=v(y,x)+dt*rhs_v(u,v,h,nx,ny,dx,dy,g,f);
	%h(y,x)=h(y,x)+dt*rhs_h(u,v,h,nx,ny,dx,dy);

    %Centered in time (Leap frog scheme)
    if(t==1) %first time step Euler forward
        u_1=u; v_1=v; h_1=h;
        u(y,x)=u(y,x)+dt*rhs_u(u,v,h,nx,ny,dx,dy,g,f);
        v(y,x)=v(y,x)+dt*rhs_v(u,v,h,nx,ny,dx,dy,g,f);
        h(y,x)=h(y,x)+dt*rhs_h(u,v,h,nx,ny,dx,dy);
        u_2=u; v_2=v; h_2=h;
    end
    u_2(y,x)=u_1(y,x)+2*dt*rhs_u(u,v,h,nx,ny,dx,dy,g,f);
    v_2(y,x)=v_1(y,x)+2*dt*rhs_v(u,v,h,nx,ny,dx,dy,g,f);
    h_2(y,x)=h_1(y,x)+2*dt*rhs_h(u,v,h,nx,ny,dx,dy);
    u_1=u; v_1=v; h_1=h;
    u=u_2; v=v_2; h=h_2;
    
    
	%RK4
	%urhs1=rhs_u(u,v,h,nx,ny,dx,dy,g,f);
	%vrhs1=rhs_v(u,v,h,nx,ny,dx,dy,g,f);
	%hrhs1=rhs_h(u,v,h,nx,ny,dx,dy);
	%u1=u; v1=v; h1=h;
	%u1(y,x)=u(y,x)+0.5*dt*urhs1;
	%v1(y,x)=v(y,x)+0.5*dt*vrhs1;
	%h1(y,x)=h(y,x)+0.5*dt*hrhs1;
	%urhs2=rhs_u(u1,v1,h1,nx,ny,dx,dy,g,f);
	%vrhs2=rhs_v(u1,v1,h1,nx,ny,dx,dy,g,f);
	%hrhs2=rhs_h(u1,v1,h1,nx,ny,dx,dy);
	%u2=u; v2=v; h2=h;
	%u2(y,x)=u(y,x)+0.5*dt*urhs2;
	%v2(y,x)=v(y,x)+0.5*dt*vrhs2;
	%h2(y,x)=h(y,x)+0.5*dt*hrhs2;
	%urhs3=rhs_u(u2,v2,h2,nx,ny,dx,dy,g,f);
	%vrhs3=rhs_v(u2,v2,h2,nx,ny,dx,dy,g,f);
	%hrhs3=rhs_h(u2,v2,h2,nx,ny,dx,dy);
	%u3=u; v3=v; h3=h;
	%u3(y,x)=u(y,x)+dt*urhs3;
	%v3(y,x)=v(y,x)+dt*vrhs3;
	%h3(y,x)=h(y,x)+dt*hrhs3;
	%urhs4=rhs_u(u3,v3,h3,nx,ny,dx,dy,g,f);
	%vrhs4=rhs_v(u3,v3,h3,nx,ny,dx,dy,g,f);
	%hrhs4=rhs_h(u3,v3,h3,nx,ny,dx,dy);
	%u(y,x)=u(y,x)+dt*(urhs1/6+urhs2/3+urhs3/3+urhs4/6);
	%v(y,x)=v(y,x)+dt*(vrhs1/6+vrhs2/3+vrhs3/3+vrhs4/6);
	%h(y,x)=h(y,x)+dt*(hrhs1/6+hrhs2/3+hrhs3/3+hrhs4/6);
    
 	%add diffusion
	u=diffusion_ftcs(u,nx,ny,ka,dx,dy,dt);
	v=diffusion_ftcs(v,nx,ny,ka,dx,dy,dt);
	h=diffusion_ftcs(h,nx,ny,ka,dx,dy,dt);

	%boundary condition: cyclic: a(1)=a(end-1); outflow: a(1)=a(2);
	u(1,:)=u(end-1,:); 
	v(1,:)=v(end-1,:);
	h(1,:)=h(end-1,:);
	u(:,1)=u(:,end-1);
	v(:,1)=v(:,end-1);
	h(:,1)=h(:,end-1);
	u(end,:)=u(2,:);
	v(end,:)=v(2,:);
	h(end,:)=h(2,:);
	u(:,end)=u(:,2);
	v(:,end)=v(:,2);
	h(:,end)=h(:,2);

	% plot results
	surf(h(y,x),'linestyle','none'); view(-30,30);
	caxis(h_base+([-p_amp p_amp]/5)); %colorbar;
	hold on;
	quiver3(h(y,x),u(y,x),v(y,x),zeros(ny,nx),'k');
	axis([2 nx+1 2 ny+1 h_base-p_amp h_base+p_amp]);
	xlabel('x'); ylabel('y');
	title(['frame ' num2str(t)]);
	hold off

	pause(0.001)
	    
	if(sum(any(h>1e10))>0)
		disp('model exploded'); break;
	end
end
