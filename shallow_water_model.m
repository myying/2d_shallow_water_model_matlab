clear all; close all;

colormap(textread('BlueRed.rgb')/255);

%parameters
dx=100000; dy=dx;   %grid spacing (m)
nx=101; ny=101;     %number of grid points
dt=100;             %time step (s)
nt=1000;            %number of time steps

cx=0; cy=0;     %advection velocity
hbar=8000;      %base level of height
amp=100;        %hill amplitude
hcx=51; hcy=51; %center location of hill
hw=20;          %hill width

% physics parameters
g=9.8;      %gravity
f=0.0001;   %Coriolis
ka=100;     %diffusion coef

nifcor=1;   %if 1, include Coriolis
nifwind=1;  %if 1, include other terms besides Coriolis in u,v eqns
nifdif=1;   %if 1, include diffusion
nifad=1;    %if 1, include advection
balanced_wind=1; %if 1, set initial wind in geostrophic balance with height.

numerics=2; % 1 = FT CS
            % 2 = CT CS
            % 4 = RK4 in time CS
bc_type=1;  % 0 = cyclic
            % 1 = reflective; use with initial u,v=0
       
show_animation=1;      %if 1, show animation; turn this off when profiling
animation_stride=10;   %if =n, skip every n frames
animation_delay=0.001; %delay in between frame in seconds
show_wind=1;           %if 1, plot wind vectors on top of height

% initial condition
x=2:nx+1; y=2:ny+1;
time=0; %time in seconds
[j,i]=meshgrid(1:ny+2,1:nx+2);
h(1:ny+2,1:nx+2)=hbar+amp*exp(-((j'-hcy).^2+(i'-hcx).^2)/((hw/4)^2));
u(1:ny+2,1:nx+2)=cx;
v(1:ny+2,1:nx+2)=cy;
if(balanced_wind==1) %geostrophic wind (ug,vg)
  u(y,x)=-(g/f)*(h(y+1,x)-h(y-1,x))/(2*dy);
  v(y,x)= (g/f)*(h(y,x+1)-h(y,x-1))/(2*dx);
end

% model integration in time
t_start=cputime;
for n=1:nt
  if(numerics==1) %FTCS (Euler)
    ua=u; va=v; ha=h;
    ua(y,x)=u(y,x)+dt*rhs_u_cs(u,v,h,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    va(y,x)=v(y,x)+dt*rhs_v_cs(u,v,h,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    ha(y,x)=h(y,x)+dt*rhs_h_cs(u,v,h,nx,ny,dx,dy,nifad);
    u=ua; v=va; h=ha;
 
  elseif(numerics==2) %CTCS (leapfrog)
    if(n==1) %first time step Euler forward
      u_1=u; v_1=v; h_1=h;
      ua=u; va=v; ha=h;
      ua(y,x)=u(y,x)+dt*rhs_u_cs(u,v,h,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
      va(y,x)=v(y,x)+dt*rhs_v_cs(u,v,h,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
      ha(y,x)=h(y,x)+dt*rhs_h_cs(u,v,h,nx,ny,dx,dy,nifad);
      u_2=ua; v_2=va; h_2=ha;
    end
    u_2(y,x)=u_1(y,x)+2*dt*rhs_u_cs(u,v,h,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    v_2(y,x)=v_1(y,x)+2*dt*rhs_v_cs(u,v,h,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    h_2(y,x)=h_1(y,x)+2*dt*rhs_h_cs(u,v,h,nx,ny,dx,dy,nifad);
    u_1=u; v_1=v; h_1=h;
    u=u_2; v=v_2; h=h_2;
  
  elseif(numerics==4) %RK4 CS
    urhs1=rhs_u_cs(u,v,h,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    vrhs1=rhs_v_cs(u,v,h,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    hrhs1=rhs_h_cs(u,v,h,nx,ny,dx,dy,nifad);
    u1=u; v1=v; h1=h;
    u1(y,x)=u(y,x)+0.5*dt*urhs1;
    v1(y,x)=v(y,x)+0.5*dt*vrhs1;
    h1(y,x)=h(y,x)+0.5*dt*hrhs1;
    urhs2=rhs_u_cs(u1,v1,h1,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    vrhs2=rhs_v_cs(u1,v1,h1,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    hrhs2=rhs_h_cs(u1,v1,h1,nx,ny,dx,dy,nifad);
    u2=u; v2=v; h2=h;    
    u2(y,x)=u(y,x)+0.5*dt*urhs2;
    v2(y,x)=v(y,x)+0.5*dt*vrhs2;
    h2(y,x)=h(y,x)+0.5*dt*hrhs2;
    urhs3=rhs_u_cs(u2,v2,h2,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    vrhs3=rhs_v_cs(u2,v2,h2,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    hrhs3=rhs_h_cs(u2,v2,h2,nx,ny,dx,dy,nifad);
    u3=u; v3=v; h3=h;
    u3(y,x)=u(y,x)+dt*urhs3;
    v3(y,x)=v(y,x)+dt*vrhs3;
    h3(y,x)=h(y,x)+dt*hrhs3;
    urhs4=rhs_u_cs(u3,v3,h3,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    vrhs4=rhs_v_cs(u3,v3,h3,nx,ny,dx,dy,g,f,nifcor,nifwind,nifad);
    hrhs4=rhs_h_cs(u3,v3,h3,nx,ny,dx,dy,nifad);
    u(y,x)=u(y,x)+dt*(urhs1/6+urhs2/3+urhs3/3+urhs4/6);
    v(y,x)=v(y,x)+dt*(vrhs1/6+vrhs2/3+vrhs3/3+vrhs4/6);
    h(y,x)=h(y,x)+dt*(hrhs1/6+hrhs2/3+hrhs3/3+hrhs4/6);
  end
  
  %add diffusion
  if(nifdif==1)
    u=diffusion_ftcs(u,nx,ny,ka,dx,dy,dt);
    v=diffusion_ftcs(v,nx,ny,ka,dx,dy,dt);
    h=diffusion_ftcs(h,nx,ny,ka,dx,dy,dt);
  end

  %boundary condition
  if(bc_type==0) %cyclic
    u(1,:)=u(end-1,:); u(:,1)=u(:,end-1); u(end,:)=u(2,:); u(:,end)=u(:,2);
    v(1,:)=v(end-1,:); v(:,1)=v(:,end-1); v(end,:)=v(2,:); v(:,end)=v(:,2);
    h(1,:)=h(end-1,:); h(:,1)=h(:,end-1); h(end,:)=h(2,:); h(:,end)=h(:,2);
    
  elseif(bc_type==1) %reflective
    u(:,1)=0; u(:,end)=0;
    v(1,:)=0; v(end,:)=0;
    h(1,:)=h(2,:); h(end,:)=h(end-1,:); h(:,1)=h(:,2); h(:,end)=h(:,end-1); 
  end  

  time=time+dt;

  % plot results
  if(show_animation==1 && mod(time,animation_stride*dt)==0)
    surf(h(y,x),'linestyle','none'); view(-30,60);
    caxis(hbar+([-amp amp]/5)); %colorbar;
    hold on;
    if(show_wind==1)
      quiver3(h(y,x),u(y,x),v(y,x),zeros(ny,nx),'k');
    end
    axis([2 nx+1 2 ny+1 hbar-amp hbar+amp]);
    xlabel('x'); ylabel('y'); zlabel('z');
    hold off
    title(['t = ' num2str(time) ' s'])
    pause(animation_delay)
  end

  if(sum(any(h>1e10))>0)
    disp(['model exploded at t=' num2str(time) 's (' num2str(n) ' steps)']); break; 
  end
end

t_end=cputime;
disp(['model integration takes ' num2str(t_end-t_start) ' s']);
