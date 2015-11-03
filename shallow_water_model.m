clear all;

% grid and numeric parameters
nx=37;
ny=37;
dt=100;
dx=100000; dy=dx;
nt=1000;

% physics parameters
g=9.8;
f=0; %0.0001;
ka=0; %0.0001*dx^2;

numerics=4; % 0 = forward in time
            % 1 = centered in time
            % 4 = Runge-Kutta 4th order in time
bc_type=2;  % 0 = cyclic
            % 1 = reflective
            % 2 = radiation
            
show_animation=1;
animation_stride=1;
animation_delay=0.001;            

% initial condition
time=0;
hbar=8000;
amp=100;
[j,i]=meshgrid(1:ny+2,1:nx+2);
%h(1:ny+2,1:nx+2)=hbar+amp*peaks(nx+2)/10;  %perturb h
h(1:ny+2,1:nx+2)=hbar+amp*exp(-((j'-20).^2+(i'-20).^2)/(5^2));
u(1:ny+2,1:nx+2)=0; %no wind
v(1:ny+2,1:nx+2)=0;
%u(y,x)=-(g/f)*(h(y+1,x)-h(y-1,x))/(2*dy); %geostrophic wind (ug,vg)
%v(y,x)= (g/f)*(h(y,x+1)-h(y,x-1))/(2*dx);

% model integration in time
for n=1:nt

  x=2:nx+1; y=2:ny+1;
  u0=u; v0=v; h0=h; %save last time step
  
  if(numerics==0) %FT (Euler)
    u(y,x)=u(y,x)+dt*rhs_u(u,v,h,nx,ny,dx,dy,g,f);
    v(y,x)=v(y,x)+dt*rhs_v(u,v,h,nx,ny,dx,dy,g,f);
    h(y,x)=h(y,x)+dt*rhs_h(u,v,h,nx,ny,dx,dy);
 
  elseif(numerics==1) %CT (leapfrog)
    if(n==1) %first time step Euler forward
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
  
  elseif(numerics==4) %RK4
    urhs1=rhs_u(u,v,h,nx,ny,dx,dy,g,f);
    vrhs1=rhs_v(u,v,h,nx,ny,dx,dy,g,f);
    hrhs1=rhs_h(u,v,h,nx,ny,dx,dy);
    u1=u; v1=v; h1=h;
    u1(y,x)=u(y,x)+0.5*dt*urhs1;
    v1(y,x)=v(y,x)+0.5*dt*vrhs1;
    h1(y,x)=h(y,x)+0.5*dt*hrhs1;
    urhs2=rhs_u(u1,v1,h1,nx,ny,dx,dy,g,f);
    vrhs2=rhs_v(u1,v1,h1,nx,ny,dx,dy,g,f);
    hrhs2=rhs_h(u1,v1,h1,nx,ny,dx,dy);
    u2=u; v2=v; h2=h;
    u2(y,x)=u(y,x)+0.5*dt*urhs2;
    v2(y,x)=v(y,x)+0.5*dt*vrhs2;
    h2(y,x)=h(y,x)+0.5*dt*hrhs2;
    urhs3=rhs_u(u2,v2,h2,nx,ny,dx,dy,g,f);
    vrhs3=rhs_v(u2,v2,h2,nx,ny,dx,dy,g,f);
    hrhs3=rhs_h(u2,v2,h2,nx,ny,dx,dy);
    u3=u; v3=v; h3=h;
    u3(y,x)=u(y,x)+dt*urhs3;
    v3(y,x)=v(y,x)+dt*vrhs3;
    h3(y,x)=h(y,x)+dt*hrhs3;
    urhs4=rhs_u(u3,v3,h3,nx,ny,dx,dy,g,f);
    vrhs4=rhs_v(u3,v3,h3,nx,ny,dx,dy,g,f);
    hrhs4=rhs_h(u3,v3,h3,nx,ny,dx,dy);
    u(y,x)=u(y,x)+dt*(urhs1/6+urhs2/3+urhs3/3+urhs4/6);
    v(y,x)=v(y,x)+dt*(vrhs1/6+vrhs2/3+vrhs3/3+vrhs4/6);
    h(y,x)=h(y,x)+dt*(hrhs1/6+hrhs2/3+hrhs3/3+hrhs4/6);
  end
  
  %add diffusion
  u=diffusion_ftcs(u,nx,ny,ka,dx,dy,dt);
  v=diffusion_ftcs(v,nx,ny,ka,dx,dy,dt);
  h=diffusion_ftcs(h,nx,ny,ka,dx,dy,dt);

  %boundary condition
  if(bc_type==0) %cyclic
    u(1,:)=u(end-1,:); v(1,:)=v(end-1,:); h(1,:)=h(end-1,:);
    u(:,1)=u(:,end-1); v(:,1)=v(:,end-1); h(:,1)=h(:,end-1);
    u(end,:)=u(2,:); v(end,:)=v(2,:); h(end,:)=h(2,:);
    u(:,end)=u(:,2); v(:,end)=v(:,2); h(:,end)=h(:,2);
  elseif(bc_type==1) %reflective
  end  

  time=time+dt;

  % plot results
  if(show_animation==1 && mod(time,animation_stride*dt)==0)
    surf(h(y,x),'linestyle','none'); view(-30,30);
    caxis(hbar+([-amp amp]/5)); %colorbar;
    hold on;
    quiver3(h(y,x),u(y,x),v(y,x),zeros(ny,nx),'k');
    axis([2 nx+1 2 ny+1 hbar-amp hbar+amp]);
    xlabel('x'); ylabel('y');
    hold off
    title(['t = ' num2str(time) ' s'])
    pause(animation_delay)
  end

  if(sum(any(h>1e10))>0)
    disp(['model exploded at t=' num2str(time) 's (' num2str(n) ' steps)']); break; 
  end
end
