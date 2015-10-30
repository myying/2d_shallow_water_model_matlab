  S H A L L O W   W A T E R   M O D E L
         by Yue Ying   2015

  State variables:  u, v, h
  Dimensions: nx,ny; the actual state is stored in (2:nx+1, 2:ny+1)
    with 1 layer of fictitious variables to help calculate gradients
    along the boundary.

  Model equations:
    dh/dt = -u*dh/dx -v*dh/dy -h*(du/dx+dv/dy) +ka*(d2/dx2+d2/dy2)h
    du/dt = -u*du/dx -v*du/dy -g*dh/dx +f*v +ka*(d2/dx2+d2/dy2)u
    dv/dt = -u*dv/dx -v*dv/dy -g*dh/dy -f*u +ka*(d2/dx2+d2/dy2)v

  Rossby deformation radius l_R = sqrt(g*h_base)/f
    if initial h perturbation has scale larger than l_R, the initial (u,v)=(ug,vg)
    will remain in balance with the height field. Otherwise, u,v,h will undergo 
    geostrophic adjustment.


