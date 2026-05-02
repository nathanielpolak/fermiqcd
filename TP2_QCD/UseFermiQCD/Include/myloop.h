
/// Compute various average (in space-time sens) loop on plane mu-nu:
/// this is well defined even if one link has zero length

/// averaged wilson loops are real by definition
/// so one can take the real part before averaging
//
// this one is averaged over all sites of the lattice
mdp_real my_average_loop(gauge_field &U,int mu,int lmu,int nu,int lnu) {
  int nc=U.nc,l;
  mdp_matrix A1(nc,nc),A2(nc,nc);
  mdp_real tmp=0,loop_at_x;
  site x(U.lattice()),loopx(U.lattice());
  forallsites(x)
 {
  A1=A2=mdp_identity(nc);
  loopx=x;
  for(l=0;l<lmu;l++) {A1=A1*U(loopx,mu);loopx=loopx+mu;} //links along mu
  for(l=0;l<lnu;l++) {A1=A1*U(loopx,nu);loopx=loopx+nu;} //links along nu
  loopx=x;
  for(l=0;l<lnu;l++) {A2=A2*U(loopx,nu);loopx=loopx+nu;} //links along nu
  for(l=0;l<lmu;l++) {A2=A2*U(loopx,mu);loopx=loopx+mu;} //links along mu
  loop_at_x = real(trace(A1*hermitian(A2)));
  tmp+=loop_at_x;
  mdp.add(tmp);}
  return tmp/(U.lattice().nvol_gl*U.nc);
}
// same but with time origin  fixed
mdp_real my_average_loop(gauge_field &U,int mu,int lmu,int nu,int lnu,int x0) {
  int nc=U.nc,l;
  mdp_matrix A1(nc,nc),A2(nc,nc);
  mdp_real tmp=0,loop_at_x;
  site x(U.lattice()),loopx(U.lattice());
  mdp_real space_vol=U.lattice().nx[1]*U.lattice().nx[2]*U.lattice().nx[3];
  forallsites(x)
 {if(x(0)==x0%U.lattice().nx[0])
 {
  A1=A2=mdp_identity(nc);
  loopx=x;
  for(l=0;l<lmu;l++) {A1=A1*U(loopx,mu);loopx=loopx+mu;} //links along mu
  for(l=0;l<lnu;l++) {A1=A1*U(loopx,nu);loopx=loopx+nu;} //links along nu
  loopx=x;
  for(l=0;l<lnu;l++) {A2=A2*U(loopx,nu);loopx=loopx+nu;} //links along nu
  for(l=0;l<lmu;l++) {A2=A2*U(loopx,mu);loopx=loopx+mu;} //links along mu
  loop_at_x = real(trace(A1*hermitian(A2)));
  tmp+=loop_at_x;
  mdp.add(tmp);
  }
  }
  return tmp/(space_vol*U.nc);
}
// this is special. only for tests
// 2D wilson loop with one corner fixed  Inelegant: one should pass the coordinate as a vector
mdp_real my_local_loop(gauge_field &U,int mu,int lmu,int nu,int lnu,int x0,int x1)
 {
  int nc=U.nc,l;
  mdp_matrix A1(nc,nc),A2(nc,nc);
  mdp_real tmp=0,loop_at_x;
  site x(U.lattice()),loopx(U.lattice());
  x.set(x0,x1);
  A1=A2=mdp_identity(nc);
  loopx=x;
  for(l=0;l<lmu;l++) {A1=A1*U(loopx,mu);loopx=loopx+mu;} //links along mu
  for(l=0;l<lnu;l++) {A1=A1*U(loopx,nu);loopx=loopx+nu;} //links along nu
  loopx=x;
  for(l=0;l<lnu;l++) {A2=A2*U(loopx,nu);loopx=loopx+nu;} //links along nu
  for(l=0;l<lmu;l++) {A2=A2*U(loopx,mu);loopx=loopx+mu;} //links along mu
  loop_at_x = real(trace(A1*hermitian(A2)));
  tmp+=loop_at_x;
  mdp.add(tmp);
  return tmp/(U.nc);
}
// wilson loop with one corner fixed
mdp_real my_local_loop(gauge_field &U,int mu,int lmu,int nu,int lnu,int x0,int x1,int x2,int x3)
 {
  int nc=U.nc,l;
  mdp_matrix A1(nc,nc),A2(nc,nc);
  mdp_real tmp=0,loop_at_x;
  site x(U.lattice()),loopx(U.lattice());
  x.set(x0,x1,x2,x3);
  A1=A2=mdp_identity(nc);
  loopx=x;
  for(l=0;l<lmu;l++) {A1=A1*U(loopx,mu);loopx=loopx+mu;} //links along mu
  for(l=0;l<lnu;l++) {A1=A1*U(loopx,nu);loopx=loopx+nu;} //links along nu
  loopx=x;
  for(l=0;l<lnu;l++) {A2=A2*U(loopx,nu);loopx=loopx+nu;} //links along nu
  for(l=0;l<lmu;l++) {A2=A2*U(loopx,mu);loopx=loopx+mu;} //links along mu
  loop_at_x = real(trace(A1*hermitian(A2)));
  tmp+=loop_at_x;
  mdp.add(tmp);
  return tmp/(U.nc);
}
// polyakov loop is a priori complex. This one is averaged over space position
mdp_complex my_polyakov_loop(gauge_field &U) {
  int nc=U.nc,l;
  int mu=0,lmu=U.lattice().nx[0];
  mdp_matrix A1(nc,nc);
  mdp_complex tmp=0.,loop_at_x;
  site x(U.lattice()),loopx(U.lattice());
  mdp_real space_vol=U.lattice().nx[1]*U.lattice().nx[2]*U.lattice().nx[3];
  forallsites(x)
 {if(x(0)==0)
 {
  A1=mdp_identity(nc);
  loopx=x;
  for(l=0;l<lmu;l++) {A1=A1*U(loopx,mu);loopx=loopx+mu;} //links along time
  loop_at_x = trace(A1);
  tmp+=loop_at_x;
  mdp.add(tmp);
  }
  }
  return tmp/(space_vol*U.nc);
}
//This one has a fixed space position
mdp_complex my_polyakov_loop(gauge_field &U,int x1,int x2,int x3) {
  int nc=U.nc,l;
  int mu=0,lmu=U.lattice().nx[0];
  mdp_matrix A1(nc,nc);
  mdp_complex tmp=0.,loop_at_x;
  site x(U.lattice()),loopx(U.lattice());
  x.set(0,x1,x2,x3);
//
  A1=mdp_identity(nc);
  loopx=x;
  for(l=0;l<lmu;l++) {A1=A1*U(loopx,mu);loopx=loopx+mu;} //links along time
  loop_at_x = trace(A1);
  tmp+=loop_at_x;
  mdp.add(tmp);
//
  return tmp/(U.nc);
}

