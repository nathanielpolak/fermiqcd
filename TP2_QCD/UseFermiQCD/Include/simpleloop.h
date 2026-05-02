
/// Compute various average (in space-time sens) loop on plane mu-nu:  
/// this is well defined even if one link has zero length
//
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
mdp_real my_local_loop(gauge_field &U,int mu,int lmu,int nu,int lnu) {
  int nc=U.nc,l;
  mdp_matrix A1(nc,nc),A2(nc,nc);
  mdp_real tmp=0,loop_at_x;
  site x(U.lattice()),loopx(U.lattice());
  x.set(1,1,1,1);
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
  return tmp/(U.lattice().nvol_gl*U.nc);
} 

