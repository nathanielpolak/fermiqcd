
/// Compute various average (in space-time sens) loop on plane mu-nu:  
/// this is well defined even if one link has zero length
mdp_complex my_complex_loop(gauge_field &U,int mu,int lmu,int nu,int lnu) {
  int nc=U.nc,l;
  mdp_matrix A1(nc,nc),A2(nc,nc);
  mdp_complex tmp=0.,loop_at_x;
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
  loop_at_x = trace(A1*hermitian(A2));
  tmp+=loop_at_x;
  mdp.add(tmp);}
  return tmp/(U.lattice().nvol_gl*U.nc);
} 
//
mdp_complex my_complex_loop(gauge_field &U,int mu,int lmu,int nu,int lnu,int x0,int x1,int x2,int x3) {
  int nc=U.nc,l;
  mdp_matrix A1(nc,nc),A2(nc,nc);
  mdp_complex tmp=0.,loop_at_x;
  site x(U.lattice()),loopx(U.lattice());
  x.set(x0,x1,x2,x3);
  A1=A2=mdp_identity(nc);
  loopx=x;
  for(l=0;l<lmu;l++) {A1=A1*U(loopx,mu);loopx=loopx+mu;} //links along mu
  for(l=0;l<lnu;l++) {A1=A1*U(loopx,nu);loopx=loopx+nu;} //links along nu
  loopx=x;
  for(l=0;l<lnu;l++) {A2=A2*U(loopx,nu);loopx=loopx+nu;} //links along nu
  for(l=0;l<lmu;l++) {A2=A2*U(loopx,mu);loopx=loopx+mu;} //links along mu
  loop_at_x = trace(A1*hermitian(A2));
  tmp+=loop_at_x;
  mdp.add(tmp);
  return tmp/(U.lattice().nvol_gl*U.nc);
} 






// averaged wilson loops are real, so one can take the real part before averaging
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



//  estimate overlap, after timne evolution "timelong",
//  between the straigth link along direction 1 of length xlong
//  and a link with same ends but deformed in the 2,3 directions
mdp_real shift_loop(gauge_field &U,int timelong,int xlong) {
	int nc=U.nc,l;
	int mu=0,lmu=timelong; 
	int nu=1,lnu=xlong;
	mdp_matrix A1(nc,nc),A2(nc,nc);
	mdp_real tmp=0,loop_at_x;
	site x(U.lattice()),loopx(U.lattice());
	forallsites(x)
	{
		A1=A2=mdp_identity(nc);
		loopx=x;
		for(l=0;l<lmu;l++) {A1=A1*U(loopx,mu);loopx=loopx+mu;}  //links along mu=0 unchanged
		A1=A1*U(loopx,2);loopx=loopx+2;                         //shift in direction 2
		for(l=0;l<lnu;l++) {A1=A1*U(loopx,nu);loopx=loopx+nu;}
		loopx=x;
		for(l=0;l<lnu;l++) {A2=A2*U(loopx,nu);loopx=loopx+nu;} //links along nu at t=0,unchamged
		for(l=0;l<lmu;l++) {A2=A2*U(loopx,mu);loopx=loopx+mu;} //links along mu unchanged
		A2=A2*U(loopx,2);                                      //shift
		loop_at_x = real(trace(A1*hermitian(A2)));
		tmp+=loop_at_x;
		mdp.add(tmp);}
		return tmp/(U.lattice().nvol_gl*U.nc);
} 
//
// a loop with straigth links along time (=timelong)
// at the initial (i) and final (f) times we superpose links along 1
 // shifted by 1 in the 2,3 amd -2,-3 direction with some amplitude
 // in the hope that this has better overlap with the physical string
mdp_real my_fat_loop(gauge_field &U,int timelong,int xlong,mdp_real mixing=0.) 
{       if(mixing==0) return my_average_loop(U,0,timelong,1,xlong);
	int nc=U.nc,l;
	int mu=0,lmu=timelong; 
	int nu=1,lnu=xlong;
	mdp_matrix timelink0(nc,nc),timelink1(nc,nc);
	mdp_matrix f_link(nc,nc),f_linkup2(nc,nc),f_linkdown2(nc,nc),f_linkup3(nc,nc),f_linkdown3(nc,nc),f_fatlink(nc,nc);
	mdp_matrix i_link(nc,nc),i_linkup2(nc,nc),i_linkdown2(nc,nc),i_linkup3(nc,nc),i_linkdown3(nc,nc),i_fatlink(nc,nc);
	mdp_real tmp=0,loop_at_x;
	site x(U.lattice()),loopx(U.lattice()),xplust(U.lattice()),xplusx(U.lattice()); 
	
	forallsites(x)
	{xplust=x;for(l=0;l<lmu;l++) {xplust=xplust+mu;} //define corner x+lmu*0
	 xplusx=x;for(l=0;l<lnu;l++) {xplusx=xplusx+nu;} //define corner x+lnu*1
	//
	//   time links
	timelink0=mdp_identity(nc);
	loopx=x;
	for(l=0;l<lmu;l++) {timelink0=timelink0*U(loopx,mu);loopx=loopx+mu;}  //link along t at x
	//
	timelink1=mdp_identity(nc);
	loopx=xplusx;
	for(l=0;l<lmu;l++) {timelink1=timelink1*U(loopx,mu);loopx=loopx+mu;}  
	timelink1=hermitian(timelink1);                                         //link along t at x+lnu*1
	//
	//   final space links
	loopx=xplust;
	f_link=mdp_identity(nc);
	for(l=0;l<lnu;l++){f_link=f_link*U(loopx,1);loopx=loopx+nu;}           //link along 1 at final time in plane
	//
	f_linkup2=U(xplust,2);loopx=xplust+2;
	for(l=0;l<lnu;l++){f_linkup2=f_linkup2*U(loopx,1);loopx=loopx+nu;}                             
	f_linkup2=f_linkup2*U(loopx,-1,2);                                    //link along 1 at final time shifted up 2
	//
	f_linkdown2=U(xplust,-1,2);loopx=xplust-2;
	for(l=0;l<lnu;l++){f_linkdown2=f_linkdown2*U(loopx,1);loopx=loopx+nu;}                             
	f_linkdown2=f_linkdown2*U(loopx,2);                                    //link along 1 at final time shifted down 2
	//
	f_linkup3=U(xplust,3);loopx=xplust+3;
	for(l=0;l<lnu;l++){f_linkup3=f_linkup3*U(loopx,1);loopx=loopx+nu;}                             
	f_linkup3=f_linkup3*U(loopx,-1,3);                                    //link along 1 at final time shifted up 3
	//
	f_linkdown3=U(xplust,-1,3);loopx=xplust-3;
	for(l=0;l<lnu;l++){f_linkdown3=f_linkdown3*U(loopx,1);loopx=loopx+nu;}                             
	f_linkdown3=f_linkdown3*U(loopx,3);                                    //link along 1 at final time shifted down 3
	//
	// linear superposition
	f_fatlink=f_link+mixing*(f_linkup2+f_linkup3+f_linkdown2+f_linkdown3);
	//
	//   initial space links
	loopx=x;
	i_link=mdp_identity(nc);
	for(l=0;l<lnu;l++){i_link=i_link*U(loopx,1);loopx=loopx+nu;}           //link along 1 at initial time in plane	
	//
	i_linkup2=U(x,2);loopx=x+2;
	for(l=0;l<lnu;l++){i_linkup2=i_linkup2*U(loopx,1);loopx=loopx+nu;}                             
	i_linkup2=i_linkup2*U(loopx,-1,2);                                    //link along 1 at final time shifted up 2
	//
	i_linkdown2=U(x,-1,2);loopx=x-2;
	for(l=0;l<lnu;l++){i_linkdown2=i_linkdown2*U(loopx,1);loopx=loopx+nu;}                             
	i_linkdown2=i_linkdown2*U(loopx,2);                                    //link along 1 at final time shifted down 2
	//
	i_linkup3=U(x,3);loopx=x+3;
	for(l=0;l<lnu;l++){i_linkup3=i_linkup3*U(loopx,1);loopx=loopx+nu;}                             
	i_linkup3=i_linkup3*U(loopx,-1,3);                                    //link along 1 at final time shifted up 2
	//
	i_linkdown3=U(x,-1,3);loopx=x-3;
	for(l=0;l<lnu;l++){i_linkdown3=i_linkdown3*U(loopx,1);loopx=loopx+nu;}                             
	i_linkdown3=i_linkdown3*U(loopx,3);                                    //link along 1 at final time shifted down 2
	//
	// linear superposition
	i_fatlink=hermitian(  i_link+mixing*(i_linkup2+i_linkup3+i_linkdown2+i_linkdown3)   );
	//
	//  average over sites
	loop_at_x=real(trace(timelink0*f_fatlink*timelink1*i_fatlink));
	tmp+=loop_at_x;
	mdp.add(tmp);}
	return tmp/(U.lattice().nvol_gl*U.nc);
} 

