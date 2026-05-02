// make fit for this z  to  get value of potential
// skip fit for this z if number of good points in t is too small (< npar=3)
//
	par[0]=1.;par[1]=1.;par[2]=1.;
	mdp_matrix covar(npar,npar);
	if(data_index+1>=npar)
	  {BaesyanLevenbergMarquardt(ftvalue,ydy, 0,data_index, par, npar, covar, fitlog);}
	potential[zloop]=par[0];
	poterror[zloop]=sqrt(abs(covar(0,0)));
	ofstream tmpfit("tmpfit.dat");
	int fitpointnumber=-1;
	for(float t=ftvalue[0];t<=ftvalue[tmax-1]+1.;t=t+0.1)
	    {fitpointnumber++;
	    tmpfit<<t<<" "<<par[0]<<" "<<par[0]+par[1]*exp(-par[2]*t)<<endl;
	    sprintf(replace_tmpfit[fitpointnumber],"g%d.s2 point %f , %f \n",zloop,t,par[0]+par[1]*exp(-par[2]*t));
	    sprintf(replace_tmpfit_as[fitpointnumber],"g%d.s1 point %f , %f \n",zloop,t,par[0]);
	           }
	tmpfit.flush();tmpfit.close();