//Fit potential according to fitpot
for(int z=0;z<zmax;z++){ potdpot[z].set(potential[z],poterror[z]);}
	   parpot[0]=1.;parpot[1]=1.;parpot[2]=1.;
	   mdp_matrix covarpot(nparpot,nparpot);
	   BaesyanLevenbergMarquardt(fzvalue,potdpot, 0,zmax-1, parpot, nparpot, covarpot, fitpot);
	   for(int n=0;n<nparpot;n++){parpoterror[n]=sqrt(abs(covarpot(n,n)));}
	   ofstream tmppot("tmppot.dat");
	   ofstream tmpfitpot ("tmpfitpot.dat");
	   for(int z=0;z<zmax;z++)
	   {tmppot<<zvalue[z]<<" "<<potential[z]<<" "<<poterror[z]<<endl;
	   }
	   for(float z=0.5;z<=fzvalue[zmax-1]+1.;z=z+0.1)
	   {tmpfitpot<<z<<" "<<parpot[0]+parpot[1]*z+parpot[2]/z<<endl;
	   }
	   tmppot.close();tmpfitpot.close();
