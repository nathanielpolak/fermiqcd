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
//PLOT  POTENTIAL
	  GracePrintf ("focus g8");
	  if(k>kupdate){GracePrintf("kill s0");GracePrintf("kill s1");}
	  GracePrintf ("s0 symbol 1");GracePrintf ("s1 symbol 0");
          GracePrintf ("s0 symbol size 0.4");
          GracePrintf ("s0 linestyle 0");GracePrintf ("s1 linestyle 1");
          GracePrintf ("s1 line color 4");
          GracePrintf ("s0 legend \"a*V(z/a)\" ");GracePrintf ("s1 legend \"A0+A1*z+A2/z\" ");
	  GracePrintf ("read xydy \"tmppot.dat\" ");GracePrintf ("read xy \"tmpfitpot.dat\" ");
	  GracePrintf (gracezmin);
	  GracePrintf (gracezmax);
	  GracePrintf ("autoscale yaxes");
          GracePrintf ("redraw");
          