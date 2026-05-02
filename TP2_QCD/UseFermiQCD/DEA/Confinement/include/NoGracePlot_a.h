//SET SCALE WITH THE STRING TENSION and plot a vs iterations
	  a_latt=sqrt(parpot[1]/string_tension);
///strange that the error covarpot[1,1] is not detected
//	  da_latt=0.5*a_latt*sqrt(abs(covarpot[1,1]))/parpot[1];
	  da_latt=0.5*a_latt*parpoterror[1]/parpot[1];
	  sprintf(isanumber1,"%0.1f",a_latt);
	  sprintf(isanumber2,"%0.1f",da_latt);
	  if(strcmp(isanumber1,"nan")!=0 && strcmp(isanumber2,"nan")!=0)
	    {
	      ofstream lattice_spacing("lattice_spacing.dat",ios::app);
	      lattice_spacing<<k<<" "<<a_latt<<" "<<da_latt<<endl;
	      lattice_spacing.close();
	    }

