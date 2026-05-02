/////////////////////////
// THINGS FOR PLOTTING//
/////////////////////////
  char gracesubtitle[100];
  char gracetitle[100];
  char gracecommand[50],gracetmin[50],gracetmax[50],gracezmin[50],gracezmax[50];
  char gracelegend[20];
  char isanumber1[5],isanumber2[5];
  int flagnumber[tmax];
  // erase lattice_spacing.dat
  ofstream lattice_spacing("lattice_spacing.dat");
  lattice_spacing<<"# lattice spacing"<<endl;
  lattice_spacing.close();
  // erase cputime.dat
  ofstream cputime("cputime.dat");
  cputime<<"#CPU update and loop"<<endl;
  cputime.close();
//
  sprintf(gracetmin,"world xmin %d",tvalue[0]-1);
  sprintf(gracetmax,"world xmax %d",tvalue[tmax-1]+1);
  sprintf(gracezmin,"world xmin %d",zvalue[0]-1);
  sprintf(gracezmax,"world xmax %d",zvalue[zmax-1]+1);
  