///////////////////////
// THINGS FOR FITTING //
///////////////////////
  int npar=3;
  float par[npar];
  float potential[zmax],poterror[zmax];            
  Measure ydy[tmax];                //         ydy[i].set(pot[i],err[i]);  
  int nparpot=3;
  float parpot[nparpot],parpoterror[nparpot];
  Measure potdpot[zmax];
  float string_tension=(450./197.33)*(450./197.33);
  float a_latt,da_latt;
  