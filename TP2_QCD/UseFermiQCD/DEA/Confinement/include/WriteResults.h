///////////////////
// WRITE RESULTS //
///////////////////
int file_index=0;
char filename[20];  //file name is incremented automatically
do{sprintf(filename,"Result_%d",file_index);file_index++;}while(file_exist(filename)==1);
ofstream result_file(filename);
//
result_file<<"Parameters for this calculation:"<<endl;
result_file<<"lattice (t,x,y,z) dimensions: "<<tsize<<" "<<xsize<<" "<<ysize<<" "<<zsize<<endl;
result_file<<"Nb of colors:"<<nc<<",  beta=2*Nc/g^2: "<<mybeta<<endl;
result_file<<"Decorrelation sweeps: "<<decorel<<", Heating sweeps: "<<nheat<<", iterations: "<<niter<<endl;
result_file<<"loop plane: "<<mu<<", "<<nu<<endl;
result_file<<"loops start at t=0 and end at: ";
for (int t=0;t<tmax;t++){result_file<<tvalue[t]<<" ";}
result_file<<endl<<"loops start at z=0 and end at: ";
for (int z=0;z<zmax;z++){result_file<<zvalue[z]<<" ";}
result_file<<endl;
//
for(int z=0;z<zmax;z++)
{result_file<<endl;
result_file<<"t  "<<"        v(t,"<<zvalue[z]<<") "<<endl;             
for (int t=0;t<tmax;t++)
{result_file<<tvalue[t]<<"  "<<averaged_loop[z][t]<<"(+/-)"<<jack_error[z][t]<<endl;}                     
}
result_file<<endl;
result_file<<"z  v(t,z) extrapolated to t=infinity"<<endl;
for (int z=0;z<zmax;z++){result_file<<zvalue[z]<<"  "<<potential[z]<<"(+/-)"<<poterror[z]<<endl;}
result_file<<endl<<"Fit of potential according to v(z)=a[0]+a[1]*z +a[2]/z"<<endl;
for(int n=0;n<npar;n++){result_file<<"a["<<n<<"]= "<<parpot[n]<<"(+/-)"<<parpoterror[n]<<endl;}
result_file<<endl;
result_file<<"lattice spacing in fm: "<<a_latt<<"(+/-)"<<da_latt<<endl;
// 
result_file.close();

