//test if a file or directory exist. Does not look inside a directory! Use cd before
//not safe due to timing overlap
int file_or_dir_exist(char *name)
{
//test if a file or directory exist. Does not look inside a directory! Use cd before
//not safe due to timing overlap
char chain[100];string line;int n;
sprintf(chain,"ls | grep ");
strcat(chain,name);
//cout<<"chain is: "<<chain<<endl;
strcat(chain," >tmpfile &");
system(chain);
ifstream tmp("tmpfile");getline(tmp,line);tmp.close();
if( line=="" ) 
{n=0;cout<<name<<" does not exist"<<endl;}
else {n=1;cout<<name<<" already exist"<<endl;}
return n;
}
//this one is safe but works only for a file
int file_exist(char *name)
{
FILE * pFile;
pFile = fopen (name,"r");
if (pFile!=NULL) {fclose (pFile);return 1;}
if (pFile==NULL) {return 0;}
}
