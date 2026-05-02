
typedef struct {
    float re,im;
} mycomplex;
typedef struct {
   mycomplex m00,m01,m10,m11;} su2;
//
int full_adress(int t,int nt,int x,int nx,int y,int ny,int z,int nz,int i,int j,int nc);
int field_adress(int,int ,int ,int ,int ,int ,int);
int CCfield_adress(int t,int nt,int x,int nx,int i,int j,int nc);
mycomplex CCc_mult(mycomplex,mycomplex);
su2 CCsu2_set(mycomplex ,mycomplex ,mycomplex ,mycomplex );
su2 CCsu2_equal(su2);
su2 CCsu2_adjoint(su2);
su2 CCsu2_prod(su2 ,su2 );
float betterCC_local_loop(mycomplex* Ut,mycomplex* Ux,int nt,int nx,int nc,int tloop,int xloop,int t0,int x0);
int make_aver_loops_with_CPU(float *RfullCCfield_t,float *IfullCCfield_t,float *RfullCCfield_x,float *IfullCCfield_x,int *box,int nc,float* CCstore_aver_loops);
//