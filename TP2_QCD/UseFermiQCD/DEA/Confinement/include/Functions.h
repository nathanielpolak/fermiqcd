//
// function for jack knife 
//
    float f1(float *x, void *a) 
    {float dt=*((float*)a);
    return log(x[0]/x[1])/dt; 
    }
//
// function for fitting the t dependence of log(W_t/W_t_1), used by Levenberg_Marquard
//
    float fitlog(float x, float *a, long ma, void *junk) {return a[0]+a[1]*exp(-a[2]*x);}
//
// function for fitting the potential
//
    float fitpot(float x, float *a, long ma, void *junk){return a[0]+a[1]*x+a[2]/x;}
//
