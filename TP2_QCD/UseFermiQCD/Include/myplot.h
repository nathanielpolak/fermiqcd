// problem here. use mysingleplot

class myplot{

private:
PLINT id1,autoy,acc;
PLINT colbox, collab, colline[4], styline[4];
PLFLT ymin, ymax, xlab, ylab;
PLFLT xmin, xmax, xjump;
char *legline[4], toplab[20];
public:
myplot();
~myplot();
void start();
void init(PLFLT,PLFLT,PLFLT,PLFLT,char[10]);
void plot(PLFLT,PLFLT);
void plot(PLFLT,PLFLT,PLFLT);
void plot(PLFLT,PLFLT,PLFLT,PLFLT);
void plot(PLFLT,PLFLT,PLFLT,PLFLT,PLFLT);
static int  HowManyPlots;
};
// constructor and destructor
myplot::myplot(){
HowManyPlots++;cout<<"create the plot number: "<<HowManyPlots<<endl;}
myplot::~myplot(){
cout<<"destroy the plot: "<<HowManyPlots<<endl;HowManyPlots--;}
void myplot::init(PLFLT xmin,PLFLT xmax,PLFLT ymin,PLFLT ymax,char title[10]){
  PLFLT xjump;
 // xjump=(xmax-xmin)/5.;
  xjump=0.1;
  cout<<"do initialise with default value, "<<title<<endl;
  if(HowManyPlots==1){
 // set output device
	  cout<<"HowManyPlots="<<HowManyPlots<<" I run plsdev"<<endl;
  plsdev("xwin");
  //set options
  plsetopt("db", "");plsetopt("np", "");
//set the size and location of the window
//plspage(1000.,1000.,500,500,100,10);
  // Initialize PLplot.
  plinit();
  pladv(0); 
  plvsta(); }
  // Specify some reasonable defaults for xmin,xmax,xjump,ymin and ymax
  // The plot will grow automatically if needed (but not shrink)
 // xmin=0.;xmax=20.;xjump=1.;
 // ymin = 0.5; ymax = 1.;

  // Axes options same as plbox. Only pen 3 is used
  colbox = 1;
  collab = 3;
  styline[0] = colline[0] = 1;        // pens color and line style
  styline[1] = colline[1] = 2;
  styline[2] = colline[2] = 3;
  styline[3] = colline[3] = 4;

  legline[0] = "";                         // pens legend
  legline[1] = "";
  legline[2] = "";
  legline[3] = "";

  xlab = 0.; ylab = 1.;     // legend position

  autoy = 1;  // autoscale y
  acc = 1;    // don't scrip, accumulate

 
//  define the strip plot and get the id number
  plstripc(&id1, "bcnst", "bcnstv",
             xmin, xmax, xjump, ymin, ymax,
             xlab, ylab,
             autoy, acc,
             colbox, collab,
             colline, styline, legline,
             "x axis", "y axis", title); 
   
  // je sais pas trop a quoi ca set:  Let plplot handle errors from here on

  plsError(NULL, NULL);
  cout<<"plot initialisation done, got id1="<<id1<<endl;
}
void myplot::plot(PLFLT x,PLFLT y0){
plstripa(id1,0,x,y0);pleop();
}
void myplot::plot(PLFLT x,PLFLT y0,PLFLT y1){
plstripa(id1,0,x,y0);plstripa(id1,1,x,y1);pleop();
}
void myplot::plot(PLFLT x,PLFLT y0,PLFLT y1,PLFLT y2){
plstripa(id1,0,x,y0);plstripa(id1,1,x,y1);plstripa(id1,2,x,y2);pleop();
}
void myplot::plot(PLFLT x,PLFLT y0,PLFLT y1,PLFLT y2,PLFLT y3){
plstripa(id1,0,x,y0);plstripa(id1,1,x,y1);plstripa(id1,2,x,y2);plstripa(id1,3,x,y3);pleop();
}
int myplot::HowManyPlots=0;
