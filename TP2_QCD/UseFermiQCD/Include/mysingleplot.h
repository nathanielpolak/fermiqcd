// problem here. use mysingleplot

class mysingleplot{

private:
PLINT test;
PLINT id1,autoy,acc,legtest;
PLINT colbox, collab, colline[4], styline[4];
PLFLT ymin, ymax, xlab, ylab;
PLFLT xmin, xmax, xjump;
char *legline[4], toplab[20];
public:
mysingleplot();
~mysingleplot();
void start();
void setlegend(char[10],char[10],char[10],char[10]);
void init(PLFLT,PLFLT,PLFLT,PLFLT,char[10]);
void plot(PLFLT,PLFLT);
void plot(PLFLT,PLFLT,PLFLT);
void plot(PLFLT,PLFLT,PLFLT,PLFLT);
void plot(PLFLT,PLFLT,PLFLT,PLFLT,PLFLT);
};
// constructor and destructor
mysingleplot::mysingleplot(){
cout<<"create a single plot"<<endl;
legtest=0;
cout<<"legtest="<<legtest<<endl;
}
mysingleplot::~mysingleplot(){
	cout<<"destroy the plot: "<<endl;plend();}

void mysingleplot::setlegend(char leg0[10]="",char leg1[10]="",char leg2[10]="",char leg3[10]=""){
	cout<<"set legend"<<endl;
	legline[0]=leg0;
	legline[1]=leg1;
	legline[2]=leg2;
	legline[3]=leg3;
	legtest=1;
}
	
void mysingleplot::init(PLFLT xmin,PLFLT xmax,PLFLT ymin,PLFLT ymax,char title[10]){
  cout<<"do initialise with default value, "<<title<<endl;	
  PLFLT xjump;
  xjump=0.1;
 // set output device
  cout<<"set device to xwin"<<endl;
  plsdev("xwin");
  //set options
  cout<<"set options"<<endl;
  plsetopt("db", "");plsetopt("np", "");
//set the size and location of the window
//plspage(1000.,1000.,500,500,100,10);
  // Initialize PLplot.
  cout<<"call plinit"<<endl; 
  plinit();
  pladv(0); 
  plvsta(); 
  // Axes options same as plbox.
  colbox = 1;
  collab = 3;
  styline[0] = colline[0] = 1;        // pens color and line style
  styline[1] = colline[1] = 2;
  styline[2] = colline[2] = 3;
  styline[3] = colline[3] = 4;
  cout<<"legende test= "<<legtest<<endl;
  if(legtest==0){
  cout<<"legend has not been defined by setlegend, use empty string"<<endl;
  legline[0] = "";                         // pens legend
  legline[1] = "";
  legline[2] = "";
  legline[3] = "";}

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
	     "Nb of iterations", "a*Pot(r)", title); 
   
  // je sais pas trop a quoi ca set:  Let plplot handle errors from here on

  plsError(NULL, NULL);
  cout<<"plot initialisation done, got id1="<<id1<<endl;
}
void mysingleplot::plot(PLFLT x,PLFLT y0){
plstripa(id1,0,x,y0);pleop();
}
void mysingleplot::plot(PLFLT x,PLFLT y0,PLFLT y1){
plstripa(id1,0,x,y0);plstripa(id1,1,x,y1);pleop();
}
void mysingleplot::plot(PLFLT x,PLFLT y0,PLFLT y1,PLFLT y2){
plstripa(id1,0,x,y0);plstripa(id1,1,x,y1);plstripa(id1,2,x,y2);pleop();
}
void mysingleplot::plot(PLFLT x,PLFLT y0,PLFLT y1,PLFLT y2,PLFLT y3){
plstripa(id1,0,x,y0);plstripa(id1,1,x,y1);plstripa(id1,2,x,y2);plstripa(id1,3,x,y3);pleop();
}
