//PLOT CPU TIME
	  ofstream cputime("cputime.dat",ios::app);
	  cputime<<k<<" "<<t_update<<" "<<t_loop<<" "<<endl;
	  cputime.close();
	  GracePrintf ("focus g6");
	  if(k>kupdate)
	   {GracePrintf("kill s0");GracePrintf("kill s1");}
	  GracePrintf ("s0 symbol 0");
	  GracePrintf ("s0 linestyle 1");
	  GracePrintf ("s0 line color 2");
	  GracePrintf ("s0 legend \"update\"");
	  GracePrintf ("s1 symbol 0");
	  GracePrintf ("s1 linestyle 1");
	  GracePrintf ("s1 line color 4");
	  GracePrintf ("s1 legend \"loop\"");
	  GracePrintf("autoscale onread xyaxes");
	  GracePrintf ("read nxy \"cputime.dat\" ");
	  GracePrintf ("redraw");
	  