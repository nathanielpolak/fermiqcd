//PLOT log(loop_t/loop_t-1) VERSUS t   AND SHOW FIT FOR THIS Z
	if(zloop<=5)
	{
	  sprintf(gracecommand,"focus g%d",zloop);
	  sprintf(gracelegend,"s0 legend \"z=%d\" ",zvalue[zloop]);
	  GracePrintf (gracecommand);
	  if(k>kupdate){GracePrintf("kill s0");GracePrintf("kill s1");GracePrintf("kill s2");}

	  GracePrintf ("s0 symbol 1");
          GracePrintf ("s0 symbol size 0.4");
          GracePrintf ("s0 symbol color 1");
          GracePrintf ("s0 linestyle 0");
	  GracePrintf("s0 type xydy");
          GracePrintf (gracelegend);
          GracePrintf ("s1 symbol 0");
          GracePrintf ("s1 linestyle 3");
          GracePrintf ("s1 line color 2");
          GracePrintf ("s2 symbol 0");
          GracePrintf ("s2 linestyle 1");
          GracePrintf ("s2 line color 2");
//
/// this way of plotting by reading a file has been replaced by directly  piping each line to solve the timing conflicts in grace.
//	  int errorxydy=GracePrintf ("read xydy \"tmp.dat\" ");
//	  if(errorxydy<=0) {logfile<<"*******  pb in read xydy, error code: "<<errorxydy<<"\n";}
//	  int errornxy=GracePrintf ("read nxy \"tmpfit.dat\" ");
//	  if(errornxy<=0) {logfile<<"******* pb in read nxy, error code: "<<errornxy<<" \n";}
/// this is the new way.
/// this is just to create the set S0
	  for(int index=0;index<=data_index;index++)
	  {
	  GraceCommand(replace_tmp[index]);
	  }
/// now make the real plot
//
	  for(int index=0;index<=data_index;index++)
	  {
	  GraceCommand(replace_tmp_x[index]);
	  GraceCommand(replace_tmp_y[index]);
	  GraceCommand(replace_tmp_dy[index]);
	  }
//
/// plot the fit
	  for (int pt=0;pt<=fitpointnumber;pt++)
	  {
	  GraceCommand(replace_tmpfit[pt]);
	  GraceCommand(replace_tmpfit_as[pt]);
	  }
//
	  GracePrintf (gracetmin);
	  GracePrintf (gracetmax);
	  GracePrintf ("autoscale yaxes");
	  GracePrintf ("redraw");
	 }
