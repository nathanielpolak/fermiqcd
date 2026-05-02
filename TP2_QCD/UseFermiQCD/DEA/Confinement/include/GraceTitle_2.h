    tfinal=time(0);
    if(systemchoice==0)sprintf(gracetitle,"title \"CPU:%dsecs, k=%d<%d\"",tfinal-tinitial,k,nheat);
    if(systemchoice==1)sprintf(gracetitle,"title \"GPU:%dsecs, k=%d<%d\"",tfinal-tinitial,k,nheat);
    if(systemchoice==2)sprintf(gracetitle,"title \"CPP:%dsecs, k=%d<%d\"",tfinal-tinitial,k,nheat);
    sprintf(gracesubtitle,"subtitle \"Heating\"");
    GracePrintf("focus g2");
    GracePrintf (gracetitle);
    GracePrintf (gracesubtitle);
    GracePrintf ("redraw");