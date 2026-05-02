	tfinal=time(0);
         if(systemchoice==0)sprintf(gracetitle,"title \"CPU: %d secs, k=%d<%d\"",tfinal-tinitial,k,niter); if(systemchoice==1)sprintf(gracetitle,"title \"GPU: %d secs, k=%d<%d\"",tfinal-tinitial,k,niter);
         if(systemchoice==2)sprintf(gracetitle,"title \"CPP: %d secs, k=%d<%d\"",tfinal-tinitial,k,niter);
         GracePrintf("focus g2");
         GracePrintf(gracetitle);
         GracePrintf ("redraw");