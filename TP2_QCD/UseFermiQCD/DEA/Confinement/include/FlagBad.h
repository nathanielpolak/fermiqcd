	  sprintf(isanumber1,"%0.1f",averaged_loop[zloop][tloop]);
	  sprintf(isanumber2,"%0.1f",jack_error[zloop][tloop]);
	  flagnumber[tloop]=0;
	  if(strcmp(isanumber1,"nan")==0 || strcmp(isanumber2,"nan")==0) {flagnumber[tloop]=1;}