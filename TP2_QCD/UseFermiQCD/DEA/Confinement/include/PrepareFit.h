//prepare fit, do not use bad points
	      if(flagnumber[tloop]==0)
	          {data_index++;
	          ydy[data_index].set(averaged_loop[zloop][tloop] , jack_error[zloop][tloop]);
	          ftvalue[data_index]=tvalue[tloop];
	          tmp<<tvalue[tloop]<<" "<<averaged_loop[zloop][tloop]<<" "<<jack_error[zloop][tloop]<<endl;
/// it does not hurt to specify the graph number but it was not necessary because the focus is set on the goo graph at the time of piping
	sprintf(replace_tmp[data_index],"g%d.s0 point %d , %f \n",zloop,tvalue[tloop],averaged_loop[zloop][tloop]);
	sprintf(replace_tmp_x[data_index],"g%d.s0.x[%d]=%d\n",zloop,data_index,tvalue[tloop]);
	sprintf(replace_tmp_y[data_index],"g%d.s0.y0[%d]=%f\n",zloop,data_index,averaged_loop[zloop][tloop]);
	sprintf(replace_tmp_dy[data_index],"g%d.s0.y1[%d]=%f\n",zloop,data_index,jack_error[zloop][tloop]);
	           }
		logfile<<"iteration: "<<k<<" zloop: "<<zloop<<" tloop: "<<tloop<<" jack result: "<<averaged_loop[zloop][tloop]<<" "<<jack_error[zloop][tloop]<<"\n";