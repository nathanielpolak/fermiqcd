
void StartGrace(){
    if (GraceOpen(2048) == -1) {
        fprintf (stderr, "Can't run Grace. \n");
        exit (EXIT_FAILURE);
    }
    /* Send some initialization commands to Grace */
//
    GracePrintf ("page size 700 500");
    GracePrintf ("arrange(3,3,0.1,0.,0.)");
//g0
    GracePrintf( "with g0");
    GracePrintf("autoscale onread none");
    GracePrintf("xaxis  tick major 1");
    GracePrintf("xaxis  ticklabel on");
    GracePrintf("xaxis  ticklabel place opposite");
    GracePrintf("xaxis tick place both");
    GracePrintf("xaxis  ticklabel char size 0.750000");
    GracePrintf("yaxis tick place normal");
    GracePrintf("yaxis  ticklabel char size 0.750000");
//set legend for g0
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.33,0.9");
//g1
    GracePrintf( "with g1");
    GracePrintf("autoscale onread none");
    GracePrintf("xaxis  tick major 1");
    GracePrintf("xaxis  ticklabel on");
    GracePrintf("xaxis  ticklabel place opposite");
    GracePrintf("xaxis tick place both");
    GracePrintf("xaxis  ticklabel char size 0.750000");
    GracePrintf("yaxis tick place normal");
    GracePrintf ("xaxis label place opposite");
    GracePrintf ("xaxis label \"t\"");
    GracePrintf("yaxis  ticklabel char size 0.750000");
//set legend for g1
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.73,0.9");
//g2
    GracePrintf( "with g2");
    GracePrintf("autoscale onread none");
    GracePrintf("xaxis  tick major 1");
    GracePrintf("xaxis  ticklabel on");
    GracePrintf("xaxis  ticklabel place opposite");
    GracePrintf("xaxis tick place both");
    GracePrintf("xaxis  ticklabel char size 0.750000");
    GracePrintf("yaxis tick place normal");
    GracePrintf("yaxis  ticklabel char size 0.750000");
//set legend for g2
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 1.13,0.9");
//g3
    GracePrintf( "with g3");
    GracePrintf("autoscale onread none");
    GracePrintf("xaxis  tick major 1");
    GracePrintf("xaxis  ticklabel off");
    GracePrintf("xaxis tick place normal");
    GracePrintf("yaxis tick place normal");
    GracePrintf ("yaxis label \"                [Log W(t,z)-Log W(t',z)]/(t'-t)\"");
    GracePrintf("yaxis  ticklabel char size 0.750000");
//set legend for g3
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.33,0.633");
 //g4
    GracePrintf( "with g4");
    GracePrintf("autoscale onread none");
    GracePrintf("xaxis  tick major 1");
    GracePrintf("xaxis  ticklabel off");
    GracePrintf("xaxis tick place normal");
    GracePrintf("yaxis tick place normal");
    GracePrintf("yaxis  ticklabel char size 0.750000");
//set legend for g4
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.73,0.633");
//g5
    GracePrintf( "with g5");
    GracePrintf("autoscale onread none");
    GracePrintf("xaxis  tick major 1");
    GracePrintf("xaxis  ticklabel off");
    GracePrintf("xaxis tick place normal");
    GracePrintf("yaxis tick place normal");
    GracePrintf("yaxis  ticklabel char size 0.750000");
//set legend for g5
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 1.13,0.633");
//g6
    GracePrintf( "with g6");
    GracePrintf("autoscale onread xyaxes");
    GracePrintf("xaxis tick place normal");
    GracePrintf("yaxis tick place normal");
    GracePrintf ("xaxis label \"iterations\"");
    GracePrintf("xaxis  ticklabel char size 0.750000");
    GracePrintf ("yaxis label \"CPU time\"");
    GracePrintf("yaxis  ticklabel char size 0.750000");
//set legend for g6
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.1,0.366");
//g7
    GracePrintf( "with g7");
    GracePrintf("autoscale onread xyaxes");
    GracePrintf("xaxis tick place normal");
    GracePrintf("yaxis tick place normal");
    GracePrintf ("xaxis label \"iterations\"");
    GracePrintf("xaxis  ticklabel char size 0.750000");
    GracePrintf("yaxis  ticklabel char size 0.750000");
//set legend for g7
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 0.5,0.366");
//g8
    GracePrintf( "with g8");
    GracePrintf("autoscale onread none");
    GracePrintf("xaxis  tick major 1");
    GracePrintf("xaxis tick place normal");
    GracePrintf("xaxis  ticklabel char size 0.750000");
    GracePrintf("yaxis tick place normal");
    GracePrintf ("xaxis label \"z\"");
    GracePrintf("yaxis  ticklabel char size 0.750000");
//set legend for g8
    GracePrintf ("legend on");
    GracePrintf ("legend loctype view");
    GracePrintf ("legend 1.05,0.174");
}
void SaveGrace(){
//   save file.agr and close GRACE
    if (GraceIsOpen()) {
        /* Tell Grace to save the data */
        GracePrintf ("saveall \"Confine.agr\"");
        /* Flush the output buffer and close Grace */
        GraceClose();
        exit (EXIT_SUCCESS);
    } else {
        exit (EXIT_FAILURE);
    }
}
void EndGrace(){
//   Close GRACE or pipe
    if (GraceIsOpen()) {
        /* Flush the output buffer and close Grace */
        GraceClose();
        // Close pipe but keep grace open: ne pas oublier de quitter grace sinon  les process
        //s'accumulent
        //GraceClosePipe();
        /* We are done */
        exit (EXIT_SUCCESS);
    } else {
        exit (EXIT_FAILURE);
    }
}
