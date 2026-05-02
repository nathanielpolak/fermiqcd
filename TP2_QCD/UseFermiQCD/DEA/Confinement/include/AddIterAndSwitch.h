//  add iterations interactively
//
	if(k==niter and niter<maxiter)
           {
            system("xmessage -print -center -buttons stop,+10?,+100?,+500?,+1000?  want to add iterations?>message ");
            ifstream tmpmess("message");
            string answer;
            getline(tmpmess,answer);
            tmpmess.close();
            cout<<answer<<endl;
            int increase;
            if(answer=="stop") increase=0;
            if(answer=="+10?") increase=10;
            if(answer=="+100?") increase=100;
            if(answer=="+500?") increase=500;
            if(answer=="+1000?") increase=1000;
            niter=niter+increase;
            if(niter>maxiter) niter=maxiter;
//
            system("xmessage -print -center -buttons  CPP,CPU,GPU  choose system:>newmessage ");
            ifstream newtmpmess("newmessage");
            string newanswer;
            getline(newtmpmess,newanswer);
            newtmpmess.close();
            cout<<"system: "<<newanswer<<endl;
           if(newanswer=="CPU") systemchoice=0;
           if(newanswer=="GPU") systemchoice=1;
           if(newanswer=="CPP") systemchoice=2;
           printf(" now using %d",systemchoice);
           }