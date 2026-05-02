//  
            system("xmessage -print -center -buttons overwrite,no overwrite Confine.agr ?>message");
            ifstream tmpmess("message");
            string answer;
            getline(tmpmess,answer);tmpmess.close();
            if(answer=="overwrite") SaveGrace();
            if(answer=="no") EndGrace();