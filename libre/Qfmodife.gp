\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\ Nouvelle version de myqfsolve...
\\ Améliorée, meilleures notations, design épuré et minimaliste :)
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\Auteur : Pierre Castel -> pierre.castel@math.unicaen.fr

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\                               \\
\\  VERSION DU 07 / 02 / 2011    \\
\\                               \\ 
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\ /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
\\ /!\
\\ /!\ Ce fichier requiert qfsolve.gp (www.math.unicaen.fr/~simon/qfsolve.gp) et 
\\ /!\ smooth_things.gp (??)
\\ /!\
\\ /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
print("");
print("===== mynewqfsolve.gp =====");
print("- completemat(Q,pow,coeff,prim=0,newsgn=1,detQ=matdet(Q))");
print("- myQfsolve(Qf,pow=1,prim=0)");
print("");


print("Chargement de smooth_things.gp");
\r smooth_things.gp
print("Chargement de qfsolve.gp (Simon)");
\r qfsolve.gp
print("Chargement de poll_sch.gp");
\r poll_sch.gp
print("Chargement de gramschmidt.gp");
\r gramschmidt.gp
print("Chargement de minim_p.gp");
\r minim_p.gp

DEBUGLEVEL_myqfsolve=0;

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
\\  DEBUT DES SOURCES                
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\Q = matrice à compléter
\\pow = exposant pour les bornes
\\coeff = borne pour la taille des coefficients du vecteur complétion
\\prim = drapeau pour savoir si on teste avec un "isprime()" lorsque prim=0
\\      ou alors avec un "ispseudoprime()" lorsque prim=1.
\\newsgn = signe du determinant de la matrice renvoyée.
\\         permet de gérer les cas du type [1,xx] ou [xx,1].

{completemat(Q,pow,coeff,prim=0,newsgn=1,detQ=matdet(Q))=
    local(detn,n,Qt,bor,QP1,i,j,flag,count,Y,H,Hred,nsp_Hred,Hred2,ispr,z,p,m,test_mat,GS,Yt,val,factest,dettest);
    
    if(DEBUGLEVEL_myqfsolve>0,
        print("Entree dans completemat avec : ");
        print("Q=",Q);
        print("pow=",pow);
        print("coeff=",coeff);
        print("prim=",prim);
        print("sgn=",newsgn);
        print("detQ=",detQ);
    );
    
    \\Précalculs
    detn = detQ;
    n = length(Q);
    Qt = matadjoint(Q)~;\\"Dois-je transposer ?--> Qt est symétrique !!!"
    
    \\Borne pour la taille des coefficients
    bor = abs(coeff);
    
    \\Contruction de la future matrice complétée:
    QP1=matrix(n+1,n+1);
    for(i=1,n,
        for(j=1,n,
            QP1[i,j]=Q[i,j];
        );
    ); 

    \\Initialisation de certaines variables
    flag = 0; \\drapeau pour le while
    count = 0;\\compteur du nombre d'essais

    \\Boucle principale :
    \\on tire un vecteur aleatoirement, puis on 
    \\fait ce qu'il faut pour que le determinant de QP1
    \\soit "smooth" ou bien "almost_smooth"
    while(!flag,
        count++;

        \\choix du vecteur
        Y=vector(n,i,random(bor)-bor>>1)~;
        H=Y~*Qt*Y;
        print1(count," ");
        if(DEBUGLEVEL_myqfsolve>0,
            print("Essai de completion:",count);
            print("Y=",Y);
            print("H=",H);
        );
        
        \\changement de signe selon la valeur de la signature
        if(newsgn==1,
            Hred = H%(detn^pow);
            if(Hred==0,next(););
            nsp_Hred = non_smooth_part(Hred);
            z=(H-Hred)/detn;
        ,\\sinon on doit ajouter un -
            Hred = lift(Mod(H,detn^pow))-abs(detn^pow);
            if(Hred==0,next(););
            nsp_Hred = non_smooth_part(Hred);
            z=(H-Hred)/detn;
        );\\fin condition signature 
        \\on élimine déjà certaines valeurs
        v2 = valuation(Hred,2);
        if(v2%2!=1,next(););
        Hred2 = Hred>>v2;
        if(prim==0,
            \\On utilise un "isprime()"
            if( abs(nsp_Hred)==1 || ispr=isprime(abs(nsp_Hred)),
                if(abs(nsp_Hred)!=1, addprimes(abs(nsp_Hred)));
                
                flag=1;
            );\\fin if "isprime"
        ,\\sinon; On utilise un "ispseudoprime()"
            if( abs(nsp_Hred)==1 || ispr=ispseudoprime(abs(nsp_Hred)),
                if(abs(nsp_Hred)!=1, addprimes(abs(nsp_Hred)));
                
                flag=1;
            );\\fin if "ispseudoprime"
        );\\fin du if prim
    );\\fin Boucle principale (while)
    
    \\On a un "bon" vecteur Y, donc on complète la matrice
    for(i=1,n,
        QP1[n+1,i]=Y[i];
        QP1[i,n+1]=Y[i];
    );
    QP1[n+1,n+1]=z;
    test_mat=matrix(n-1,n-1);
    for(i=1,n-1,
        for(j=1,n-1,
            test_mat[i,j]=Q[i,j];
        );
    );
    return(QP1);
}\\fin completemat

\\Donne une solution de l'équation X~*Q*X=0, où Q est la matrice d'une forme quadratique
\\ en 5 variables ou plus. Si prim=0, teste la primalité avec isprime(.), si prim=1, 
\\teste la primalité avec un ispseudoprime(.). pow est l'exposant pour la borne des coeffs

{myqfsolve(Qf,pow=1,prim=0)=
    local(flag2,Gmin,Min,Q,flagmin,SNF,dimker,pker,sgn,sgtmp,newsgn,detQ,borne,n,QP1,S1,B1,SW1,QP1nb,LLLg1,QP1nb2,B2,SW2,QP1nb3,CGB1,QM2,S2,B3,SW3,QM2nb,LLLg2,QM2nb2,B4,SW4,QM2nb3,CGB2,CGB2P2,CGB,DB,CB,E,F,S,St);
    if(DEBUGLEVEL_myqfsolve>0,
     print("Entree dans myQfsolve avec:");
     print("Q=",Q);
     print("pow=",pow);
     print("prim=",prim);
    );
    if(abs(matdet(Qf))<=10000,return(qfsolve(Qf)));
    Q=Qf;
    \\Eventuelle minimisation
    Min=is_minimisable(Q);
    Gmin=matid(length(Q));
    if(Min!=0,flagmin=1;);
    while(Min!=0,
        Q=Min[1];
        Gmin=Gmin*Min[2];
        Min=is_minimisable(Q);
    );
    if(valuation(matdet(Q),2)>=1,
        Min=minim_en_2(Q);
        Q=Min[1];
        Gmin=Gmin*Min[2];
        flagmin=1;
    );
    \\précalculs
    detQ = matdet(Q);
    if(abs(detQ)<10^3,
        my (s = Gmin * qfsolve(Q));
        return (s / content(s))
    );
    sgn = qfsign(Q);
    borne = abs(detQ)^pow;
    n=length(Q)+1;
    
    \\initialisation
    if(sgn[1]*sgn[2]==0,
        print("Il n'y a pas de solutions. La signature est : ",sgn);
        return(0);
    );
    \\pour que dans la suite, on ait dans le cas le pire toujours un - à ajouter
    if(sgn[1]<sgn[2]    ,
        Q=-Q;
        sgtmp=sgn[1];
        sgn[1]=sgn[2];
        sgn[2]=sgtmp;
        detQ = -detQ;
    );

    \\choix du signe a ajouter dans la future matrice
    newsgn=1;
    if(sgn[2]==1,
        newsgn=-1;
    );

    \\complétion de la matrice Q
    flag2=0;
    while(!flag2,
      flag2=1;
      QP1 = completemat(Q,pow,borne,prim,newsgn,detQ);
      
      \\Calcul d'une solution
      S1 = qfsolve(QP1);
      S1 *= denominator(S1);
      
      \\S1 vérifie S1~*QP1*S1=0
      \\on s'en sert pour construire une base
      B1 = completebasis(S1);
      
      \\cette procédure met le vecteur à la fin 
      \\de la base, donc on le remet au début
      \\Matrice du swap :
      SW1=matid(n);
      SW1[1,1]=SW1[n,n]=0;
      SW1[1,n]=SW1[n,1]=1;

      \\swap des vecteurs
      B1*=SW1;
      \\On écrit la matrice de QP1 dans la nouvelle base :
      QP1nb=B1~*QP1*B1;  
      \\cette matrice a un 0 en haut à gauche on la réduit avec LLLgoon
      LLLg1=LLLgoon(QP1nb);
      QP1nb2=LLLg1[1];
      B2=LLLg1[2];   
      
      \\le bloc à extraire est au milieu,
      \\on va swaper les vecteurs pour plus de facilité
      \\matrice du swap
      SW2=matid(n);
      SW2[n,2]=SW2[2,n]=1;
      SW2[n,n]=SW2[2,2]=0; 
      
      \\swap
      QP1nb3=SW2~*QP1nb2*SW2; 
      \\le changement de base de cette partie est:
      CGB1=B1*B2*SW2;
      
      \\QP1nb3 est diagonale par blocs avec un 2x2 en haut à gauche
      \\on extrait le bloc restant
      QM2=matrix(n-2,n-2);
      for(i=3,n,
          for(j=3,n,
              QM2[i-2,j-2]=QP1nb3[i,j];
          );
      );
      
      \\On fait le même genre de choses mais avec QM2
      S2=qfsolve(QM2);
      if(type(S2)=="t_INT",
          print("Pas de solution en step 2 en ",S2);
          flag2=0;
      );
    );\\fin while
    
    \\prévention en cas d'erreur
    if(type(S2)=="t_INT",
        print("Q=",Q);
        print("QM2=",QM2);
        error("Pas de solution à la seconde résolution !");
    );
    \\mise en coordonnées entières
    S2 *= denominator(S2);
    
    \\S2 vérifie S2~*QM2*S2=0
    \\on s'en sert pour construire une base
    B3=completebasis(S2);
    \\cette procédure met le vecteur à la fin 
    \\de la base, donc on le remet au début
    SW3=matid(n-2);
    SW3[1,1]=SW3[n-2,n-2]=0;
    SW3[1,n-2]=SW3[n-2,1]=1;
    
    \\swap des vecteurs
    B3*=SW3;
    
    \\On écrit la matrice de QP1 dans la nouvelle base :
    QM2nb=B3~*QM2*B3;
    \\cette matrice a un 0 en haut à gauche on la réduit avec LLLgoon
    LLLg2=LLLgoon(QM2nb);
    QM2nb2=LLLg2[1];
    B4=LLLg2[2];
    
    \\on va swaper les vecteurs pour plus de facilité
    SW4=matid(n-2);
    SW4[n-2,2]=SW4[2,n-2]=1;
    SW4[n-2,n-2]=SW4[2,2]=0;
    QM2nb3=SW4~*QM2nb2*SW4;
    \\le changement de base de cette partie est
    CGB2=B3*B4*SW4;
    \\on va remettre tous les changements de base ensemble
    CGB2P2=matid(n);
    for(i=3,n,
        for(j=3,n,
            CGB2P2[i,j]=CGB2[i-2,j-2];
        );
    );
    CGB=CGB1*CGB2P2;
    \\la matrice finale
    QP1F=CGB~*QP1*CGB;
    \\la matrice est prête, on calcule les solutions
    DB=QP1F;
    CB=CGB;
    n--;
    
    \\2 vecteurs isotropes orthogonaux :
    E=CB[,1];
    F=CB[,3];
    \\Calcul de la solution finale
    S=vector(n);
    if(E[n+1]==0,
        for(i=1,n,
            S[i]=E[i];
        );    
        S=S~;
        if(flagmin, S=Gmin*S;);
        S/=content(S);
        return(S);
    );
    if(F[n+1]==0,
        for(i=1,n,
            S[i]=F[i];
        );    
        S=S~;
        if(flagmin, S=Gmin*S;);
            S/=content(S);
        return(S);
    );
    
    \\Préparation de la solution:
    St=E*F[n+1]-F*E[n+1];
    
    \\extraction de la partie du haut
    for(i=1,n,
        S[i]=St[i];
    );
    \\transposition du vecteur
    S=S~;
    
    \\On remet le vrai signe si il a été changé
    if(flagmin, S=Gmin*S;);
    S/=content(S);
    \\Vérification
    if(DEBUGLEVEL_myqfsolve>0,
        printp("Vérification : S~*Q*S = ",S~*Qf*S);
    );
    if(S~*Qf*S!=0,error("Q=",Qf,"\nS=",S,"\nBad Solution !!!"));
    return(S);
}\\fin myQfsolve


addhelp(completemat,"completemat(Q,pow,coeff,prim=0,newsgn=1,detQ=matdet(Q)) : procedure de completion de la forme quadratique Q, pow=exposant de la borne de choix, coeff= taille des coefficients pour la completion,si prim=0 utilise un test de primalite lors de la completion si prim=1 utilise un test de pseudo-primalite (=0 par defaut), newsgn=signe a ajouter a la future forme. Renvoie la forme completee Qp ");
addhelp(myQfsolve,"myQfsolve(Qf,pow=1,prim=0) : algorithme de resolution de l'equation X~*Q*X=0, Q=forme quadratique entiere de dimension 5, pow=exposant de la borne de choix pour la procedure de completion completemat (=1 par defaut), si prim=0 utilise un test de primalite lors de la completion si prim=1 utilise un test de pseudo-primalite (=0 par defaut). Renvoie une solution S sous forme d'un vecteur");
