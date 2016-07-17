\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
\\
\\
\\
\\
\\
\\
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\Auteur : Pierre Castel -> pierre.castel@math.unicaen.fr

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\                               \\
\\  VERSION DU 31 / 01 / 2011    \\
\\                               \\ 
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

\\Fichier contenant des programmes en rapport à la friabilité

\\On place en variable globale la base des facteurs premiers
\\ainsi que le produit de ces nombres
\\la base des facteurs est Bf et leur produit Pf
\\posons nf la taille de la base

global(nf,Bf,Pf);
nf=1;
Bf=primes(nf);
Pf=prod(i=1,nf,Bf[i]);


\\maintenant on ecrit plusieurs procedures utiles
\\rapport a la friabilite

\\un test qui renvoie 1 si N est Bf-friable, 0 sinon
{is_smooth(N)=
    local(g);
    if(N==0,
        print("le nombre à tester est nul !");
        return(1);
    );
    g = gcd(N,Pf);
    while(g!=1,
        N\=g;
        g=gcd(N,g);
    );
    return(abs(N)==1);
}

\\un test de presque friabilité
\\ie : un friable x un premier qui n'est pas dans la base:
{is_almost_smooth(N)=
    local(g);
    if(N==0,
        print("le nombre à tester est nul !");
        return(1);
    );
    g = gcd(N,Pf);
    while(g!=1,
        N\=g;
        g=gcd(N,g);
    );
    return(abs(N)==1 || isprime(abs(N)));
}





\\une fonction qui renvoie la partie friable de N sur Bf
{smooth_part(N)=
    local(g,sp);
    sp=1;
    if(N==0,
        print("le nombre à tester est nul !");
        return(1);
    );
    g = gcd(N,Pf);
    while(g!=1,
        N\=g;
        sp=g*sp;
        g=gcd(N,Pf);
    );
    return(sp);
}

\\enfin, la complémentaire, la fonction qui renvoie la partie
\\non friable de N sur Bf
{non_smooth_part(N)=
    local(g,nsp);
    if(N==0,
        print("le nombre à tester est nul !");
        return(1);
    );
    g = gcd(N,Pf);
    while(g!=1,
        N\=g;
        g = gcd(N,Pf);
    );
    return(N);
}

addhelp(nf,"nf : nombre de nombres premiers de la base Bf. Pour changer la base, changer nf dans le fichier puis le recharger");
addhelp(Bf,"Bf : vecteur contenant les premiers de la base de friabilite. Pour changer la base, voir nf");
addhelp(is_smooth,"is_smooth(N) : renvoie 1 si N est Bf-friable, 0 sinon.");
addhelp(is_almost_smooth,"is_almost_smooth(N) : renvoie 1 si N est Bf-presque friable (ie: la partie non friable est un nombre premier), 0 sinon. ");
addhelp(smooth_part,"smooth_part(N) : donne la partie friable de N sur Bf.");
addhelp(non_smooth_part,"non_smooth_part(N) : donne la partie non friable de n sur Bf.");








