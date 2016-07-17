# Programme_GP

Charger le programme charger.gp dans Pari/GP pour charger tous les programmes. Dans le dossier programme ce sont mes programmes, dans le dossier libre ce sont certains programmes de Pierre Castel et de Denis Simon.

La fonction Alea(n,m) permet de générer une matrice symétrique aléatoire de taille n avec des coefficients entre -m et m.

La fonction Qfn(A,B) permet de calculer un zéro rationnel commun à deux formes quadratiques, dont les matrices symétriques sont A et B, de manière probabiliste.
Pour n=13, le programme est rapide car on utilise un programme de Pierre Castel qui calcule un zéro d'une forme quadratique sans factoriser. Pour n <= 11 ce n'est pas le cas, donc pour des coefficients trop grand le temps de calcul est long.

La fonction qf13(A,B) est similaire à Qfn mais calcule le zéro de manière déterministe, donc en construisant la base pas à pas ce qui prend beaucoup plus de temps. Pour le moment je n'utilise pas le programme de Castle pour faire les changements de bases, c'est pour cela que le programme est long. Je vais changer \c ca sous peu.

