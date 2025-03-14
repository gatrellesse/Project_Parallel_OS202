# Project_Parallel_OS202

#Étudiants

Gabriel BAPTISTA TRELLESSE

Guilherme GELMI SALVO

Rian RADECK SANTOS COSTA

## Première étape

Lien vers le source: https://github.com/gatrellesse/Project_Parallel_OS202/tree/premiere

```bash=
lscpu
```

| **CPU(s)** | **L1d** | **L1i** | **L2** | **L3** |
|:----------:|:-------:|:-------:|:------:|:------:|
|     16     |  384KB  |  256KB  |  10MB  |  20MB  |

Pour ce faire, nous avons dû transformer le dictionnaire contenant les case brûlante en un vecteur. Ceci est nécessaire car la directive openmp ne supporte pas les conteneurs sans « itérateurs aléatoires ». Nous décomposons la boucle de répétition en 3 autres boucles. La première est chargée de calculer quelles nouvelles cellules sont affectées par le feu et d'atténuer celles qui brûlent déjà. Pour éviter les « race conditions », nous stockons les cellules qui ont commencé à brûler dans un vecteur local que nous transformons ensuite en un vecteur global à l'aide d'une session critique. Les deux autres calculent les changements dans les cartes en fonction des maisons qui brûlent.

Nous vérifions que la solution est toujours la même en utilisant un fichier txt contenant l'état de chaque maison pour chaque itération. Pour que la solution parallèle soit toujours la même que la solution séquentielle, nous avons apporté une petite modification au code séquentielle, de sorte que son résultat ne dépende plus de l'ordre d'accès aux cases brûlantes.

Pour calculer les accélérations pour cette question et les questions suivantes, nous avons utilisé un script bash qui calcule la moyenne des valeurs d'un certain nombre d'exécutions pour chaque configuration de thread. 

![speedup_analysis](https://github.com/gatrellesse/Project_Parallel_OS202/blob/premiere/projet/src/speedup_analysis.png?raw=true)

| ENV_VALUE | n_iterations | avg_global_time | avg_time_update | avg_time_affichage |
|-----------|--------------|-----------------|-----------------|--------------------|
| 1         | 1209         | 12.4741200000   | 0.0010559902    | 0.0092466526       |
| 2         | 1209         | 11.9475800000   | 0.0008583218    | 0.0090070397       |
| 3         | 1209         | 11.6225000000   | 0.0006830563    | 0.0089171733       |
| 4         | 1209         | 11.3838600000   | 0.0005821306    | 0.0088222867       |
| 5         | 1209         | 11.4101800000   | 0.0005322900    | 0.0088932199       |
| 6         | 1209         | 11.8211800000   | 0.0005197844    | 0.0092456349       |
| 7         | 1209         | 12.0260000000   | 0.0005352615    | 0.0093990784       |
| 8         | 1209         | 12.3651800000   | 0.0005293415    | 0.0096853244       |
| 9         | 1209         | 12.5437800000   | 0.0005118092    | 0.0098505788       |
| 10        | 1209         | 12.9939400000   | 0.0004932867    | 0.0102412704       |
| 11        | 1209         | 13.2976800000   | 0.0004748270    | 0.0105111460       |
| 12        | 1209         | 13.4689400000   | 0.0004696056    | 0.0106576697       |
| 13        | 1209         | 13.6782800000   | 0.0004564326    | 0.0108439811       |
| 14        | 1209         | 14.2681000000   | 0.0004877520    | 0.0113008417       |
| 15        | 1209         | 14.9260200000   | 0.0007174380    | 0.0116153505       |
| 16        | 1209         | 18.8089400000   | 0.0023888039    | 0.0131554469       |

On constate que l'accélération du temps de mise à jour augmente avec le nombre de threads utilisés, mais qu'à partir d'un certain nombre de threads, les accès synchronisés aux zones critiques entraînent probablement un surcoût. En outre, l'amélioration de l'accélération globale est faible car le temps d'affichage est relativement plus long que le temps de mise à jour.

## Deuxième étape

Lien vers le source: https://github.com/gatrellesse/Project_Parallel_OS202/tree/deuxieme

Après avoir initialisé l'environnement MPI, nous utilisons le rang des processus pour différencier les tâches qui leur sont assignées. Nous utilisons un vecteur tampon et une répartition non bloquante afin de pouvoir calculer l'étape suivante pendant la répartition.

| ENV_VALUE | n_iterations | avg_global_time | avg_time_update | avg_time_affichage | avg_time_for |
|-----------|--------------|-----------------|-----------------|--------------------|--------------|
| 1         | 1211         | 11.3218000000   | 0.0009974322    | 0.0092888561       | 0.0000000000 |

![time_analysis](https://github.com/gatrellesse/Project_Parallel_OS202/blob/deuxieme/projet/src/time_analysis.png?raw=true)

En analysant les résultats obtenus, on remarque que le temps total global moyen par itération (avg_global_time) a diminué par rapport au programme séquentiel, passant de 12.47412 à 11.32180, soit une amélioration de **9,24 %** du temps d'exécution.

Les résultats ne sont pas meilleurs, car l'affichage constitue un goulot d'étranglement, étant beaucoup plus lent que les calculs. Ainsi, le calcul de l'itération suivante reste en attente dans un buffer jusqu'à la fin de l'affichage.

## Troisième étape

Lien vers le source: https://github.com/gatrellesse/Project_Parallel_OS202/tree/troisieme

Pour cette étape, nous avons simplement fusionné le code des étapes précédentes et ajusté les commandes d'exécution.

![speedup_analysis](https://github.com/gatrellesse/Project_Parallel_OS202/blob/troisieme/projet/src/speedup_analysis.png?raw=true)

On peut observer dans les deux graphes ci-dessus qu'en ajoutant la parallélisation avec OpenMP, nous avons obtenu de bonnes accélérations du temps de mise à jour (calcul). Cependant, l'accélération globale reste la même, ce qui confirme que la contrainte principale provient du temps d'affichage, qui est plus long que le calcul. 

Ainsi, même si nous rendons les calculs beaucoup plus rapides, cela ne changera pas le temps d'exécution global, car celui-ci restera limité par la lenteur de l'affichage. 

On peut également remarquer que l'accélération diminue avant d'atteindre 16 cœurs, ce qui correspond une légère saturation de la bande passante mémoire avec l’augmentation de nombre de cœurs.

## Quatrième étape

Lien vers le source: https://github.com/gatrellesse/Project_Parallel_OS202/tree/quatrieme

La dernière étape consiste à définir un groupe de communication pour les processus de rank ≠ 0, responsables du calcul. Chaque processus est chargé d'une tranche du terrain et doit envoyer et recevoir des cellules fantômes situées aux frontières entre les tranches, où le feu se propage. En calculant les accélérations globales et l'avancement dans le temps, nous avons obtenu les résultats suivants :

![speedup_analysis](https://github.com/gatrellesse/Project_Parallel_OS202/blob/quatrieme/projet/src/speedup_analysis.png?raw=true)

On a obtenu des résultats très similaires à l'étape précédente, mais il y a quelques observations à faire. Par exemple, la courbe de croissance de l'accélération du temps de mise à jour (calcul) est plus linéaire qu'avec OpenMP. 

Cela s'explique par le fait que MPI utilise des systèmes distribués, où chaque cœur possède sa propre mémoire. Cela permet d'éviter les problèmes de contention et de cohérence de cache qui peuvent survenir avec OpenMP lorsque l'on ajoute plusieurs cœurs. 

Il est également important de noter que le speed-up maximal est atteint avec seulement 7 cœurs, alors que dans le cas précédent, un speed-up équivalent nécessitait environ 11 cœurs

Derniere remarque, notre code n'est pas assyncrono entre les 2 groupes de communications(Affichage, Calcul), ils utilisent des communications blocantes entre eux, on a essayé de faire tourner, mais il n'a pas fonctionné très bien. Pour cette application, ça ne vas pas donner des problémes reels grace au goulout d'affichage, mais dans un autre cas de simulation que les calculs sont plus couteses ça serait important de fonctionner.
