# AUTO privé : mesures des appels publics sur release

Ce résultat concerne le paquet AUTO privé figé, pas le R6 intégré. Aucun fichier livré ni binaire DEBUG n’a été modifié pour ces mesures.

Base R5r1 immuable : source `2955abf8b02c1e520ddb921d17e56614b8f4c5d391b1fb87730118881665dbe0`, natif `532386bda46046d3d21c99f49c1b863fac62bdee5c357d03cf5123bde6913d4b`. Après privé : sources du patch `ef6aa735…` sur base 2eb, release native `8eed6deb5f98078e57d4128e63cf567e03db9f8ab562c708d598e52fa913c0e8`, manifeste release `13ac37cfac07e97c8d257d9e9825bb028817edea664f92c14c12b41642b6ceb8` ; les 505 tests nouveaux passent aussi sur cette release.

Le noyau privé provient d’une base antérieure à R5. Un troisième bras, même binaire privé avec seulement l’adaptateur AUTO désactivé, isole donc le coût de ce port des autres différences de noyau. Ce contrôle est indispensable pour attribuer un gain à AUTO.

Neuf processus frais alternés sur Demeter CPU24, trois par bras ; 14 cas, huit matrices distinctes par cas et deux passages par matrice. Les choix sont certifiés avant chronométrage. Le premier passage mesure donc le premier appel RREF sur chaque matrice après cette préparation, sans prétention de processus entièrement froid. Le second passage réutilise ces huit matrices. Les temps ci-dessous sont les médianes des trois médianes de processus, second passage ; construction des entrées et vérifications exactes hors chronométrage, appel public `.rref()` ou `.rref_den()` complet dans le chronométrage.

| Cas | R5 (µs) | Privé AUTO désactivé (µs) | Privé AUTO actif (µs) | R5/actif | Désactivé/actif |
|---|---:|---:|---:|---:|---:|
| QQ-CD-3-dense-rref | 299.7 | 313.4 | 316.2 | 0.948 | 0.991 |
| QQ-CD-3-dense-rref_den | 310.9 | 320.1 | 320.5 | 0.970 | 0.999 |
| QQ-CD-10-dense-rref | 957.6 | 969.7 | 705.6 | 1.357 | 1.374 |
| QQ-CD-10-dense-rref_den | 928.2 | 926.8 | 710.1 | 1.307 | 1.305 |
| QQ-CD-20-dense-rref | 3473.5 | 3470.5 | 2497.7 | 1.391 | 1.390 |
| QQ-CD-20-dense-rref_den | 3444.7 | 3436.7 | 2470.1 | 1.395 | 1.391 |
| QQ-CD-10-sparse-rref | 801.4 | 807.7 | 631.7 | 1.269 | 1.279 |
| QQ-CD-10-sparse-rref_den | 769.7 | 763.6 | 655.3 | 1.175 | 1.165 |
| QQ-FF-10-dense-rref | 4490.1 | 4473.3 | 4257.0 | 1.055 | 1.051 |
| QQ-FF-10-dense-rref_den | 4348.2 | 4352.8 | 4242.7 | 1.025 | 1.026 |
| ZZ-FF-10-dense-rref | 998.5 | 972.2 | 859.7 | 1.162 | 1.131 |
| ZZ-FF-10-dense-rref_den | 782.0 | 815.4 | 682.5 | 1.146 | 1.195 |
| QQ-GJ-10-dense-rref | 2362.6 | 2302.5 | 2271.9 | 1.040 | 1.013 |
| QQ-GJ-10-dense-rref_den | 2363.4 | 2317.8 | 2207.0 | 1.071 | 1.050 |

Le gain public le plus net est QQ/CD dense : 1,31–1,36× à 10×10 et 1,39× à 20×20. Les autres branches ont des gains plus modestes. Le petit cas QQ/CD 3×3 conserve un surcoût observé face à R5 : 5,5 % sur rref et 3,1 % sur rref_den. Le contrôle même binaire est presque neutre sur ces échantillons (+0,9 % et +0,1 %) : aucun gain n’est annoncé sur ce petit cas. Les écarts de quelques pourcents restent à considérer avec la dispersion brute ; aucune optimisation supplémentaire n’a été mêlée au paquet livré.

Les 2016 certificats et tous les choix/digests sont identiques entre les trois bras. Chaque matrice rationnelle est strictement dominante par lignes, vérifié exactement avec Fraction. Si un vecteur du noyau était non nul, sa composante de module maximum donnerait une contradiction par cette dominance ; les matrices sont donc inversibles. Le résultat doit être exactement I, ou den*I avec den non nul pour rref_den, avec tous les pivots. Les domaines, formats, dénominateurs, entrées et deux passages sont également vérifiés ; aucune tolérance numérique ni temps d’échec/unsupported compté comme gain.

Données brutes et commandes : `raw/*.json`, `raw/*.stderr`, `raw/commands.json`, `raw/worker.py`, `raw/run.py` ; hashes dans `raw-manifest.json`. Les chemins de packages, versions Python, affinité et hashes natifs sont enregistrés et contrôlés pour chaque processus. `summary.json` conserve les temps des deux passages. Fenêtre terminée et hôte rendu au parent.

Le même rejeu est préparé pour R6 intégré, uniquement après réception de son snapshot/natif et fin de sa suite complète, dans une nouvelle fenêtre calme. Aucune conclusion R6 ne découle des chiffres privés ci-dessus.
