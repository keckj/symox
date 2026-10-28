# AUTO intégré R6r1 — rejeu public release

Snapshot accepté : source `80c98adae719e1cdfc2e7621a234d0ade781cdb3bc26c4285360263b0ca719ae`, natif `cc6d021087f8c3d8404011ecf59947d071dfa7ef9bb865c2ae7bd44431cc09ab`. Suite complète du parent : 50416 passes, 55 skips, 4 xfails, 50475 nœuds exact-once. Mesure lancée uniquement ensuite, dans la fenêtre calme autorisée ; aucune source modifiée.

Même worker que la publication privée : 14 cas publics, huit matrices distinctes, deux passages par matrice, trois processus par bras sur CPU24. R5r1 constitue la baseline ; les deux autres bras utilisent exactement le même binaire R6r1 avec AUTO désactivé/actif. 2016 certificats rationnels exacts et tous les choix/digests sont identiques. Les 3018 hashes R6r1 ont été vérifiés avant/après. Préflight charge/processus et services ordinaires conservés dans `run-r6r1.log`.

Temps du second passage, médianes des trois médianes de processus. L’appel public complet est chronométré ; construction et certificats sont hors chronométrage. Le premier passage suit la certification des choix : aucune affirmation de processus entièrement froid. Les données des deux passages restent dans les JSON bruts.

| Cas | R5 (µs) | R6 AUTO désactivé (µs) | R6 AUTO actif (µs) | R5/R6 | R6 désactivé/actif |
|---|---:|---:|---:|---:|---:|
| QQ-CD-3-dense-rref | 314.3 | 296.7 | 305.5 | 1.029 | 0.971 |
| QQ-CD-3-dense-rref_den | 306.8 | 297.3 | 304.3 | 1.008 | 0.977 |
| QQ-CD-10-dense-rref | 943.3 | 929.1 | 708.3 | 1.332 | 1.312 |
| QQ-CD-10-dense-rref_den | 976.1 | 932.1 | 706.5 | 1.382 | 1.319 |
| QQ-CD-20-dense-rref | 3446.5 | 3435.7 | 2534.5 | 1.360 | 1.356 |
| QQ-CD-20-dense-rref_den | 3422.9 | 3452.4 | 2507.1 | 1.365 | 1.377 |
| QQ-CD-10-sparse-rref | 785.3 | 775.5 | 636.1 | 1.234 | 1.219 |
| QQ-CD-10-sparse-rref_den | 768.1 | 778.4 | 632.0 | 1.215 | 1.232 |
| QQ-FF-10-dense-rref | 4438.6 | 4440.2 | 4166.4 | 1.065 | 1.066 |
| QQ-FF-10-dense-rref_den | 4343.1 | 4286.1 | 4069.8 | 1.067 | 1.053 |
| ZZ-FF-10-dense-rref | 964.9 | 952.2 | 856.6 | 1.127 | 1.112 |
| ZZ-FF-10-dense-rref_den | 777.8 | 778.3 | 681.2 | 1.142 | 1.143 |
| QQ-GJ-10-dense-rref | 2305.1 | 2302.1 | 2145.9 | 1.074 | 1.073 |
| QQ-GJ-10-dense-rref_den | 2298.8 | 2284.9 | 2138.8 | 1.075 | 1.068 |

QQ/CD dense gagne environ 1,33–1,38× à 10×10 et 1,36× à 20×20 face à R5 ; le contrôle même binaire confirme environ 1,31–1,38×. Sparse 10×10 gagne 1,21–1,23× ; autres branches 1,06–1,14×. Le petit 3×3 montre un coût du dispatch AUTO dans le contrôle même binaire : 296,7→305,5µs(+3,0 %) et 297,3→304,3µs(+2,3 %). Les chiffres face à R5 fluctuent de quelques pourcents ; aucun gain n’est revendiqué sur ce petit cas. Les surcoûts et la publication privée antérieure restent visibles ; aucune optimisation n’a été mêlée à la livraison pour les masquer.

La preuve complète de chaque résultat suit la publication privée : dominance stricte rationnelle par lignes, donc inversibilité, sortie exactement I ou den*I avec den non nul, pivots complets, domaine/format et entrées inchangés. Aucune tolérance, erreur ou résultat unsupported chronométré comme gain.

Bruts/commandes/identités : `results/*.json`, `results/*.stderr`, `results/commands.json`, `results/identity.json`. Hashes : `results-manifest.json`. Agrégat : `summary-r6r1.json`. Worker inchangé : `worker.py` ; pilote : `run.py`. La fenêtre est terminée et Demeter rendu au parent.
