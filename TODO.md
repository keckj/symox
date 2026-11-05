# TODO — Capacités mathématiques et parité SymPy

Mis à jour le **28 octobre 2026**, à partir de `main` à **`5eb39d6`**
(dernier commit du 24 octobre), avec les lots des 26–28 octobre dans l’arbre de travail.
Référence : **SymPy 1.14.0**. Dernier lot validé sur Demeter le 28 octobre.

Dernière suite complète : binaire `e2b10b0c`, sources `cc617a4a`
(**56 311 tests**, Demeter) ; 762 tests Rust, 1 test d’intégration
(651 cas LLL) et 5 doctests réussis, clippy 1.99 validé, bindings compris.
Il reste quatre xfail connus. Les produits sparse et
les noyaux DDM de multiplication et
Gauss-Jordan sont intégrés. Calcul sparse, fraction-free, `log_to_atan`, payloads
EX Float0, dénominateurs QQ et déterminant DDM sont validés dans leurs domaines.
LRT QQ/ZZ, Berkowitz dense ZZ/QQ, la limite 22334 et les branches affines
Laplace E1 passent aussi la suite complète. Les dérivées d'ordre symbolique
restent partielles ; les grands ordres littéraux certifiés sont validés.
PolyMatrix ground QQ/corps de fractions et log_to_real rationnel de degré 1/2
sont validés. Le pilote rationnel QQ/ZZ admis, les statistiques/décisions AUTO
ZZ/QQ et Berkowitz sparse sélectionné passent aussi la suite complète ; les autres
contextes conservent leurs parcours complets. Solveur inexact stock,
arrondis CC certifiés et égalité composite typée sont validés sur `8a8a8e12`.
LLL public stock et arrondis exacts QQ/ZZ/RR sont validés dans leurs routes
admises ; le vrai premier appel LLL froid est aussi natif dans le préfixe
stock désormais prouvé. Les autres protocoles gardent leur calcul complet.
Gauss-Jordan sparse stock QQ est également natif et validé, avec métadonnées
ordonnées et identités de clés conservées ; les autres routes restent complètes.
LU stock QQ et la réduction polynomiale Risch QQ/ZZ admise sont natives et
validées ; le rejet LLL froid évite le parcours lourd sans classifier le cache.
Les substitutions triangulaires stock QQ sont aussi natives ; les autres
domaines restent à porter. Le refus Risch sans étape
et la lecture Unicode des métadonnées LLL sont corrigés dans les routes testées.
HNF ordinaire ZZ et les gardes LU/solve sans événements d’audit ajoutés sont
validés ; FFLU ZZ/QQ et ses quatre phases sont aussi natives et validées.
HNF modulo D stock ZZ et splitfactor Risch QQ fermé sont aussi natifs.
HNF numberfield, Smith et les autres domaines FFLU restent à porter.
Le corpus d'intégration est rejoué sur `e7b29759`, les EDO sur `570f1832`,
les exemples documentés sur `576a9f76`.

Objectif : atteindre au moins les capacités mathématiques de la référence,
avec les mêmes domaines acceptés, l'arithmétique exacte, la complétude des
solutions, la précision demandée et un comportement d'échec explicite.
La couverture des noms d'API et les performances d'exemples choisis ne
prouvent pas la parité. **La parité algorithmique globale n'est pas certifiée.**

Ce fichier suit le travail restant. Les bilans des vagues 1 à 17 et les
anciennes listes sont conservés dans l'[archive du TODO](docs/history/TODO-before-2026-10-26.md).
Les preuves, sondes et corrections intentionnelles sont dans le
[contrat de parité](docs/algorithmic-parity.md). Les extensions et sous-classes
Python restent un chantier distinct, hors de ce TODO mathématique.

## Reprise sur Demeter (état au 30 octobre)

Main intègre la vague 38 (branche `integrate-wave38-20261030`). Dernière
validation : suite complète **60 020 réussis, 0 échec** ("SUITE CLEAN",
60 085 collectés et exécutés, 682 fichiers), cargo 852 + 1 + 5, clippy
correctness propre. Les corpus intégration et ODE n'ont pas été rejoués sur
la vague 38 (le score ODE 303/279/274/315 ci-dessus est périmé : les agents
mesurent 315/285/267 sur 21068f42 avec un script reconstruit).

Vague 38 (calcul en Rust, Python limité à l'encodage et la reconstruction) :
- arr2mat : `convert_array_to_matrix`, ses passes, `_EditArrayContraction`
  et `_ArgE` natifs (`arrayexpr_a2m.rs`) ; tests SymPy des tableaux 151/201,
  corpus 517/527 ; 86/527 cas retombent encore sur la version vendue
  (lacune) ; `_chain_array_to_matrix` des dérivées matricielles reste en
  Python (lacune). `Trace(MatMul)` n'est plus canonicalisé (aligné sur
  SymPy, validé le 30 octobre).
- arrindex : `array_derive`, `convert_indexed_to_array` et le plan
  d'indices de `convert_array_to_indexed` natifs (`arrayexpr_index.rs`),
  corpus 162/168.
- odepart : ensemble d'essai et résolution des coefficients indéterminés,
  variation des paramètres (`ode/undetermined.rs`, `ode/variation.rs`) ;
  correctif de parité `dotprodsimp` (gardé, validé le 30 octobre, mais
  vop_013 dépasse 6 Go : lacune). simplify/trigsimp, évaluation des
  intégrales et enveloppes Euler restent en Python.
- smallgaps : déterminant Bareiss Float (77/77), checksol 17799, seconde
  branche puissance-de-E (`ops.rs`), unités, correctif clippy.
- limgaps : corrections limites/séries en Rust, 1118/1293 asserts, déclins
  Gruntz 92 → 90 ; `factorial(x)/x` et `gamma(x)/x` en -oo laissés non
  évalués (SymPy donne 0 ; écart validé le 30 octobre).
- prdelog : cas de base de `is_log_deriv_k_t_radical_in_field` en Rust
  (`integrate/prde_log.rs`, 938 appels, 0 écart) ; les cas non de base
  restent en Python (lacune).

Interrompu, à reprendre : ratintex (log_to_atan sur EX, RootSum, exroots),
worktree `~/symox-ci/wt/ratintex-20261030k`, travail non commité et non
validé, branche `feat-ratintex-20261030k` encore à 21068f42.

Écarts volontaires en attente de validation utilisateur : Markov spec7 1/6,
conditionnement multi-temps markovrust (22/28 cas où SymPy se trompe) et
recherche de période ; radicalrust réseaux/ConditionSet ; holonomic y0 ;
derivsym ; seriesrest ; Derivative(Array) dans tracearr et arrayrust ;
ratint 1/(x**8+16) ; limite unilatérale Piecewise ; limit(atan(a*x)/x, x, oo)
= 0 ; x**(n*y+3) à n=0 ; Mod dans solveset ; décalage Laplace ;
piecewise_fold w23rust ; limit(subfactorial(x)/factorial(x), x, oo) laissé
non évalué (SymPy donne 0, valeur vraie 1/e) ; règle Abs/sign(c) de
limitsrest r2 identique à SymPy même pour y complexe ; trois asserts de
`tests/test_solver_checks.py` alignés sur SymPy (liste vide True,
UserWarning, flags inconnus ignorés).

Lacunes ouvertes principales (détail dans `docs/algorithmic-parity.md`) :
- ODE : simplification et intégrales des solutions particulières en
  Python ; vop_013 au-delà de 6 Go ; collect à exposants partiels ; parité
  d'erreur du casus irreducibilis ; factorable, 2F1, Bessel, nth_algebraic,
  lie_group.
- prde : is_log_deriv_k_t_radical_in_field hors cas de base, parametric_log_deriv,
  limited_integrate ; surcoût d'admission de splitfactor.
- ratint/exroots (ratintex, 30 octobre) : TypeError QQ<sqrt(2)> corrigé
  pour un seul corps (corps différents et domaines paramétriques encore en
  TypeError) ; log_to_atan sur EX en Rust (choix du quotient le plus court porté en Rust) ; ordre
  d'impression RootSum corrigé (constructeur RootSum encore en Python) ;
  roots cas 20 corrigé. Restent : formes cubiques (simplify natif dans
  roots_quadratic sur EX), binômes à radicaux imbriqués (hypothèses
  is_extended_real/is_positive, écart en attente de validation), ordre des
  étapes d'expand.
- Limites : 175/1293 asserts encore différents ; 90 déclins Gruntz.
- Matrices : perf lambdarr gamma/poly2.
- Tableaux : 50/201 tests SymPy de tableaux, repli vendu de
  `convert_array_to_matrix` (86/527), `_chain_array_to_matrix` en Python.

Procédure de validation (leçons du 30 octobre) : n'accepter que "SUITE
CLEAN" (le compte de lots "passed" incluait des lots en échec) ; le build
doit sortir en succès avant de copier `lib_core.so` (un build tué en
silence laissait l'ancien `.so`) ; doctests parallèles à relancer en
mono-thread ("Couldn't compile") ; tests Python en worktree avec
`PYTHONPATH=python` ; jamais de `pkill` large.

## État des preuves et prochaine étape

Le critère de clôture du tri de l'audit, décidé le 22 octobre, reste : ne
laisser que des différences de forme, des corrections intentionnelles ou des
lacunes explicitement documentées comme non prises en charge. Documenter une
lacune clôt son tri ; cela ne donne pas la capacité manquante et ne prouve
pas l'équivalence avec SymPy.

| Preuve | Dernier état consigné | Portée |
| --- | --- | --- |
| Suite complète du 28 octobre, ports Risch/LLL/HNF (48 lots sur Demeter) | **56 252 réussis / 55 ignorés / 4 xfail / 0 échec inattendu**, sur 56 311 tests collectés | Tous les 54 428 anciens identifiants et les motifs skip/xfail préservés ; exactement une fois |
| Suite complète du 28 octobre, FFLU (48 lots sur Demeter) | **54 369 réussis / 55 ignorés / 4 xfail / 0 échec inattendu**, sur 54 428 tests collectés | Tous les anciens tests préservés ; extras/compilateurs actifs ; exactement une fois |
| Suite complète du 28 octobre, HNF/audit (48 lots sur Demeter) | **52 504 réussis / 55 ignorés / 4 xfail / 0 échec inattendu**, sur 52 563 tests collectés | Extras inclus ; collection/exécution exactement une fois ; compilateurs actifs |
| Suite complète du 28 octobre, solve/Risch/LLL (48 lots sur Demeter) | **52 218 réussis / 55 ignorés / 4 xfail / 0 échec inattendu**, sur 52 277 tests collectés | Extras inclus ; collection et exécution exactement une fois ; compilateurs actifs |
| Suite complète du 28 octobre, LU/Risch/LLL (48 lots sur Demeter) | **51 948 réussis / 55 ignorés / 4 xfail / 0 échec inattendu**, sur 52 007 tests collectés | Extras inclus ; collection et exécution exactement une fois ; compilateurs actifs |
| Suite complète du 28 octobre, GJ sparse QQ (48 lots sur Demeter) | **51 606 réussis / 55 ignorés / 4 xfail / 0 échec inattendu**, sur 51 665 tests collectés | Extras inclus ; collection et exécution exactement une fois ; compilateurs actifs |
| Suite complète du 28 octobre, LLL/arrondis/protocoles (48 lots sur Demeter) | **51 366 réussis / 55 ignorés / 4 xfail / 0 échec inattendu**, sur 51 425 tests collectés | Extras inclus ; collection et exécution exactement une fois ; compilateurs Fortran et Rust actifs |
| Suite complète du 27 octobre, rationnel/AUTO/sparse (48 lots sur Demeter) | **50 416 réussis / 55 ignorés / 4 xfail / 0 échec inattendu**, sur 50 475 tests collectés | Extras inclus ; collection et exécution complètes vérifiées ; compilateurs Fortran et Rust actifs |
| Rust, ports Risch/LLL/HNF du 28 octobre | 762 unités + 1 intégration (651 cas LLL) + 5 doctests réussis ; 1 ignoré ; clippy core et bindings validés | Exécution fraîche complète ; 364 entrées natives, contrat Python LLL et TSV compilé compris |
| Rust, FFLU du 28 octobre | 755 tests unitaires + 1 intégration (651 cas LLL) + 5 doctests réussis ; 1 ignoré ; clippy core et bindings validés | Exécution fraîche complète ; 352 entrées natives vérifiées, TSV compilé compris |
| Rust, HNF/audit du 28 octobre | 751 tests unitaires + 1 test d’intégration (651 cas) + 5 doctests réussis ; 1 ignoré ; clippy core et bindings validés | Nouvelle exécution complète ; 350 entrées natives vérifiées, TSV compilé compris |
| Rust, solve/Risch/LLL du 28 octobre | 747 tests unitaires + 1 test d’intégration (651 cas) + 5 doctests réussis ; 1 ignoré ; clippy core et bindings validés | Nouvelle exécution complète ; 347 entrées natives, dont le TSV compilé ; aucun résultat natif réutilisé |
| Rust, LU/Risch/LLL du 28 octobre | 743 tests unitaires + 1 test d’intégration (651 cas) + 5 doctests réussis ; 1 ignoré ; clippy core et bindings validés | Nouvelle exécution complète ; 345 entrées natives dont le TSV compilé ; aucun résultat natif réutilisé |
| Rust, GJ sparse QQ du 28 octobre | 737 tests unitaires + 1 test d’intégration (651 cas) + 5 doctests réussis sur Demeter ; 1 ignoré ; clippy core et bindings validés | Nouvelle exécution complète ; 340 entrées natives vérifiées, dont le TSV compilé ; aucun résultat Rust réutilisé |
| Rust, LLL/arrondis/protocoles du 28 octobre | 732 tests unitaires + 1 test d’intégration (651 cas) + 5 doctests réussis sur Demeter ; 1 ignoré ; clippy core et bindings validés | Nouvelle exécution complète ; 340 entrées natives vérifiées, dont le TSV compilé ; aucun résultat Rust réutilisé |
| Rust, rationnel/AUTO/sparse du 27 octobre | 725 tests + 5 doctests réussis sur Demeter ; 1 ignoré ; clippy 1.99 core et bindings validés | Résultats du premier lot réutilisés après correction de deux fichiers Python ; mêmes 335 entrées natives et binaire vérifiés dans les deux snapshots |
| Suite complète, vague 17 (`f425cf0`, 24 octobre) | 44 338 tests Python ; `cargo test` et clippy validés | Historique, avant le dernier commit |
| Corpus EDO, rejeu du 26 octobre (`570f1832`) | 315 classifications / 281 indices / 263 résolutions par défaut identiques sur 315 ; 25 observations limitées en ressources | Politique stricte différente de la vague 17 ; écarts de forme et complétude encore à trier |
| Corpus d’intégration, R11 corrigé du 28 octobre (`e7b29759`) | **1 164 identiques / 36 différences revues / 10 doubles retenues**, sans erreur/timeout | 1 210 cas par moteur exactement une fois ; 2 420 appels ; tous les arbres natifs identiques à R2 ; 21 contrôles de paramètres et 10 E1 réussis |
| Corpus d'intégration, rejeu du 28 octobre (`c3f413bd`) | 1 164 identiques / 20 branches conditionnelles certifiées / 5 équivalences / 10 corrections E1 / 1 lacune affine nonzero / 10 doubles retenues | 1 210 cas exécutés exactement une fois par backend ; les 4 anciennes branches complexes sont corrigées ; 1 200 autres arbres natifs identiques à cb5 |
| Rejeu complet des exemples, 26 octobre (`576a9f76`) | 2 719 observations ordinaires identiques / 2 636 par mots-clés, sur 2 916 objets par mode ; 278 objets différents à trier | Toutes les entrées exécutées ; différences de forme, protocole ou capacité à distinguer |
| Correctifs solveur / `Integral`, dernier commit | 342 tests ciblés réussis sur l'extension release corrigée, selon les preuves du showcase | Validation ciblée ; pas une nouvelle suite complète |

Les anciens totaux « 1 363 différences », « 155 fonctionnalités manquantes »
et « 167 écarts d'impression » ne sont plus des compteurs du travail actuel.
Le [rejeu du 26 octobre](audits/sympy-1.14.0/documentation-2026-10-26/README.md)
fournit les observations actuelles sur `576a9f76` ; le tri mathématique des
278 objets reste ouvert. Voir l'[audit versionné](docs/sympy-parity-audit.md), le
[tri de la vague 8](audits/sympy-1.14.0/triage-wave8.md) et les
[preuves du showcase](examples/showcase_v4/data/evidence/native-validation.json).

- [x] Rejouer tous les exemples ordinaires et les appels par mots-clés sur
  la release identifiée `576a9f76` : [sorties complètes et couverture](audits/sympy-1.14.0/documentation-2026-10-26/README.md).
- [x] Isoler sur les sorties archivées 11 objets (32 observations) dont les
  différences sont uniquement des permutations d’ensembles affichés :
  [tri reproductible du 27 octobre](audits/sympy-1.14.0/documentation-triage-2026-10-27/README.md).
  Tous leurs exemples des deux modes sont vérifiés ; ce n’est pas un rejeu actuel.
- [ ] Trier les 267 autres objets du rejeu, avec preuves de domaine, précision
  et complétude ; actualiser les verdicts après les ports suivants.
- [x] Rejouer la suite complète après `5eb39d6` : build release, 48 lots sur
  Demeter, contrôles Rust et clippy ; sources, binaire, commandes et résultats
  liés dans les [preuves rationnel/AUTO/sparse du 27 octobre](audits/sympy-1.14.0/validation-2026-10-27-matrix-rational.json),
  les [preuves PolyMatrix/log_to_real du 27 octobre](audits/sympy-1.14.0/validation-2026-10-27-polymatrix-real-log.json),
  les [preuves grands ordres du 27 octobre](audits/sympy-1.14.0/validation-2026-10-27-literal-derivatives.json),
  les [preuves domaines inexacts du 27 octobre](audits/sympy-1.14.0/validation-2026-10-27-inexact.json),
  les [preuves LRT/limites du 27 octobre](audits/sympy-1.14.0/validation-2026-10-27-lrt-limits.json),
  les [preuves dérivées/algèbre du 27 octobre](audits/sympy-1.14.0/validation-2026-10-27-derivative.json),
  les [preuves calcul sparse/FF](audits/sympy-1.14.0/validation-2026-10-26-calculus.json),
  les [preuves de réduction rationnelle](audits/sympy-1.14.0/validation-2026-10-26-rational.json),
  les [preuves produits/matrices](audits/sympy-1.14.0/validation-2026-10-26-products.json),
  les [preuves stockage/solveurs](audits/sympy-1.14.0/validation-2026-10-26-storage.json),
  les [preuves des corrections natives](audits/sympy-1.14.0/validation-2026-10-26-corrections.json),
  les [preuves d’algèbre native](audits/sympy-1.14.0/validation-2026-10-26-algebra.json),
  les [preuves des noyaux natifs](audits/sympy-1.14.0/validation-2026-10-26-kernels.json),
  du [suivi natif](audits/sympy-1.14.0/validation-2026-10-26-native.json)
  et du [premier lot](audits/sympy-1.14.0/validation-2026-10-26.json).
- [x] Rejouer séparément les corpus EDO et d’intégration sur `570f1832`,
  avec référence indépendante et absence de SymPy dans les calculs natifs :
  [entrées, sorties et tri](audits/sympy-1.14.0/corpus-2026-10-26/README.md).
  Les formes EDO restent à trier ; le suivi des intégrales à paramètres est
  actualisé par le rejeu du 27 octobre ci-dessous.
- [x] Rejouer les 1 210 intégrations sur `cb5c0c98`, avec les domaines et
  substitutions des 20 nouvelles branches certifiées : [sorties et tri](audits/sympy-1.14.0/corpus-2026-10-27/README.md).
  La portée complexe de l'ancienne preuve Laplace est corrigée dans ce tri ;
  le correctif produit E1 est validé sur `daef4243`.
- [x] Rejouer les 1 210 intégrations sur `daef4243` : [nouvelle archive et preuves](audits/sympy-1.14.0/corpus-r2-2026-10-27/README.md).
  Les 4 anciennes lacunes complexes sont corrigées, les 1 200 autres arbres
  natifs restent identiques à cb5 ; 21 contrôles conditionnels et 10 E1 passent.
- [x] Rejouer les 1 210 intégrations sur la release `532386bd` :
  [corpus complet du lot PolyMatrix/log_to_real](audits/sympy-1.14.0/corpus-r5-2026-10-27/README.md).
  Les 1 210 arbres natifs sont identiques à R2 ; 21 contrôles conditionnels et
  10 E1 passent. La lacune affine nonzero et les doubles retenues restent suivies.
- [x] Rejouer les 1 210 intégrations sur la release `cc6d0210` :
  [corpus complet du lot rationnel/AUTO/sparse](audits/sympy-1.14.0/corpus-r6-2026-10-27/README.md).
  Les 1 210 arbres natifs sont identiques à R2 ; 21 contrôles conditionnels et
  10 E1 passent. La lacune affine nonzero et les doubles retenues restent suivies.
- [x] Rejouer les 1 210 intégrations sur la release `c3f413bd` :
  [corpus complet du lot LU/Risch/LLL](audits/sympy-1.14.0/corpus-r9-2026-10-28/README.md).
  Les 1 210 arbres natifs restent identiques à R2 ; 21 contrôles conditionnels
  et 10 E1 passent, sans erreur/timeout. Lacune affine et doubles retenues ouvertes.
- [x] Rejouer le corpus d’intégration sur le lot R10 corrigé : 1 210 entrées
  par moteur, 2 420 appels sans erreur ni timeout ; sorties natives identiques
  au précédent rejeu, 21 contrôles de paramètres et 10 E1 réussis. L’écart affine
  non nul et les 10 observations tenues restent ouverts.
  [Résultats complets](audits/sympy-1.14.0/corpus-r10-2026-10-28/README.md).
- [x] Rejouer les 1 210 intégrations sur le lot R11 corrigé : 2 420 appels
  exactement une fois, sans erreur/timeout, runners et limites inchangés ;
  1 210 arbres natifs identiques à R2, 21 contrôles de paramètres et 10 E1.
  La lacune affine non nulle et les 10 observations tenues restent ouvertes.
  [Archive complète](audits/sympy-1.14.0/corpus-r11-2026-10-28/README.md).
- [ ] Pour chaque écart ci-dessous, ajouter ou conserver une régression
  différentielle couvrant domaines, paramètres, options et échecs ; comparer
  les ensembles complets de solutions, les multiplicités et la précision réelle.

## Priorité : moteur Rust et performances vérifiées

Priorités actives : porter les algorithmes mathématiques restants en Rust,
fermer les lacunes de domaine/complétude/précision et valider chaque lot par
la suite complète sur Demeter. Les observations de coûts précédentes restent
archivées ; elles ne placent pas les microcoûts des gardes avant ces ports.

Le README exige un moteur de calcul Rust avec une API Python. L'arbre actuel
contient aussi des algorithmes mathématiques en Python : leur utilisation
d'expressions natives ne suffit pas à respecter cette architecture. Leur portage
doit conserver les domaines, replis, précisions et solutions existants.

**Règle de travail pour chaque lot (consigne du 29 octobre)** : le calcul va en
Rust (`crates/symox-core`, liaison PyO3 dans `crates/symox-py`). Le Python ne
sert qu'à l'admission, à l'aiguillage vers le natif, aux enveloppes d'API et au
passage des arguments. On n'ajoute ni n'étend d'algorithme Python : un écart de
parité qui ne se corrigerait qu'en Python est consigné comme lacune, pas corrigé
en Python. Faire passer les tests de SymPy n'est pas un but en soi si le calcul
reste en Python. Chaque rapport de lot indique, correctif par correctif, si le
calcul est en Rust ou en Python ; un lot qui n'a pas porté est renvoyé.

- [x] Inventorier les calculs Python effectivement utilisés, avec entrées
  publiques, appels Rust et boucles restantes : [audit du 27 octobre](audits/architecture/python-runtime-2026-10-27/README.md).
  Le comptage reproductible hors tests/docs donne 135 304 SLOC Python et
  156 489 Rust ; il ne prouve ni parité ni temps d'exécution.
  Le [recomptage R9](audits/architecture/python-runtime-2026-10-28-r9/README.md)
  donne 136 799 SLOC Python / 159 397 Rust, mêmes exclusions et script ;
  les algorithmes Python restants demeurent une dette explicite, pas une clôture.
- [ ] Porter les pilotes et noyaux Python restants vers Rust. Les premiers
  chantiers restants sont les préfixes froids, LU hors QQ et FFLU hors ZZ/QQ, les normalforms
  et les pilotes généraux. LLL stock et les arrondis exacts QQ/ZZ/RR sont
  intégrés dans les routes admises ; leurs replis complets restent conservés. La réduction stock
  `_solve_inexact` (Heurisch) est maintenant native et validée.

  Inventaire ciblé des calculs restant effectivement en Python (non exhaustif) :

  | Famille | Calcul Python restant | Source à porter |
  | --- | --- | --- |
  | Résolution | Simplexe hors route rationnelle native et protocoles personnalisés | `solvers/simplex.py` ; pivots et préparation `QQ` natifs validés |
  | EDO | Orchestration algorithmique, systèmes et séries | `solvers/ode/ode.py`, `systems.py`, `_series_api.py` |
  | Intégration | Pilote rationnel hors route QQ/ZZ admise, `log_to_real` général (degré ≥ 3/autres domaines), réduction hors QQ natif | `integrals/rationaltools.py` ; Horowitz, atan et LRT QQ natifs validés |
  | Intégration | Risch, équations différentielles paramétriques, Heurisch | `integrals/risch.py`, `prde.py`, `heurisch.py` |
  | Polynômes | `PolyMatrix` hors corps exacts audités, replis inexact/GF hors natif | `polys/polymatrix.py`, `polys/solvers.py`, `polys/_dense/` |
  | Matrices | Préfixes froids Berkowitz, produits/GJ hors routes natives, LU hors QQ et FFLU hors ZZ/QQ, substitutions hors QQ, LLL froid hors préfixe stock prouvé, préparations et normalforms | `polys/matrices/`, `matrices/_legacy_*.py` ; DDM Bareiss déjà validé |
  | Probabilités | Conditionnement, énumérations Markov, certains moments/assemblages | `stats/_joint.py`, `_processes.py`, `_queries.py` ; moments stock déjà natifs |
  | Tableaux et matrices | Dérivées d'ordre symbolique de matrices/tableaux, dérivées de traces imbriquées, repli `piecewise_fold` des constructeurs `Sum`/`Product`/`Integral` (vague 23, écrits en Python) | `tensor/array/array_derivatives.py`, `_public_calculus.py`, `_public_core.py` |

  Les chemins natifs spécialisés déjà portés ne ferment pas ces familles.
  Les classes, signatures, conversions et protocoles d'extension restent des
  adaptateurs Python ; les tests Python ne sont pas de la dette algorithmique.
- [x] Retirer des chemins mesurés les parcours symboliques et constructions de
  domaines Python inutiles : `radsimp`, polynôme caractéristique et inversion
  rationnelle ; joindre les mesures reproductibles avant/après et les preuves
  d'équivalence des résultats dans l'[audit de performance](audits/performance/2026-10-26-native-dispatch/README.md).
- [ ] Étendre ces mesures aux calculs réellement lents ; suivre séparément les
  appels initiaux, les appels répétés et les calculs encore non pris en charge.
  Le [contrôle du sélecteur rationnel](audits/performance/2026-10-27-rational-selector/README.md)
  conserve les régressions et leur amélioration. La [mesure du lot intégré](audits/performance/2026-10-27-rational-selector/integrated-r5/README.md)
  confirme LRT vector à 2,18–3,10× face à R2 ; le petit Atan linéaire vector
  reste à **+15,9 µs (+18,7 %)**, malgré son amélioration face à R3b.
  Un résultat faux ou non calculé ne constitue pas une accélération.
- [x] Porter les moments bruts binomiaux et normaux, ainsi que l'agrégation
  polynomiale de lois indépendantes exactes, en Rust. Les tables personnalisées,
  conditions, flottants et requêtes hors budget gardent leurs chemins complets.
- [x] Porter la réduction augmentée `_solve_inexact` en Rust, avec premier
  pivot non nul et opérations séparées à la précision du domaine sélectionné,
  y compris corps/anneaux à paramètres ; conversions et replis complets conservés.
- [x] Porter Gauss-Jordan flottant `RR`/`CC` en Rust : sélection des pivots,
  calcul à la précision du domaine, conventions d'entrée/sortie, classes et
  interruptions conservés, y compris les précisions 1/139/211 bits.
- [x] Porter la préparation cartésienne des coefficients d'Heurisch en Rust,
  avec racines principales exactes, refus des coefficients inconnus et retour
  de l'objet original quand aucun coefficient n'est traité.
- [x] Mettre en cache les coefficients Spouge par précision complète et
  paramètre `a` : 16 clés / 2 Mio maximum, calcul hors verrou et flux original
  au-delà du budget. Aucun cache approché des valeurs gamma.
- [x] Porter le comptage des opérations scalaires en Rust, avec multiplicité
  arbitraire des DAG et maintien des protocoles de repli.
- [x] Certifier en Rust la couverture des racines rationnelles déjà calculées,
  pour éviter une seconde résolution ; une preuve déclinée conserve le solveur.
- [x] Porter les moments de Poisson et les pivots du simplexe rationnel en
  Rust. La préparation des contraintes/bornes/variables QQ est aussi native ;
  les coefficients non rationnels et protocoles personnalisés gardent leurs
  parcours complets. Mesures et certificats dans le suivi de performance.
- [x] Corriger les facteurs répétés tenus et contrôler les exposants stockés
  sur 32 bits sans perdre les annulations exactes ni les replis EX disponibles.
  Ce contrôle est borné aux chemins audités, pas une garantie universelle.
- [x] Corriger les degrés totaux, homogénéisations et frontières publiques
  `contains` / `reduced` auditées : sommes larges, précision et refus explicite
  d’un véritable dépassement, sans confusion avec zéro ou non-appartenance.
- [x] Introduire le stockage natif Dense/Sparse EX/CC : construction, vues,
  conversions, structure et opérations additives/scalaires, sans densifier les
  grandes coordonnées ni perdre ground/précision ou les anciens parcours.
- [x] Étendre ce stockage aux produits/carrés/puissances : regroupements et
  arrondis récursifs DMP conservés, y compris la frontière de Karatsuba EX.
- [x] Étendre ce stockage au calcul différentiel/intégral EX/CC, avec
  facteurs BigInt, regroupement des specs et précision du domaine ;
  annulations exactes avant la borne de coordonnées.
- [ ] Étendre les autres algorithmes encore denses et l’évaluation non nulle.
- [x] Raccorder l’évaluation EX en zéro au chemin public univarié :
  **186,2 → 68,8 µs, 2,7×** sur le témoin exact contrôlé ; les surcharges
  publiques et autres contextes conservent le parcours complet.
- [x] Conserver les payloads Float0 EX intérieurs/trailing dans construction,
  TC/nth, évaluation et pickle. Les zéros de tête entièrement nuls restent
  éliminés selon les conventions de référence.
- [ ] Étendre ce contrôle aux conversions implicites RR/CC vers EX, à
  l'inférence implicite et aux anciens constructeurs bruts hors route native.
- [x] Porter la multiplication DDM entière/rationnelle et Gauss-Jordan QQ ;
  [mesures et contrôles exacts](audits/performance/2026-10-26-native-products/README.md).
- [x] Porter la réduction fraction-free dense/sparse utilisée par DomainMatrix
  AUTO en ZZ/QQ ; routes FF/CD et conventions de dénominateurs préservées.
- [x] Porter la suppression des dénominateurs QQ globale et par ligne ;
  le témoin public AUTO/CD gagne **1,25×** sans changer ses sorties.
- [x] Porter le déterminant Bareiss DDM ZZ/QQ, avec mutations, alias, erreurs
  et replis conservés ; `Matrix.det()` était déjà natif.
- [x] Porter Berkowitz dense ZZ/QQ, avec coefficients exacts, conventions de
  sortie, signaux et protocoles conservés ; validé en suite complète sur `daef4243`.
  [Mesures `DDM.charpoly()`](audits/performance/2026-10-27-native-berkowitz/README.md) :
  **5,2× à 13,7×** sur 15 matrices, coefficients et Cayley-Hamilton certifiés.
- [x] Porter RREF/nullspace ground PolyMatrix sur QQ et les corps QQ(y)/ZZ(y)
  audités, utilisés par PRDE/Risch, avec les conventions d'élimination propres
  à chaque entrée, domaines, pivots et replis conservés.
- [x] Porter log_to_real pour les polynômes auxiliaires rationnels squarefree
  de degré 1/2 et représentants exacts admis ; racines, normes et atan exacts,
  sans preuve evalf ni changement des branches. Les autres cas restent complets.
  [Preuves et suite complète](audits/sympy-1.14.0/native-polymatrix-real-log-2026-10-27/README.md).
- [x] Porter le pilote rationnel QQ/ZZ admis : extraction/annulation, parties
  polynomiale et rationnelle, LRT puis logarithmes réels certifiés, avec les
  mêmes unités, représentants et constantes de branche. Toute preuve déclinée
  reprend le pilote complet ; aucun résultat partiel n’est publié.
- [ ] Étendre ces routes aux autres domaines/formes, RootSum/degrés supérieurs
  et options générales ; porter les pilotes Risch/PRDE encore Python.
- [x] Porter les statistiques et décisions AUTO stock ZZ/QQ : densité,
  longueurs binaires, PPCM et seuils binary64 identiques à CPython ; mêmes
  choix GJ/FF/CD, formats et dénominateurs, replis complets conservés.
  [Mesures publiques du lot intégré](audits/performance/2026-10-27-native-auto/README.md) :
  QQ/CD dense 20×20 **1,36×**, avec contrôle sur le même binaire et 2 016
  certificats exacts ; les petits 3×3 gardent **+2,3 à +3,0 %** de surcoût.
- [x] Porter le noyau Berkowitz sparse ZZ/QQ, sans densification, avec
  suffixes, produits matrice-vecteur et convolution Toeplitz exacts.
  Le préfixe froid réellement observé est conservé avant la conversion qui
  autorise le natif ; une récursion froide profonde peut encore rester Python.
  [Preuves et suite complète](audits/sympy-1.14.0/native-matrix-rational-2026-10-27/README.md).
- [x] Raccorder LLL entier exact aux entrées publiques DDM/DomainMatrix admises :
  calcul Rust, certificats T*A=B / déterminant / Gram-Schmidt-Lovasz, conventions
  de publication et erreurs conservées. Corriger floor/ceil exact de QQ/ZZ et
  des dyadiques RR stock ; aucun passage binary64 sur ces routes.
  [Preuves et suite complète](audits/sympy-1.14.0/native-lll-rounding-protocols-2026-10-28/README.md).
  Les neuf écarts de la référence LLL dus à l'ancien arrondi restent archivés ;
  aucune accélération n'est comptée contre un résultat faux ou une assertion.
- [x] Vérifier nativement les snapshots de fonctions des outils rationnels,
  avec cellules réelles, interruptions et reprise du vérificateur Python si
  les dépendances changent ; bootstrap all/len d'un module étranger refusé.
- [x] Porter Gauss-Jordan sparse stock QQ (`sdm_irref`, SDM/DomainMatrix GJ/AUTO)
  en Rust : pivots, annulations et arithmétique exacte ; les trois résultats,
  ordre des métadonnées et identités de clés sont conservés. Les autres
  domaines, zéros explicites et protocoles modifiés gardent le parcours complet.
  Le réimport d’un module neuf décline sa preuve plutôt que d’utiliser les
  snapshots d’un module retiré ; reload ordinaire et module d’origine préservés.
  [Preuves et suite complète](audits/sympy-1.14.0/native-sdm-irref-2026-10-28/README.md).
  [Mesures du lot intégré](audits/performance/2026-10-28-native-sdm-irref/README.md) :
  témoins 24 lignes **2,26–2,70× au premier appel**, **2,67–3,34× répétés** ;
  petits 2×3 encore **+17,3/+20,5 µs (+28,3/+26,9 %)** au premier appel.
- [x] Porter LU stock QQ de `ddm_ilu` en Rust, utilisée par DDM/DomainMatrix :
  préfixes de dimensions réellement exécutés, dépendances effectives range/min,
  publication et permutations étape par étape, alias et erreurs conservés.
  Les autres domaines/protocoles et lignes aliasées gardent leur corps complet.
- [x] Porter la boucle Risch `polynomial_reduce` QQ/ZZ admise, avec tours de
  dérivations polynomiales rationnelles, promotions de domaines q/r et certificat
  p=Dq+r. L’admission degré≤256/bits≤4096 ne réduit pas le pilote complet,
  repris pour les autres cas. Ce n’est pas un port Risch/PRDE général.
  [Preuves et suite complète](audits/sympy-1.14.0/native-lu-risch-lll-2026-10-28/README.md).
- [x] Refuser LLL froid sans parcourir la preuve lourde quand le cache ABC
  est réellement froid : lecture seule, aucune classification ni préchauffage,
  aucun appel d’un getter/alias étranger capturé avant import.
- [x] Porter les substitutions triangulaires QQ de `ddm_ilu_solve`, avec RHS
  multiples, formes/permutation/allocation y réelles, publications y/x et erreurs
  partielles conservées. Getters du frame natifs, sans résolution prématurée
  des callbacks tardifs ; autres domaines et alias non prouvés gardent le corps complet.
- [x] Refuser rapidement les candidats Risch sans étape par lecture seule
  du payload/degré et du vrai membre C d. Ce refus ne produit aucun q/r :
  l’ancien corps reste entier, y compris getters modifiés et erreurs de constructeur.
- [x] Lire les métadonnées LLL sans encoder les clés non consommées : toutes
  les clés str exactes validées avant un seul accès utile ; pas de callback
  d’égalité étranger, pas de classification ni de mutation/réutilisation du cache.
  [Preuves et suite complète](audits/sympy-1.14.0/native-solve-risch-lll-2026-10-28/README.md).
- [x] Supprimer les lectures Python de métadonnées ajoutant des événements d’audit
  dans les admissions LU/solve/HNF : cinq champs relus par APIs C publiques
  optionnelles sur CPython 3.10.12/.20 Linux x86_64, repli complet ailleurs.
  Types et toutes les clés exactes sont vérifiés avant les accès utiles :
  modules/helpers/domaines et clés étrangères n’ajoutent aucun callback.
  551 ciblés et 6 276 couplés locaux/Demeter, avant le gate RELEASE complet.
- [x] Préserver les protocoles de fermeture et le transfert scalaire officiel
  diffgeom : tuple exact avant lecture, méthodes/origines/payloads authentifiés,
  collisions et objets étrangers refusés avant callbacks ajoutés. Les Wild
  à propriétés appellent le calcul complet d’origine ; nouveau RELEASE et
  gate complet corrigé distincts du premier gate R11 échoué.
- [x] Porter HNF ordinaire ZZ en Rust : gcd, divisions exactes, signes et
  réductions de colonnes après le vrai préfixe de copie. Certificat unimodulaire
  A*V=[0|H], pivots/restes/ordre HNF, identités et publications conservés.
  [Code et preuves complètes](audits/sympy-1.14.0/native-normalforms-risch-lll-2026-10-28/README.md).
- [x] Vérifier les builtins effectivement conservés par les fonctions HNF,
  avant capture/admission : restaurer globals.__builtins__ ne change pas
  f_builtins déjà figé par FunctionType. Les fonctions étrangères reprennent
  le calcul entier avec ses callbacks/erreurs ; 4 contre-exemples avant
  correction et 5 contrôles après sont conservés, avec 562 ciblés par hôte.
- [x] Authentifier aussi le code mathématique HNF avant un nouveau bootstrap
  du garde : `_gcdex` et `add_columns` sont comparés aux templates locaux,
  avec builtins réels ; les fonctions étrangères reprennent leurs callbacks/erreurs.
  Les quatre échecs R12 sur les deux hôtes restent archivés ; un cinquième cas
  de fausse autorité appartient seulement au candidat privé rejeté. Le lot R13
  corrigé et sa suite complète couvrent les cinq régressions sans dépendance
  sémantique à une table de codes fournie par un autre module.
  [Contre-exemples R12 gelés](audits/sympy-1.14.0/hnf-bootstrap-gap-2026-10-28/REPORT.json).
- [x] Porter les quatre phases mathématiques FFLU stock ZZ/QQ en Rust :
  élimination Bareiss, publication de permutation, extraction L/U et produits D,
  avec copies, allocations, alias et conversions DDM/DomainMatrix conservés.
  Corriger aussi le calcul générique : recherche de pivots sur toutes les colonnes,
  diviseur issu de la colonne réelle du pivot précédent. Deux défauts SymPy 1.14
  sont corrigés intentionnellement, avec anciennes sorties/erreurs archivées.
  404 cas de domaines, 1 616 appels certifiés par parcours natif et par repli
  complet sur chacun des deux hôtes ; certificat Fraction PA=LD^-1U,
  rang/pivots complets, facteurs inférieurs/diagonaux inversibles et échelon.
  [Port, sources et nouvelle suite complète](audits/sympy-1.14.0/native-fflu-2026-10-28/README.md).
- [x] Porter `splitfactor` Risch QQ fermé à dérivation stock : gcd exact,
  dérivées, quotients exacts et extraction récursive des facteurs spéciaux.
  Unités, représentants/domaines et alias conservés à chaque niveau ; admission
  degré≤64/bits≤4096, ancien corps entier repris pour tous les autres cas.
  La lacune Hermite préexistante sur les représentants ZZ/QQ reste ouverte ;
  ceci ne clôt pas Risch/PRDE ni l’intégration générale.
- [x] Sélectionner le calcul LLL Rust dès le véritable premier appel froid
  admis : le préfixe réellement consommé garde son cache, ses effets ABC,
  ses erreurs et ses audits. Témoins séparés froid/refusé et froid/natif,
  1 104 triplets publics et 11 tests publics de premier appel. Admission froide
  vérifiée/activée sur CPython 3.10.12 et 3.10.20, Linux x86_64 ; ailleurs
  et pour les protocoles non prouvés, le calcul complet est repris.
- [x] Porter HNF modulo D stock ZZ en Rust : récurrence exacte, publications
  et provenance du module D conservées ; préfixes publics de domaines, rang,
  erreurs, allocation et conversion exécutés selon le parcours d’origine.
  Les multiples déterminantaux valides ont un certificat de réseau unimodulaire ;
  les autres options/protocoles déclinent vers l’ancien algorithme entier.
  [Trois ports, preuves et suite complète R13](audits/sympy-1.14.0/native-math-ports-2026-10-28/README.md).
- [ ] Porter les préfixes froids/préparations restants, GJ hors QQ admis,
  LU hors QQ et FFLU hors ZZ/QQ, substitutions hors QQ admis, HNF numberfield
  et formes de Smith sans réduire
  les domaines ni les callbacks.
  Les routes LLL froides hors préfixe stock prouvé gardent le calcul entier.
  Les [mesures historiques LLL](audits/performance/2026-10-28-native-lll/README.md)
  restent archivées : appels répétés 4/6 lignes **4,20–5,60×**, premier appel
  alors Python et **+1,9 à +21,1 % (66–287 µs)**. Aucun nouveau gain chronométré
  n’est déduit du port mathématique R13 ; les anciennes observations sont conservées.
  [Mesures contrôlées du 27 octobre](audits/performance/2026-10-27-native-derivative/README.md).

- [x] Porter la récurrence des moments gaussiens corrélés vers un graphe
  natif sans récursion Python, avec coefficients/arrondis conservés et
  contrôles des protocoles. Ordre 1200 et import de `diffgeom` couverts ;
  les empilements Vector/Tensor stock sont reconnus sans masquer les vrais
  callbacks. Le noyau conserve le repli pour les autres protocoles.
- [x] Corriger les preuves LLL de réimport/bootstraps : nouveau namespace
  `mfloor`, lookups tels que `max`, alias de sélection et getters ABC étrangers
  déclinent avant tout appel non consommé par l’ancien calcul. Valeurs, erreurs
  et callbacks sont comparés au parcours complet ; le calcul froid reste Python.
- [x] Corriger l’erreur Unicode du probe LLL sur des attributs non consommés
  contenant un surrogate, y compris métadonnées du module ABC. Treize échecs
  du précédent lot et sept contrôles sont conservés, sans perte des anciens témoins.
- [x] Supprimer le gros surcoût Risch sans étape du lot R9 : le rejeu à trois
  versions R8/R9/R10 retrouve ~55–188 µs répétés au lieu de ~772–958 µs.
  Des écarts résiduels face à R8 restent visibles (~+1,5–6,8 µs sur trois contrôles).
  [Mesures et certificats exacts](audits/performance/2026-10-28-native-risch-r10/README.md).
- [ ] Réduire les surcoûts LLL répétés : le correctif R10 varie encore de
  ~3,6 % plus lent à ~7 % plus rapide que R9 selon le cas ; pas de récupération
  générale annoncée. [Tous les cas](audits/performance/2026-10-28-native-lll-r10/README.md).
- [x] Mesurer séparément HNF, LU et solve entre R10 corrigé et R11 corrigé :
  3 456 certificats exacts, 432 processus frais, trois campagnes séquentielles
  après les suites complètes et arrêt des autres calculs des agents.
  [HNF](audits/performance/2026-10-28-native-hnf-r11/README.md) : rapports R10/R11 0.642–15.471× initiaux, 0.553–18.482× répétés ; 13 médianes plus lentes conservées.
  [LU](audits/performance/2026-10-28-native-lu-r11/README.md) : rapports R10/R11 0.576–0.819× initiaux, 0.347–0.831× répétés ; 36 médianes plus lentes conservées.
  [Solve QQ](audits/performance/2026-10-28-native-solve-r11/README.md) : rapports R10/R11 0.593–2.012× initiaux, 0.376–0.752× répétés ; 33 médianes plus lentes conservées.
  Un rapport supérieur à 1 indique R11 plus rapide ; tous les coûts absolus,
  régressions et commandes restent publiés. Aucun gain global n’est déduit.
- [ ] Réduire les surcoûts observés LU/solve du lot R11 en conservant les
  audits, callbacks et replis complets : les 18 cas LU sont plus lents
  au premier appel (+42,5 à +64,0 µs ; +22,1 à +73,6 %) et en répété
  (+31,5 à +45,5 µs ; +20,4 à +188,4 %).
  Les 18 cas solve répétés restent aussi plus lents (+73,8 à +92,2 µs ;
  +33,0 à +165,8 %), ainsi que 15 premiers appels sur 18 (+75,5 à
  +96,7 µs ; +28,8 à +68,6 %).
  Les 13 médianes HNF plus lentes restent publiées avec les gains des autres cas.
- [ ] Réduire les premiers appels LU/solve encore plus chers sur certains cas.
  La substitution QQ native gagne jusqu’à ~4,9× en répété mais deux premiers
  appels sparse restent plus lents (~+7 et +74 µs).
  [Solve QQ](audits/performance/2026-10-28-native-solve-r10/README.md).
  Les observations précédentes R9 restent archivées séparément :
  [LU](audits/performance/2026-10-28-native-lu-r9/README.md),
  [LLL](audits/performance/2026-10-28-native-lll-r9/README.md),
  [Risch](audits/performance/2026-10-28-native-risch-r9/README.md).
- [ ] Optimiser les premiers appels gamma et les nouveaux niveaux de précision
  au-delà de ce cache. Le raccourci Gauss convergent retiré reste désactivé ;
  toute nouvelle route doit conserver précision et performances mesurées.

## Acquis à conserver

Ces cases ferment les anciens chantiers dans la portée des tests et des
sondes documentés ; elles ne déclarent pas une couverture exhaustive.

- [x] Corps algébriques généraux, éléments primitifs, plongements et isolation
  certifiée des racines complexes ; opérations polynomiales avec coefficients
  algébriques et `CRootOf` dans les cas audités.
- [x] Réduction exacte avec objets racines, extraction et regroupement des
  radicaux ; `RootSum` (`evalf`, `N`, `latex`, `subs`) et conversions
  d'`AlgebraicNumber`. L'ancienne mention « Non commencé » est obsolète.
- [x] Anneaux imbriqués, fractions gaussiennes, `DomainMatrix`/`DDM`/`SDM`,
  polynômes sur Z/nZ ; groupes de Galois de degrés 5 et 6.
- [x] Domaines `RR`/`CC` explicites et inférés avec générateurs, précision
  conservée, PGCD multivarié sur `RR[y]`, `sqf_list` gaussien multivarié.
- [x] Racines numériques proches et formes factorisées sur `EX` dans les cas
  du [corpus dédié](tests/test_nroots_precision_parity.py).
- [x] Intégration Meijer G, options Risch/heurisch/manual, `NonElementaryIntegral`,
  transformées et `AccumBounds` ; conditions de convergence et de support
  corrigées dans les cas audités.
- [x] Familles complexes de solutions exp/trigo, exponentielles imbriquées,
  coefficients trigonométriques symboliques, Lambert, radicaux, équations
  modulaires et `solve_triangulated(extension=True)` dans les familles testées.
- [x] `solve(sec(2*x) - 3)` et `solve(csc(3*x) + 4)` ; le blocage du
  `nonlinsolve` rationnel à racines cubiques est corrigé.
- [x] `limit_seq` avec `dominant`, sommes imbriquées et déterminisme ;
  asymptotiques gamma/Bessel et nombreuses limites sur les coupures.
- [x] Modules d'EDO/EDP et systèmes différentiels de référence, Frobenius,
  options auditées et replis natifs conservés. L'ancien blocage de 40 minutes
  de l'heuristique de Lie n'est plus l'état courant.
- [x] `evalf(strict=True, maxn=...)`, fonctions spéciales multiprécision,
  grands exposants binaires, précision près des pôles, choix des coupures
  de Meijer G et propagation immédiate des interruptions.
- [x] Séquence de simplification de référence, `powdenest(polar=True)`,
  opérations non commutatives et corrections de `Piecewise` dans les cas audités.
- [x] Matrices flottantes et symboliques : routes de Jordan, inversion
  `RR`/`CC`, tests de pivots, inférence des domaines, dérivées matricielles
  et tableaux ; gardes explicites des opérations non commutatives non supportées.
- [x] Budgets de ressources audités en vague 17 : corrections des faux
  résultats à épuisement de budget, dont Wang, séries, `Mod`, `checkodesol`
  et `decompogen`. Les limites de capacité restent suivies ci-dessous.
- [x] Dernier commit : signes des candidats algébriques par encadrements
  certifiés, préimages quadratiques exactes sans explosion polaire et
  protocole de simplification d'`Integral` conservant l'intégrale retenue.

## 1. Corps algébriques, polynômes et domaines

- [x] Inférer `CC` sans générateurs explicites : `Poly(x + 0.5*I)` garde
  seulement `x` comme générateur et conserve la précision des coefficients.
- [x] Convertir les coefficients de `CC` vers `QQ_I`, y compris les domaines
  inférés : la conversion explicite suit la convention à 53 bits de SymPy.
  Les calculs dans un `CC` de précision supérieure gardent leur précision.
- [x] Faire passer `cancel(..., domain=RR/CC)` par la route inexacte :
  expressions, tuples, options et domaines composites dans les familles testées.
  Corriger aussi les pertes de précision des coefficients symboliques `RR` et
  de la reconstruction exacte utilisée par l’arithmétique `CC`.
- [x] Étendre `cancel`, PGCD, cofactors et PPCM aux domaines `CC` et composites
  complexes à plusieurs générateurs, via les algorithmes gaussiens Rust ;
  contenus, tuples, monômes, zéros et précision 140 bits couverts par les
  [régressions multivariées](tests/test_complex_multivariate_cancel_parity.py).
- [x] Étendre annulation, PGCD, cofactors et PPCM multivariés aux corps de
  fractions réels `RR(y)`, avec contenus, unités et précision conservés dans
  les familles testées. Les normalisations résiduelles restent suivies ci-dessous.
- [ ] Aligner la normalisation PGCD/PPCM et `sqf_list` sur `RR[y]`, ainsi que
  le domaine de sortie d'`integrate` (`RR[y]` / `RR(y)`). Les cas documentés
  de normalisation sont égaux à une unité près ou ont les mêmes coefficients.
- [ ] Requalifier les derniers écarts de réduction / `nsimplify` avec radicaux
  et `CRootOf`, ainsi que les conversions et erreurs de domaines résiduelles,
  après les vagues `algfield`, `algroot`, `rrcc` et `ccdom`.

Preuves : [algfield](docs/algorithmic-parity.md#algebraic-fields-root-objects-and-inexact-domains-algfield-2026-10-23),
[algroot](docs/algorithmic-parity.md#algroot-wave-16-root-objects-and-algebraic-fields),
[ccdom](docs/algorithmic-parity.md#ccdom-wave-17-explicit-cc-rry-gcd-gaussian-composites).
La factorisation sur l'extension `(-2)**(1/12)` est corrigée par `reslimits` :
ne pas rouvrir cet ancien écart de `ccdom`.

## 2. Intégration symbolique, transformées et distributions

- [x] Corriger la primitive de `1/(1 + 3*cos(x)**2)` qui sautait de `-pi`
  et perdait la croissance par période depuis `5eb39d6`. Heurisch retrouve
  les paires conjuguées à coefficients numériques ; les tests vérifient
  dérivée exacte, accroissement par période et continuité, sur trois familles.

- [x] Valider en suite complète le port Rust de `ratint_ratpart` en QQ/ZZ ;
  34 tests ciblés et 96 certificats différentiels exacts réussis. Les autres
  domaines et tailles hors route native gardent le calcul complet existant.
  [Mesures avant/après](audits/performance/2026-10-26-native-rational/README.md) :
  gains de 4,5× à 20,6× sur les cas rationnels contrôlés.
- [x] Porter `log_to_atan` QQ/ZZ en conservant les représentants euclidiens
  et les constantes de branches ; les wrappers holonomes stock sont reconnus
  sans ignorer les vraies surcharges. Environ **2×** sur les cas Euclide mesurés.
- [x] Porter et valider `ratint_logpart` LRT QQ/ZZ : représentants, domaines,
  signes, wrappers stock et vrais callbacks couverts par la suite complète.
  [Mesures release appariées](audits/performance/2026-10-27-native-lrt/README.md) :
  **2,26× à 3,23×** sur cinq familles, 144 comparaisons de paires exactes.
- [ ] Étendre LRT aux domaines restants sans retirer leurs calculs complets.
- [x] Donner des branches certifiées à 20 des 21 anciennes intégrales à
  paramètres retenues, avec l'intégrale d'origine en défaut et substitutions
  partielles sûres. Les options/méthodes explicites gardent leurs parcours.
- [ ] Compléter la pente affine supposée non nulle et les autres familles
  non décidées ; préserver la spécialisation zéro et les pôles réels.
- [x] Corriger les branches complexes des règles Laplace réciproques et
  logarithmes affines déjà gardées : `exp(s)*expint(1,s)` donne la continuation
  correcte de `1/(t+1)` sur `Re(s)>0`. Preuve analytique, décalages et moments
  contrôlés ; autres branches générales encore à compléter.
- [ ] Étendre la couverture résiduelle de Risch/heurisch et des transformées,
  en distinguant primitive équivalente, conditions différentes et calcul absent.
- [ ] Réduire la croissance des expressions de `ratint(1/(x**3 + x + 1))` :
  relevé du 26 octobre sur `146f6917`, **82 704 caractères contre 16 653**
  pour SymPy. L'écart persiste ; [mesure identifiée](audits/sympy-1.14.0/ratint-growth-2026-10-26/ratint-growth-native.json).
  La taille du texte seule ne prouve ni erreur mathématique ni équivalence.
- [ ] Compléter les calculs dépendant des limites encore non décidées, dont
  `integrate(atan(a*x)/x**2, (x, 1, oo))`, avec leurs conditions paramétriques.
- [ ] Étendre les distributions à paramètres, ordres dérivés et supports plus
  généraux lorsque les conditions de validité sont établies.
- [ ] Documenter et tester systématiquement la convention de masse d'un delta
  en borne finie (demi-masse dans les cas relevés, contre 0 chez SymPy), et
  l'écart d'erreur de `series(Heaviside(x - 1), x, oo)`.

Le [rejeu du 27 octobre](audits/sympy-1.14.0/corpus-2026-10-27/README.md)
certifie 20 branches conditionnelles et conserve une lacune affine nonzero.
Les dix anciens écarts Laplace se répartissent en six corrections réelles de
préfacteur et quatre lacunes complexes natives encore à corriger. Recopier
la formule de référence, elle-même fausse, ne les résout pas. Conserver les
intégrales retenues quand les pôles ou la convergence ne sont pas décidés.
Les produits de distributions tels que
`DiracDelta(x)**2` restent retenus, leur produit n'étant pas défini en général.

Preuves : [intrest4](docs/algorithmic-parity.md#intrest4-wave-13-integration-corpus-leftovers),
[defparam](docs/algorithmic-parity.md#defparam-wave-17) et
[tests des distributions](tests/test_defparam_parity.py).

## 3. Résolution, paramètres et complétude des solutions

- [x] Conserver toutes les préimages lorsque l’inverse trigonométrique est
  exact, dont `solve(sec(2*x) + 1, x) = [-pi/2, pi/2]` et `csc(2*x) - 2`.
- [x] Étendre les multiples entiers de `sec`/`csc`, dont `sec(4*x) - 3`,
  en gardant toutes les racines de l’unité et les racines doubles. Le filtrage
  réel et les signes de logarithmes utilisent des encadrements certifiés.
- [x] Étendre `solve(sec/csc(...))` aux arguments affines décalés et aux
  pentes rationnelles dans la portée native. Les `2*p` préimages (ou `p`
  racines doubles) sont les représentants finis de `solve` ; les familles
  périodiques complètes relèvent de `solveset`. Pentes symboliques et tailles
  non représentables gardent leurs replis, sans liste partielle déclarée complète.
- [x] Inverser les chaînes exponentielles affines à second membre symbolique
  dans Rust : toutes les branches entières, préimages nulles exclues,
  coefficients nuls et conditions d'origine conservés. Voir les
  [régressions de familles complètes](tests/test_symbolic_nested_exp_parity.py).
- [x] Préserver équations et gardes lors des substitutions partielles
  sous garde explicite de non-finitude : `exp(exp(x)+a)=y`, `y=0` puis
  `a=-oo`, permutations de plusieurs paramètres et domaines étendus. Le
  strict xfail ciblé est retiré ; NaN conserve une équation tenue.
- [ ] Étendre cette protection aux substitutions algébriques et gardes
  arbitraires au-delà du contexte conjonctif explicitement audité.
- [ ] Compléter les cas de paramètres dégénérés et les composantes singulières
  des systèmes. La résolution générique est déjà portée ; elle ne décrit
  pas automatiquement toutes les spécialisations des coefficients.
- [ ] Étendre les familles transcendantes et les représentations en ensembles
  au-delà des familles auditées, avec restrictions réelles et intersections.
- [ ] Poursuivre l'audit des options de `solve`, `solveset` et `roots`, et la
  vérification exacte des candidats avec branches, pôles et exclusions.
  Le signe algébrique fermé est maintenant certifié ; un cas indécidable
  doit rester conditionnel ou non pris en charge.
- [ ] Compléter les inégalités, tests d'inclusion et calculs d'analyse réelle
  encore indécidables là où la référence sait conclure.

Preuves : [holo-rest](docs/algorithmic-parity.md#holonomic-multiple-derivatives-triangulated-extension-and-parametric-systems-holo-rest-2026-10-23),
[solve/EDO](docs/algorithmic-parity.md#solve-and-ode-leftovers-ode-solve-2026-10-23)
et [tests de vérification des candidats](tests/test_solve_checks.py).
Les anciens écarts Lambert, radicaux produits, pôles symboliques et racines
cubiques trigonométriques ont des correctifs ultérieurs : les réévaluer avec
leurs tests avant de les remettre dans cette liste.

## 4. Limites, séries et équations différentielles

- [x] Décider la limite de l'issue 22334 :
  `(n + 1)**k / ((n + 1)**(k + 1) - n**(k + 1))` pour `n -> oo`, `k > 0`
  donne `1/(k + 1)`. Assertion initiale conservée, xfail retiré, suite complète verte.
  La normalisation passe d'abord par une échelle neutre et un certificat rationnel
  lorsque l'argument tourne ; le faux fini trouvé en revue reste rejeté.
- [ ] Étendre les limites symboliques encore retenues, notamment certains
  rapports de différences décalées et racines ; une réponse indécidée reste explicite.
- [ ] Étendre les développements à l'infini de `floor`, `ceiling` et `Mod`
  avec coefficients `AccumBounds`, sans annulation abusive de `AB - AB`.
- [ ] Compléter les limites encore retenues : annulations dominantes de
  `Ei`/`erf`, `uppergamma(x, 1)/gamma(x)`, `subfactorial`, reste de Stirling
  de `log(factorial(x))`, `harmonic(n, m + 1)` et certaines limites bilatérales.
- [ ] Compléter les séries d'une fonction appliquée à un `Piecewise` et les
  séries en l'infini avec conditions indécidables ; localiser correctement
  le nouveau point limite dans `O(x**2).subs(x, cos(x))`.
- [ ] Réconcilier les derniers indices et formes du corpus EDO avec les
  correctifs ultérieurs des solveurs et de `simplify`. Les solutions
  séparables/Riccati auparavant implicites sont maintenant résolues ;
  les différences de forme seules ne sont pas des solutions manquantes.

Preuves : [limrest2](docs/algorithmic-parity.md#limrest2-wave-14),
[limite 22334](tests/test_limwrong_parity.py),
[Piecewise](docs/algorithmic-parity.md#pw-nc-wave-15-piecewise-and-noncommutative-leftovers)
et [corpus EDO](docs/algorithmic-parity.md#ode-corpus-forms-and-errors-ode-forms-2026-10-23).
Le `limit_seq` imbriqué donnant 4 et l'étape `dominant` sont déjà portés et
couverts par [leurs régressions](tests/test_limitseq_dominant_parity.py).

## 5. Précision et évaluation numérique

- [x] Arrondir exactement les éléments stock QQ/ZZ et les littéraux finis RR
  selon leur valeur stockée, y compris grandes mantisses et exposants :
  floor/ceil ne passent plus par binary64. Surcharges float/index et erreurs
  des autres routes conservées ; la capacité entière reste explicitement bornée.
  [Preuves et suite complète](audits/sympy-1.14.0/native-lll-rounding-protocols-2026-10-28/README.md).
- [x] Corriger les arrondis de multiplication/division CC littéraux finis :
  cellules d’arrondi certifiées par comparaisons dyadiques exactes, y compris
  grands écarts d’exposants et précisions 1/53/139/140/211 bits. Les témoins
  140/211 bits révèlent aussi un produit SymPy mal arrondi, qui n’est pas copié.
  [Preuves et validation complète](audits/sympy-1.14.0/native-inexact-2026-10-27/README.md).
- [x] Rendre réflexive l’égalité des composites RR/CC typés du même domaine,
  sans modifier hachage, coefficients ni précision ; `(a/one)==a` est rétabli.
- [ ] Réconcilier les wrappers Python de ring et les éléments typés natifs,
  ainsi que l’ordre non numérique ; l’égalité entre représentations distinctes
  et l’évaluation numérique générale restent des chantiers séparés.
- [x] Évaluer `N(hyper([1, 1], [2], 1.0))` en `oo`, et certifier la divergence
  des `pF_(p-1)` à paramètres rationnels positifs, excès nul et argument 1.
  Simplifier `uppergamma(-k, 0.0)` en `zoo` pour `k` entier strictement positif.
  Les ordres indécidables, les coupures et les interruptions sont conservés.
- [x] Certifier nativement les pôles `pF_(p-1)(1)` à paramètres rationnels,
  excès non positif et signe gamma exact, dont `hyper([1/3], [], 1)` ;
  corriger les faux finis et propager les pôles dans les compositions numériques.
  Paramètres terminants/interdits et arguments voisins gardent leurs routes.
- [ ] Étendre la certification au-delà des familles rationnelles testées,
  notamment paramètres complexes, symboliques et irrationnels. La garde contre
  les faux pôles dus à l'arrondi reste requise.
- [ ] Étendre l'audit des fonctions multivaluées au-delà des points testés,
  en mesurant les chiffres corrects, les annulations et le côté des coupures.
- [ ] Réexaminer la borne de stockage des exposants binaires (65 536 bits)
  et le coût d'affichage près de cette borne. `exp(2**70000).evalf()` reste
  symbolique ; SymPy n'a pas terminé en 15 minutes dans la sonde documentée.
- [ ] Étudier les lenteurs d'`elliptic_pi` près des coupures et les dernières
  différences de chiffres des fonctions spéciales ; distinguer les cas où
  les deux moteurs ne terminent pas des écarts de capacité avérés.

Preuves : [branchcut](docs/algorithmic-parity.md#branchcut-wave-16-analytic-continuation-on-cuts-branch-points-and-poles),
[grands exposants](docs/algorithmic-parity.md#arbitrary-size-binary-exponents-bigexp2-wave-8-2026-10-23).

## 6. Simplification, matrices et expressions non commutatives

- [ ] Compléter les ordres symboliques de `Derivative` (issue 22210) :
  représentation native exacte et opérations gardées validées ; formules
  symboliques, transformations génériques et protocoles spéciaux restent partiels.
  La substitution par un Float ne doit pas devenir un ordre entier exact.
- [x] Calculer les grands ordres littéraux Int/BigInt certifiés : affine
  exp/sin/cos/sinh/cosh, puissances affines et dérivées au-delà/au sommet du
  degré polynomial, avec coefficients numériques fermés finis exacts.
  Les autres fonctions/domaines gardent le parcours complet ; aucune preuve
  déclinée ne devient zéro. [Preuves et suite complète](audits/sympy-1.14.0/native-literal-derivatives-2026-10-27/README.md).
- [x] Conserver les grandes puissances de recettes dans substitutions, pickle,
  expand et CSE ; Pow/func/int/doit ordinaires restent complets et num_pow
  inchangé. Exports SymPy canoniques ordinaires et callbacks rétablis, sans
  retirer les assertions qui ont détecté les deux régressions.
- [x] Corriger le 0/0 radical de `radsimp` (régression SymPy issue 5934,
  précédemment mal étiquetée 15119) : preuve algébrique exacte du zéro et
  `nan`, y compris sous wrappers et options. Cela ne clôt pas tous les
  calculs non finis ni les zéros indécidables.
- [ ] Résoudre les crashs `hyperexpand` plan/ushift et les lenteurs des
  opérateurs de décalage documentés ; consolider les autres écarts de
  `cse`, `fu`, `radsimp`, `powsimp`, `trigsimp` et `nsimplify`.
- [x] Faire suivre à `eigenvects()` le polynôme caractéristique simplifié des
  matrices contenant des identités symboliques, comme pour `eigenvals()`.
- [ ] Compléter les formes de dérivées de traces imbriquées et les règles
  de simplification des tableaux, contractions et fonctions élément par élément.
- [ ] Requalifier les dernières formes de QR, pseudo-inverse, exponentielle
  de matrice exacte et résolution paramétrique ; préserver précision,
  spécialisations et décisions de nullité dans les chemins concernés.
- [ ] Étendre les opérations non commutatives encore explicitement refusées
  quand une méthode correcte existe ; conserver les gardes des algorithmes
  qui supposent la commutativité et les corrections intentionnelles via LU.

Preuves : [simplification](docs/algorithmic-parity.md#simp-rest6-wave-15),
[matparam](docs/algorithmic-parity.md#matparam-wave-16),
[dérivées matricielles](docs/algorithmic-parity.md#mat-rest8-wave-14),
[non commutatif](docs/algorithmic-parity.md#pw-nc-wave-15-piecewise-and-noncommutative-leftovers).

## 7. Fonctionnalités spécialisées et compatibilité de sortie

- [ ] Représenter entièrement `Idx` dans les indices composés tels que `i+1`,
  avec leurs bornes. L'identité native des références à `Idx` directs et leur
  sérialisation sont corrigées ; la promotion historique en simple label dans
  un indice composé reste distincte de ce correctif.
- [ ] Stats : états `Str` de `DiscreteMarkovChain` et paramètres
  `MatrixSymbol` de dimension symbolique pour `MultivariateNormal`.
- [ ] Holonome : domaines flottants, puissances symboliques ou négatives,
  réécriture `erfc` et conditions initiales singulières restantes.
- [ ] Fonctions de matrices : construction directe de `sin(MatrixSymbol(...))`
  encore refusée, alors que la substitution dans `sin(x)` est prise en charge.
- [ ] `lambdify` : accès au source via `inspect.getsource` d'une
  `CompiledFunction` ; distinguer les différences de types scalaires NumPy.
- [ ] Impression et erreurs : alias des corps algébriques, noms de types,
  objets Python en plusieurs lignes, derniers chiffres, noms de `Dummy`,
  formes de `nsolve` et messages d'erreur. Reprendre les restes du rejeu
  `audit15` ; exclure l'ordre aléatoire des ensembles et les sorties de backends
  dépendant de l'environnement des bugs mathématiques.

Preuves : [audit9](audits/sympy-1.14.0/triage-wave9.md),
[audit15](docs/algorithmic-parity.md#audit15-wave-15),
[holo-rest](docs/algorithmic-parity.md#holonomic-multiple-derivatives-triangulated-extension-and-parametric-systems-holo-rest-2026-10-23).
Les groupes de Galois de degrés 5/6, `g(mat)` en Julia/Octave, `to_meijerg`,
`fps` symbolique et les anciens lots c4 ne sont plus des tâches ouvertes.

## 8. Capacités limitées par les ressources et performances

- [ ] Accélérer `count_roots(x**5000 - x - 1)` : plus de 60 s contre
  13,7 s pour SymPy dans la sonde documentée.
- [ ] Lever le blocage de `factor_list(x**300 + x + 1)` par estimation du
  travail avant recherche, là où SymPy termine en environ 9 s.
- [x] Réduire les intermédiaires des grands résultants :
  `resultant(expand((x + 1)**200), expand((x - 2)**150 + 1), x)` atteint
  encore la limite de stockage de la chaîne de sous-résultants. Fait en
  vague 20 (polycap3) : PRS dense exacte sur ZZ, 3,9 s contre 16,7 s pour
  SymPy avec `includePRS=True`.
- [ ] Étendre les bornes de `cyclotomic_poly(20001, x)` et
  `decompogen(x**4100 + x**2, x)` là où la référence termine.
- [ ] Évaluer exactement au-delà des plafonds actuels : `bernoulli(5002)`,
  `legendre(2001, Rational(1, 3))`, `fibonacci(1000001)`, `catalan(100001)`
  et `factorial(100001)`.
- [ ] Réduire les épuisements mémoire de l'intégration rationnelle et de la
  simplification de ses primitives ; conserver les cas documentés, dont
  `integrate(1/(x**5 + x + 1)**3, x)`, même quand SymPy ne termine pas non plus.
- [x] Reprendre l'inventaire des budgets des solveurs et des racines :
  `solve.rs` et `poly/roots.rs` n'ont été que documentés dans `reslimits`.
  Vérifier les chemins d'épuisement restants et tout nouveau plafond. Fait
  en vague 20 (polycap3) : aucun plafond ne refuse une entrée que SymPy
  calcule ; reste l'écart de vitesse des tours `exp` imbriquées.

Preuves : [inventaire et cas ouverts de reslimits](docs/algorithmic-parity.md#resource-limits-wave-17).
Les timings ci-dessus sont des mesures historiques à reproduire, sans délai
mural ajouté au moteur. Un budget épuisé doit lever une erreur explicite,
laisser le calcul non évalué ou conserver un repli complet ; il ne prouve
jamais l'absence de solution, l'irréductibilité ou une valeur nulle.

## Corrections intentionnelles et entretien des preuves

Conserver les résultats mathématiquement corrects lorsque SymPy 1.14 se
trompe : familles de solutions supplémentaires, exclusions de pôles,
intégrales divergentes retenues, inégalités périodiques, signes des
factorisations gaussiennes, matrices non commutatives, `RootSum.subs`,
séries de `Piecewise` et précision près des pôles. Les décisions utilisateur
et leurs justifications restent dans l'archive et le contrat détaillé.

- [ ] Dans le nouveau tri, distinguer résultat faux, solution incomplète,
  précision insuffisante, domaine absent, erreur différente, forme équivalente,
  correction intentionnelle et différence d'environnement.
- [ ] Réconcilier les passages historiques contradictoires du contrat détaillé
  avec les correctifs ultérieurs ; ce TODO est l'index actif des tâches.
- [ ] Mettre à jour ce fichier, le contrat et les preuves à chaque écart fermé,
  sans annoncer une parité au-delà des domaines effectivement validés.

Vérification initiale du nettoyage, avant ce lot : **284 tests ciblés réussis le 26 octobre 2026**
sur l'extension locale existante (sans reconstruction) : `limitseq_dominant`,
`holo_rest`, `nroots_precision`, `ccdom`, `matparam`, `solve_checks`,
`solve_decomposed_radicals` et `general_simplify_controls`. Ce contrôle ne
remplace ni le rejeu complet de l'audit ni la suite complète sur une build neuve.

Dix sondes comparatives supplémentaires ont reproduit les écarts de `sec`,
exponentielle imbriquée symbolique, ordre symbolique de dérivée, inférence de
`CC`, conversion `CC -> QQ_I`, `cancel` sur `RR`, singularité de `hyper`,
`sin(MatrixSymbol)` et limite harmonique. Elles ne constituent pas un audit
exhaustif. L'extension utilisée est `python/symox/_core.abi3.so`, SHA-256
`edb43bd43ec222c8ab6f08a074d5958157ac67b47d333a5550262a4d8d8ebeed`.

Premier lot du 26 octobre : [algorithmes, portée et écarts restants](docs/algorithmic-parity.md#follow-up-on-2026-10-26-trig-preimages-inexact-domains-and-numerical-poles).

Validation du premier lot : **178 nouvelles régressions**, toutes incluses dans la
suite complète Demeter ; aucun test collecté omis. Binaire natif SHA-256
`8323eddaa68e4c1f2c55cff6933c84fa3c4af204ed56eabf4e25a4d35bec1733`. Les
[preuves](audits/sympy-1.14.0/validation-2026-10-26.json) distinguent résultats
définitifs, skips, xfail, environnement et écarts de capacité encore ouverts.
