# TODO — Capacités mathématiques et parité SymPy

Mis à jour le **6 novembre 2026**, sur la branche `integrate-wave42-20261106`
(`main` à **`e4b1d011`**, vague 41, plus les trois lots de la vague 42).
Référence : **SymPy 1.14.0**. Dernier lot validé sur Demeter le 6 novembre
(voir « Reprise sur Demeter » ; le paragraphe ci-dessous décrit la suite du
28 octobre).

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

## Reprise sur Demeter (état au 6 novembre)

**Résultat faux prioritaire résolu** (vague 41, lot unrad) :
`solve(-2*y/sqrt(y**3) - C1 - x, y)` rendait `[]` (absence de solution
affirmée à tort) ; il rend maintenant `[4/(C1**2 + 2*C1*x + x**2)]` comme
SymPy, par la branche unrad native de `_solve`, atteinte avant les étapes
`_tsolve` comme dans la référence. Les cas EDO fact_13 et fact_16 sont
débloqués. La réserve de la vague 41 (radicaux de l'inconnue mêlés à
d'autres fonctions de celle-ci, `sqrt(x)*exp(x)`, exposants symboliques,
qui passaient par les replis Python appelant `_tsolve` avant unrad) est
levée par la vague 42 (lot tsolve) : repli unrad puis `_tsolve` natifs,
dans l'ordre de SymPy ; le dernier WRONG_EMPTY du corpus tsolve est
corrigé. Restent en Python la chaîne établie après un NotImplementedError
natif et deux crochets (`nsimplify`, fin d'`equals`) ; voir « Lacunes
ouvertes principales ».

La branche `integrate-wave42-20261106` (main `e4b1d011` + vague 42 :
tsolve, ode3, prde3, fusionnés dans cet ordre) est validée, à avancer dans
main. Validation (`XLA_PYTHON_CLIENT_PREALLOCATE=false`) sur `dfc19715`
(fusion des trois lots, code identique au commit TODO) : suite complète
**60 482 réussis, 0 échec** ("SUITE CLEAN", 6 xfail, 45 sautés,
60 533 collectés et exécutés, 693 fichiers ; preuves dans
`~/symox-ci/evidence/integ42-20261106`), soit exactement l'union des trois
lots (+27, +73 et +28 tests collectés et 3 fichiers sur les 60 405 et 690
de main ; 2 xfail stricts levés par prde3, 8 → 6) ; cargo test
871 + 1 + 0 + 5, clippy correctness propre ; corpus EDO rejoué sur la
branche intégrée : **315/289/270**, verdicts identiques ligne à ligne à
ceux de main (`integ41-20261106`) ; corpus tsolve rejoué : 106 SAME,
6 EQUIV, 2 SUPERSET, 1 DECLINED, 0 WRONG_EMPTY, comparaison identique
ligne à ligne à celle du lot ; corpus unrad rejoué : 87 SAME, 2 EQUIV,
1 EXTRA, identique au lot (preuves dans `.../integ42-20261106/solve`).
Interactions vérifiées : `nth_algebraic` (ode3) appelle le `solve` public
modifié par tsolve (lignes algeb du corpus EDO et 73 tests ode3
inchangés) ; l'ordre global des générateurs de `construct_domain` (prde3)
ne change aucun verdict des corpus EDO, tsolve et unrad ni aucun test.
Conflits de fusion limités au texte de `TODO.md` (fusion de prde3 :
entrées vague 42 d'ode3 et de prde3 gardées toutes deux ; lacunes solve de
tsolve et prde de prde3, les plus récentes) ; `docs/algorithmic-parity.md`
fusionné sans conflit (sections ajoutées bout à bout). Aucun correctif
d'intégration nécessaire. Le corpus d'intégration n'a pas été rejoué
depuis la vague 37.

La branche `integrate-wave41-20261106` (vague 41 : unrad, accumb,
arrcomp) est dans main (`e4b1d011`) : suite 60 352 réussis (60 405
collectés, 690 fichiers, 8 xfail), cargo test 867 + 1 + 0 + 5, corpus EDO
315/289/270.

Vague 43 (en cours ; calcul en Rust, Python limité à l'admission,
l'aiguillage et les enveloppes) :
- matcher (`feat-matcher-20261106`) : filtreur de référence natif
  (`crates/symox-core/src/matches.rs`, liaison `_core._reference_match`),
  portage du protocole `match` de SymPy 1.14 (`Basic.matches`,
  `Wild.matches`, `AssocOp._matches_commutative`, `_combine_inverse` de
  `Mul`/`Add`, `Pow.matches`, `Mul.as_base_exp`), mêmes parcours et donc
  même choix entre plusieurs correspondances valides ; il décline
  (repli) plutôt que de répondre autrement.
  - `Expr.match` public (et `matches` sans `repl_dict`) : natif, repli
    sur `python/symox/_matching.py` quand le natif décline.
  - filtreurs EDO natifs `second.rs`, `liouville.rs`, `reducible.rs`
    (auparavant `pattern::pattern_match`) et `bessel.rs` (auparavant le
    rappel Python `Expr.match`, retiré de `Hooks`, `ode_classify.rs` et
    `_classify_native.py`) : filtreur de référence natif, déclin = non
    modélisé.
  - corpus différentiel (`audits/performance/2026-11-06-matcher`, 2311
    appels : 1563 EDO, 248 `test_match.py`, 500 aléatoires) : natif 2226
    identiques à SymPy, 0 différent, 85 déclinés ; `Expr.match` public
    2309 → 2311 identiques (les deux `d**2` sur `y**2*cos(4)**2`
    corrigés par `Mul.as_base_exp`) ; `pattern_match` 2059/252. Aucune
    des 99 différences EDO de `pattern_match` ne vient d'un motif de
    `second.rs`, `liouville.rs`, `reducible.rs` ou Airy (914 appels) :
    pas de réponse EDO fausse observée avant le lot. Temps cumulés :
    SymPy 25,9 s, Python 14,8 s, natif 1,9 s ; les filtreurs EDO natifs
    passent d'environ 0,06 ms à 1 ms par appel (protocole de référence).
  - corpus EDO 315/289/270, verdicts identiques ligne à ligne à
    `integ42-20261106` ; suite **60 490 réussis, 0 échec** ("SUITE
    CLEAN", 60 541 collectés, 694 fichiers, 6 xfail ; course 3 de
    `~/symox-ci/evidence/matcher-20261106/suite-run3` ; les courses 1 et
    2, sous charge 40, avaient perdu le shard 23 sur le cgroup de 6 Go et
    deux tests sur délai, sans échec du lot) ; cargo test
    872 + 1 + 0 + 5 (+1 : `matches::tests::reference_choices`), clippy
    correctness propre.
- limits2 (`feat-limits2-20261106`) : écarts limites/séries à SymPy
  réduits en Rust (`audits/performance/2026-11-06-limits2/`). Asserts
  des fichiers de test SymPy identiques : 1139 → 1211 sur 1293 (aseries 5
  → 5, demidovich 53 → 54, formal 59 → 59, gruntz 103 → 103, limits 533 →
  538, lseries 5 → 5, nseries 165 → 172, order 172 → 205, series 44 →
  70), 72 corrigés, aucun cassé ; déclins Gruntz 87 → 84. Natif :
  - formes de `Order` (termes de l'Add développé, facteurs à `x` lié
    écartés, produit de fonctions gardé sur un pôle), facteurs non
    commutatifs absorbés (`Mul.flatten`), `conjugate`, `subs` par une
    valeur portant son propre ordre, renormalisation des ordres par
    `expand` ;
  - Subs/Derivative/Integral gardés : `_eval_nseries`, boucle
    `Expr.series`, `Derivative._eval_as_leading_term`
    (`series/held.rs`) ;
  - Taylor générique des fonctions indéfinies (Dummy `xi`), acoth (forme
    log en ±1 et sur la coupure), terme dominant à exposants réels
    symboliques ;
  - limites par un Dummy positif (`limitinf`), valeur prouvée nulle par
    les hypothèses, `heuristics` de référence avec bornes
    (`limits/heuristics.rs`), `Limit` à variable non symbole et ses
    `free_symbols`.
  Python (admission et API seulement) : `Order.as_expr_variables`,
  `Expr._eval_as_leading_term`, `lseries` sans `x`, admission
  `has_free` de `as_leading_term`/`leadterm`. Le test
  `test_bound_calculus_compat` attendait NotImplementedError pour
  `Limit(f(x) + y, x, 0).subs(x, 2)` ; SymPy rend `Limit(y + f(2), 2, 0)`,
  comme symox désormais.
  Validation (`XLA_PYTHON_CLIENT_PREALLOCATE=false`) : suite complète
  **60 503 réussis, 0 échec** ("SUITE CLEAN", 6 xfail, 45 sautés,
  60 554 collectés et exécutés, 694 fichiers ; +21 de
  `tests/test_native_limits2.py` ; preuves dans
  `~/symox-ci/evidence/limits2-20261106`) ; cargo test 871 + 1 + 0 + 5 ;
  clippy correctness propre. Écarts à valider et lacunes : voir la ligne
  Limites des lacunes ouvertes et la section limits2 de
  docs/algorithmic-parity.md.

Vague 42 (calcul en Rust, Python limité à l'admission,
l'aiguillage et les enveloppes) :
- tsolve (`feat-tsolve-20261106`) : repli `_tsolve` de `_solve` natif
  (`crates/symox-core/src/solve/native_tsolve.rs`), dans l'ordre de SymPy :
  - repli unrad de `_solve` avant `_tsolve` pour les radicaux mêlés à
    d'autres fonctions de l'inconnue (`radical::solve_mixed`,
    `sqrt(x)*exp(x) - 1`, `exp(sqrt(x)) - 3`) ;
  - étapes `_tsolve` (`solve/tsolve.rs`) avec un hôte natif : `_invert`,
    factor/logcombine, `Pow` (candidats rationnels), `sin`/`cos`/`LambertW`,
    `rebuild(rewrite(exp))`, bloc bivarié (`_solve_lambert`,
    `bivariate_type`), reprise posify (drapeaux sans `bivariate`) ;
  - `_vsolve` imbriqués natifs, sans filtrage par les hypothèses de
    l'inconnue comme la référence (`u` positif de `bivariate_type`,
    posify : `x**x - x` garde `-1`) ;
  - étape multi-générateurs de `_solve` (`_as_base_q`, substitution de la
    fonction la moins profonde) ; logcombine limité aux termes en `x`
    (`(x**3)**(x/2) + pi/2` rendait `[]`) ; partie primitive d'une
    quadratique à coefficients symboliques avant la forme close, comme
    `preprocess_roots` (`solve.rs::quadratic_primitive`) ;
  - indécidable natif : déclin et relance, jamais `[]`.
  Corpus de 115 équations (`audits/performance/2026-11-06-tsolve/`) : avant
  99 SAME, 6 EQUIV, 3 SUBSET, 1 SUPERSET, 1 WRONG_EMPTY, 5 DECLINED ; après
  106 SAME, 6 EQUIV, 2 SUPERSET, 1 DECLINED ; cas passant par la chaîne
  Python 66 → 10. Corpus unrad inchangé (87 SAME, 2 EQUIV, 1 EXTRA).
  Validation (`XLA_PYTHON_CLIENT_PREALLOCATE=false`) : suite complète
  **60 379 réussis, 0 échec** ("SUITE CLEAN", 8 xfail, 45 sautés, 60 432
  collectés, 691 fichiers ; +27 de `tests/test_native_tsolve_fallbacks.py` ;
  preuves dans `~/symox-ci/evidence/tsolve-20261106`) ; cargo test 867 + 1
  + 0 + 5 ; clippy correctness propre ; corpus EDO 315/289/270 (inchangé,
  aucune ligne modifiée). Écarts à valider : branches complexes de Lambert
  en SUPERSET (`x**3 - 3**x`, `3**cos(x) - cos(x)**3`), `_vsolve`
  imbriqués sans filtrage par hypothèses, forme de SymPy pour
  `almost_lin_01` (`tests/test_ode_solve_parity.py`) ; voir la section
  wave 42 de docs/algorithmic-parity.md.
- ode3 (`feat-ode3-20261106`) : solveurs EDO `nth_algebraic`,
  `factorable` et `2nd_linear_bessel` natifs, classification et
  construction de la solution (`crates/symox-core/src/ode/algebraic.rs`,
  `ode/bessel.rs`, liaisons `ode_nth_algebraic`, `ode_factorable_match`,
  `ode_factorable_solutions`, `ode_bessel_match`, `ode_bessel_solution`
  dans `crates/symox-py/src/ode_classify.rs`). Les appels publics
  (`solve`, `dsolve`, `simplify`, `factor`, `cancel`, `expand`,
  `Expr.match`) sont des rappels de l'hôte ; un résultat non modélisé
  décline vers le corps de référence. Exemples SymPy de `test_single.py`
  pour ces trois indices : 23/23, 20/20, 12/12 avant et après ; sur les
  315 exemples de `test_single.py`, 791 décisions natives, aucun
  déclin. `test_single.py`/`test_ode.py` copiés : 53 réussis, 20 échecs,
  identiques avant et après. Corpus EDO 315/289/270 avant et après. `lie_group` reste en
  Python (non porté dans ce lot). Validation du lot : suite complète
  **60 425 réussis, 0 échec** ("SUITE CLEAN", 8 xfail, 45 sautés,
  60 478 collectés, 691 fichiers ; +73 tests et 1 fichier :
  `tests/test_native_ode3.py` ; un premier passage avait échoué sur le
  délai de 5 s d'un sous-processus de `test_dirac_affine_integration.py`
  sous charge, réussi isolément puis au second passage) ; cargo test
  871 + 1 + 0 + 5 ; clippy correctness propre.
- prde3 (`feat-prde3-20261106`) : EDR de Risch paramétrique en Rust sur les
  niveaux tan et non linéaires (`crates/symox-py/src/prde_param.rs`) :
  `prde_no_cancel_b_small` hors niveau de base, `NotImplementedError` de
  `prde_no_cancel_b_equal` et des cas non linéaires, bras tan de
  `prde_special_denom` (`real_imag`, `recognize_log_derivative`),
  `bound_degree` tan/non linéaire, `v == 0` de `param_rischDE`, appels
  publics `prde_no_cancel_b_*` au-dessus de la base ; `PolynomialError`
  kkinv de SymPy reproduite (cas 44, xfail strict levé) ; ordre des
  générateurs `_sort_gens` (`QQ(x,t0)`) dans le noyau et, par aiguillage
  vers le natif, dans `construct_domain` ; `Poly + None` lève
  `SympifyError` (liaison Rust, xfail `risch_prim` levé). Python :
  aiguillage seulement. Mesures (`audits/performance/2026-11-06-prde3/`) :
  582 appels aléatoires tan/non linéaires, réponses natives 64 → 376,
  identiques à SymPy 555 → 573, natif 112,8 s → 12,3 s (SymPy 138,9 s) ;
  584 appels aléatoires du lot prde2, natifs 461 → 506, identiques
  536 → 554 ; 142 + 485 appels enregistrés identiques avant et après ;
  aucun écart propre au natif. Suite 60 382 réussis / 60 433 collectés
  (691 fichiers, 6 xfail, SUITE CLEAN) ; `cargo test --release`
  867 + 1 + 0 + 5, clippy correctness propre.

Vague 41 (dans main ; calcul en Rust, Python limité à l'admission,
l'aiguillage et les enveloppes) :
- unrad (`feat-unrad-20261106`) : branches radicales de `_solve` natives
  (`crates/symox-core/src/solve/radical.rs`), dans l'ordre de SymPy, avant
  `_tsolve` :
  - un seul générateur : racines du polynôme en `g` puis inversion de `g`
    (`y**(3/2) = x` a ses trois solutions) ;
  - plusieurs générateurs : `unrad` natif, avec un hôte natif, puis la
    résolution imbriquée, changement de variable compris ;
  - vérifications existantes et `check=False` ; `_unrad=False` accepté ;
    valeurs non simplifiées comme la référence (lues par l'enveloppe
    Python, aiguillage seul) ;
  - le `solveset` natif garde son `radical_step` ;
  - un cas indécis est décliné, jamais `[]`.
  Corpus de 90 équations radicales (`audits/performance/2026-11-06-unrad/`) :
  avant 63 SAME, 12 EQUIV, 7 SUBSET, 4 WRONG_EMPTY, 2 DECLINED, 1 SUPERSET,
  1 EXTRA ; après 87 SAME, 2 EQUIV, 1 EXTRA. Corpus EDO
  315/287/268 → 315/289/270 (fact_13, fact_16 identiques à SymPy).
  Validation (`XLA_PYTHON_CLIENT_PREALLOCATE=false`) : suite complète
  **60 329 réussis, 0 échec** ("SUITE CLEAN", 8 xfail, 45 sautés,
  60 382 collectés et exécutés, 688 fichiers ;
  `~/symox-ci/evidence/unrad-20261106-r2`), soit +11 tests nouveaux et
  3 xfail stricts levés ; cargo test 865 + 1 + 0 + 5, clippy correctness
  propre. Deux passages antérieurs, sous charge 30 à 56 (agents
  concurrents), ont échoué uniquement sur des délais : l'OOM du lot
  `tests/extra` à 6 Go, que le binaire de main atteint aussi, puis deux
  délais de sous-processus (fuzz 20 s, Dirac 5 s), aussi lents sur main.
- accumb (`feat-accumb-20261106`) : arithmétique `AccumBounds` native
  (`limits/accum.rs` : constructeur, opérateurs et opérateurs réfléchis,
  puissances, `_eval_is_ge`, intersection, union, appartenance, repli des
  coefficients de `Add`/`Mul`), valeurs des fonctions sur des bornes
  (`sin`, `cos`, `tan`, `cot`, `exp`, `log`, `frac`, hyperboliques en
  `±I*oo`) et marcheur des limites oscillantes (`oscillating_limit`) en
  Rust ; 873 lignes de Python retirées. `log(AccumBounds(a, b))` avec `a`
  négatif ou de signe inconnu vaut `nan` comme la référence (auparavant
  laissé tel quel, d'où `limit(log(sin(x)), x, oo)` = `log(AccumBounds(-1,
  1))` là où SymPy garde la limite). Asserts limites 1139/1293 inchangés,
  les 38 limites que le marcheur Python calculait passent par le natif,
  rejeu de 288 appels identique, déclins Gruntz 88 → 87 (la valeur
  `nan` de `log(cos(oo))` arrête la marche avant le déclin, limite gardée
  comme SymPy), limites oscillantes 1,6 à 3,7 fois plus rapides. Suite
  complète du lot : 60 331 réussis, "SUITE CLEAN" (60 387 collectés, 688
  fichiers, `~/symox-ci/evidence/accumb-20261106`) ; cargo test 867 + 1 +
  0 + 5 ; clippy correctness propre. Preuves :
  `audits/performance/2026-11-06-accumb/`.
- arrcomp (`feat-arrcomp-20261106`) : modèles natifs kronecker
  (`combine_kronecker`, identités, produit explicite, `expand`, entrée
  symbolique, conjugué), erreurs matexpr de SymPy (`as_explicit` à forme
  symbolique, NonSquareMatrixError, `2/B`, clé plate sans repli), `Eq`/`Ne`
  à dimensions prouvées différentes, attributs `MatrixElement`/`Trace`,
  `Matrix[i, j]` symbolique en `MatrixElement`, et `ArrayComprehension`
  calculé en Rust (`crates/symox-py/src/array_comprehension.rs`, classe
  Python réduite à l'enveloppe). Tests SymPy de tableaux : 165 → 182
  sur 201, sans régression. Validation (`XLA_PYTHON_CLIENT_PREALLOCATE=false`) :
  suite complète **60 322 réussis, 0 échec** ("SUITE CLEAN", 11 xfail,
  45 sautés, 60 378 collectés, 688 fichiers ; +7 de
  `tests/test_native_arrcomp.py` ; preuves dans
  `~/symox-ci/evidence/arrcomp-20261106`) ; cargo test 865 + 1 + 0 + 5,
  clippy correctness propre. Écarts délibérés à valider : voir la liste
  de la vague 41 ci-dessous et la section arrcomp de
  docs/algorithmic-parity.md.

Vague 40 (dans main ; calcul en Rust, Python limité à l'admission,
l'aiguillage et les enveloppes) :
- ode (`feat-ode-20261105`) : 2F1 natif (`ode/hypergeometric.rs` :
  `equivalence_hypergeometric`, `match_2nd_2F1_hypergeometric`,
  `get_sol_2F1_hypergeometric`), identique à la référence sur les 34 EDO
  linéaires d'ordre 2 du corpus ; `Integral(F,(u,b)).doit()` garde
  `Integral(g,(u,b))` (`integrate.rs`, `binders.rs`) ; recast de `f(x)` en
  Dummy au-delà des équations rationnelles (`solve/recast.rs`). Corpus
  315/285/267 → 315/287/268, sans régression.
- limits (`feat-limits-20261105`) : voir le paragraphe limits ci-dessous.
- prde2 (`feat-prde2-20261105`) : voir le paragraphe prde2 ci-dessous.

Vague 39 (dans main ; calcul en Rust, Python limité à l'admission,
l'aiguillage et les enveloppes) :
- ratintex (`feat-ratintex-20261030k`) : logpart sur QQ<a> (un corps, deux
  corps unifiés en QQ<b, a>), log_to_atan sur EX, ordre RootSum, roots
  cas 20, EX numer/denom de référence, ordre d'association des radicaux
  dans expand, ordre d'impression des radicaux ; `integrate(1/(x**3+x+1))`
  et `1/(x**3+2*x-1)` identiques à SymPy.
- arrexpr (`feat-arrexpr-20261105`), suite d'arr2mat :
  - `Mul`/`MatMul`/`Trace` à coefficient tableau natifs, ce qui fait passer
    le repli vendu de 86/527 à 6/527 cas ; ces 6 cas sont arr15, arr34,
    arr145, arr249, arr308 et arr336 (Hadamard, `DiagMatrix` ou transposée
    d'un opérande tableau : lacune) ;
  - corpus à 520/527 ;
  - `_chain_array_to_matrix` natif (repli Python gardé) ;
  - `_get_contraction_tuples` natif ;
  - `MatPow` imbriqués combinés, et `A/B` ;
  - ordre de `MatrixElement` et `ArrayElement` (arrindex passe de 162/168 à
    167/168 ; reste `u_int_mul`, ordre des arguments d'`ArrayDiagonal`) ;
  - `matrix_derivative` passe d'abord par le chemin tableau natif comme
    SymPy 1.14 (0,5 à 1,7 ms contre quelques µs avant : le rappel Python
    `matrix_sort_key` domine, suivi de performance) ;
  - tests SymPy des tableaux 165/201 (159 avant).
- prde (`feat-prde-20261105`) : `is_log_deriv_k_t_radical_in_field` et
  `parametric_log_deriv_heu` (donc `parametric_log_deriv`) sur les tours
  exp/primitive (`integrate/prde_tower.rs`) : 533 appels aléatoires,
  529 natifs, 4 déclins, aucun écart propre au natif ; appels des tests
  SymPy et du corpus identiques à main. Trace des Dummies du noyau de base,
  `Poly.free_symbols` des domaines composites et boucles d'admission de
  splitfactor (`risch_sf_admission.rs`) en Rust.

Détail des lots limits et prde2 de la vague 40 :
- limits (`feat-limits-20261105`) : API `sympy.series.gruntz` native
  (`limits/gruntz/api.rs`), `nsimplify` des flottants de `Limit.doit`
  (`limits/floats.rs`), étape terme dominant à exposant symbolique et
  `b**p` en `±oo` (`limits/symbolic_power.rs`), NotImplementedError de
  référence (signe, `zoo`), Bessel en `±I*oo`, coefficient `I*a`,
  `(-1)**k` (8634), hyperboliques (26991), `Limit.pow_heuristics` pour
  une base constante avec la clause `exp(±oo*t)` d'`exp.eval`
  (`constant_base_power`, `(I/2)**x` → 0, `(-2)**x` → `zoo`), `powsimp`
  avant `pow_heuristics` et une seconde passe native
  (`exp(b*x)*exp(-x)` lève l'erreur de signe de SymPy), signe d'un
  coefficient dominant développé lu après `factor_terms` : sonde
  symbolique 48/54 identiques à SymPy (15 avant) ; 1139/1293 asserts (1118 avant,
  21 corrigés, 0 cassé), test_gruntz 37/38, test_limits 158/173 ; déclins
  Gruntz 90 → 88, dont 10 lèvent maintenant l'erreur de SymPy.

- prde2 (`feat-prde2-20261105`, sur cab152e3) : `param_rischDE`,
  `param_poly_rischDE`, `prde_cancel_liouvillian` et `limited_integrate`
  au-dessus du niveau de base, sur les tours exp/primitive, en Rust
  (`symox-py/src/prde_param.rs`) : dénominateurs normal et spécial,
  contraintes linéaires, `constant_system`, `prde_spde`, `bound_degree`
  (y compris la branche primitive `db == da` qui appelle
  `limited_integrate`), b_large, b_small au niveau de base, annulation
  liouvillienne, `limited_integrate_reduce`. 584 appels aléatoires :
  0 écart propre au natif (536 identiques, 31 écarts déjà présents dans le
  corps Python, 17 dépassements SymPy à 20 s) ; appels des tests SymPy
  (142) et du corpus (485) identiques à SymPy, corpus 16,9 s contre
  37,1 s SymPy (42,7 s avant). Python : admission, aiguillage, Dummies,
  exceptions. Test cargo `cancel_include_matches_reference` et test
  splitfactor hermite_reduce alignés sur SymPy (QQ conservé depuis
  88545232) ; le test d'atteinte des noyaux de première vague
  (`test_native_prderust.py`) renvoie aussi les noyaux de tour au corps
  Python. Validation : suite **60 264 réussis, 0 échec** ("SUITE CLEAN",
  7 xfail, 45 sautés, 60 316 collectés, 685 fichiers ; +54 =
  `tests/test_native_prde_param.py`, preuves dans
  `~/symox-ci/evidence/prde2-20261106`) ; cargo test 859 + 1 + 0 + 5 ;
  clippy correctness propre.

Vague 38 (dans main) :
- arr2mat : `convert_array_to_matrix`, ses passes, `_EditArrayContraction`
  et `_ArgE` natifs (`arrayexpr_a2m.rs`). `Trace(MatMul)` n'est plus
  canonicalisé (aligné sur SymPy, validé le 30 octobre).
- arrindex : `array_derive`, `convert_indexed_to_array` et le plan
  d'indices de `convert_array_to_indexed` natifs (`arrayexpr_index.rs`).
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
  (`integrate/prde_log.rs`, 938 appels, 0 écart).

Écarts volontaires de la vague 42 en attente de validation utilisateur :
- prde3 : ordre global des générateurs de `construct_domain`
  (`polys/constructor.py._infer_domain`) : à plusieurs symboles, l'ordre
  vient maintenant des générateurs natifs de `_construct_composite`
  (`_sort_gens` de SymPy : `x`, `y`, `z` d'abord, `ZZ[x,a]`, `ZZ(x,t0)`,
  `ZZ[y,a,b]`) au lieu de `default_sort_key`, conservé seulement quand ces
  générateurs ne sont pas exactement les symboles libres ; cela change
  tout `construct_domain` à plusieurs symboles, pas seulement prde ;
  élargissement de `SympifyError` : toute opérande d'une opération
  arithmétique de `Poly` qui n'est pas une expression (et plus seulement
  `None`) lève `SympifyError(o)` comme le `_sympify` strict de
  `_polifyit`, là où l'erreur était un `TypeError` ordinaire.
- ode3 : le filtrage du `2nd_linear_bessel` natif passait par le
  `Expr.match` public (filtreur de référence porté en Python, rappelé
  depuis Rust) ; depuis la vague 43 (lot matcher) il passe par le filtreur
  de référence natif (`crates/symox-core/src/matches.rs`), sans rappel.
- tsolve : branches complexes de Lambert en plus (SUPERSET : `x**3 - 3**x`
  et `3**cos(x) - cos(x)**3` rendent aussi
  `-3*LambertW((1 ± sqrt(3)*I)*log(3)/6)/log(3)`, solutions valides ; déjà
  le cas avant le lot pour `x**3 - 3**x`) ; filtre par hypothèses abandonné
  dans les `_vsolve` imbriqués des étapes `_tsolve` natives, comme la
  référence (`checksol` seul : `u` positif de `bivariate_type`, posify où
  `x**x - x` garde `-1`), alors que la chaîne PyHost conservée filtre
  encore (`solve` public) ; forme de SymPy reprise pour `almost_lin_01`
  dans `tests/test_ode_solve_parity.py`.

Écarts volontaires de la vague 41 en attente de validation utilisateur :
- unrad : là où la référence propage une erreur (`unrad` qui lève
  PolynomialError ou IndexError, NotImplementedError du changement de
  variable), la branche native décline et les étapes suivantes répondent ;
  un repli établi peut alors donner une valeur (`sqrt(x) + root(x, 3) +
  root(x + 1, 5) - 2` : SymPy lève NotImplementedError, symox rend la racine
  `CRootOf` vérifiée du repli `_unrad_resultant`, comme avant le lot) ; le
  marquage « non simplifié » s'applique par égalité de valeur aux solutions
  de la même équation (une solution aussi produite sous forme simplifiée
  par une autre branche reste sous la forme native).
- accumb (écarts antérieurs au lot, conservés) : `limit(x*sin(x), x, oo)`,
  `exp(x)*sin(x)`, `x + x*sin(x)` en `oo` et `sin(1/x)/x` en 0 valent
  `AccumBounds(-oo, oo)`, `x*log(2 + sin(x))` vaut `AccumBounds(0, oo)`, là
  où SymPy rend `oo*sign(AccumBounds(...))` ; `limit(sin(x)*log(x), x, oo)`
  vaut `AccumBounds(-oo, oo)` et `Max(sin(x), 0)`, `x**sin(x)` en `oo` sont
  gardés, là où SymPy lève un `TypeError` interne (`_eval_nseries()` sans
  `cdir`).
- arrcomp : `combine_kronecker` sans la règle de somme
  (`KP(A, B) + KP(C, D) -> KP(A + C, B + D)`, fausse en général) et ne
  distribuant que les puissances entières (`KP(X, Y)**x` gardé ; SymPy rend
  `KP(X**x, Y**x)`, faux sur la branche principale), d'où l'échec de
  `test_KroneckerProduct_combine_mul` et `_combine_pow` ;
  `Eq(MatrixSymbol('P', n, 1), MatrixSymbol('Q', m, 1))` non évalué (SymPy
  rend False alors que `n = m` est possible) ; `ArrayComprehension.subs`
  garde la substitution sans capture (indice lié non substitué, indice
  libre du remplaçant renommé en Dummy entier ; SymPy substitue et
  capture) ; développement plafonné à 100 000 entrées (ValueError au-delà,
  plafond des tableaux denses existant) ; `SparseMatrix[i, j]` à indice
  symbolique a pour parent une matrice dense immuable (SymPy : la matrice
  creuse ; entrées et affichage identiques).

Écarts volontaires de la vague 40 en attente de validation utilisateur :
- ode : le recast rend pour une cible `f(x)` les valeurs d'un symbole
  (`sqrt(f)+f-x` résolu, il levait avant ; forme symox, pas celle de
  SymPy) ; le 2F1 natif décline plus large (Floats, racines non réelles ou
  non décidées) au profit du corps de référence.
- prde2 : trois tests réalignés sur SymPy (QQ unifié par `cancel`) ;
  raccourcis `from_polys`/`getitem` sans `unify` quand générateurs et
  domaine sont égaux ; déclin (corps Python d'origine) au lieu de
  reproduire les erreurs d'origine ; Dummy `z` privé du noyau, absent des
  résultats comparés.
- limits : `gruntz()` se replie sur la limite native (étapes heuristiques)
  quand l'algorithme de Gruntz seul échoue, là où SymPy lève ; un flottant
  sans petit rationnel à plus de 12 chiffres significatifs (ou hors
  `[1e-12, 1e15)`) est gardé tel quel dans `Limit.doit` (SymPy peut y
  trouver une forme close via `identify`) ; `I**x` en `oo` reste non
  évalué (SymPy rend `nan`).

Écarts volontaires de la vague 39 en attente de validation utilisateur :
- ratintex : `nested_binomial` de roots (binômes à radicaux imbriqués,
  hypothèses is_extended_real/is_positive) ;
- arrexpr : `matrix_derivative` décline quand une feuille scalaire dépend
  de `x` (`exp(Trace(X))`, `X[0, 1]`, `DotProduct(a, a)`), où SymPy rend un
  zéro faux, et la chaîne précédente répond ; NotImplementedError levée dès
  la conversion, et non à l'affichage, pour un `MatMul` à coefficient
  tableau ;
- prde : lire `_core.__dict__` une seule fois par appel de splitfactor
  (au lieu d'environ 395 fois) changerait le contrat d'authentification ;
  non fait, à valider.

Écarts volontaires antérieurs en attente de validation : Markov spec7 1/6,
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
  d'erreur du casus irreducibilis ; lie_group (Python, non porté par
  ode3 : lie_group_01/02/03 limites de ressources, 08/09/12 formes
  d'intégration) ; nth_algebraic, factorable et Bessel natifs depuis la
  vague 42, mais la classe `diffx` de nth_algebraic reste créée en Python
  (son inverse tourne dans le `solve` public, rappel de l'hôte) et le
  filtrage du Bessel est natif depuis la vague 43 (voir la lacune
  suivante) ; algeb_04
  (ordre des deux solutions de `solve` sur une équation à plusieurs
  générateurs : chemin de changement de variable de `_solve` dans SymPy
  contre l'étape factor de `_tsolve` natif ; domaine solve). 2F1 natif depuis la vague 40 (corpus 315/285/267 →
  315/287/268) ; fact_13/fact_16 identiques depuis la vague 41 (branche
  unrad native ; corpus 315/289/270 sur la branche intégrée) ;
  conserved_03/04 (dsolve par défaut) dépassent la mémoire dans le
  `checksol` de `_is_special_case_of` (SymPy y dépasse aussi son délai) ;
  `integrate(x/(x**2 - 1))` garde deux
  logarithmes ([détail](docs/algorithmic-parity.md#ode-2026-11-05-native-2f1-step-upper-limit-integraldoit-fx-solve-recast)).
- Filtrage (vague 43) : `Expr.match` public et filtreurs EDO natifs
  (`second.rs`, `liouville.rs`, `reducible.rs`, `bessel.rs`) passent par
  le filtreur de référence natif `matches::reference_match` (protocole
  SymPy 1.14 porté en Rust), qui décline plutôt que de répondre
  autrement ; corpus différentiel de 2311 appels (EDO, `test_match.py`,
  aléatoires) : 2226 identiques à SymPy, 0 différent, 85 déclinés. Restent
  en Python (repli du public sur `python/symox/_matching.py`) : produits
  non commutatifs, `zoo`/`oo`/`nan`, `Sum`/`Integral` liés, relationnels et
  logique, `Piecewise`, matrices, Float face à un côté sans symbole libre
  dans l'inverse multiplicatif, `collect` à plus de 3 symboles libres de
  coefficients, `WildFunction` (classe Python), sous-classes hôtes
  définissant `matches`, appels avec `repl_dict`. Le `has` natif public
  ignore toujours le sous-ensemble commutatif de `Mul`
  (`(2*x*y*z).has(x*y)` vaut False contre True dans SymPy) ; le filtreur
  natif en a sa propre version. `pattern::pattern_match` (utilisé par
  `simp.rs`, meijerint et `collect`) n'est pas le protocole de référence
  (252 différences sur le corpus) ; ses appelants n'ont pas été audités
  ([détail](docs/algorithmic-parity.md#matcher-wave-43-native-reference-matcher-2026-11-06)).
- solve (radicaux, `_tsolve`) : repli `_tsolve` natif depuis la vague 42 ;
  restent en Python `nsimplify` d'un nombre (`mpmath.identify`, crochet)
  et la fin de `Expr.equals` indécise, `force=False` et `bivariate=False`
  publics, l'inverse d'une classe de fonction inconnue de la table native,
  la chaîne Python après un NotImplementedError natif (conservée :
  `_solve_piecewise`, masquage, inversion, `_unrad_resultant`), dont le
  `_tsolve` PyHost qui filtre encore ses solutions imbriquées par
  hypothèses ; équations imbriquées que l'enveloppe Python réécrit d'abord
  (Floats, hyperboliques, `Abs`, `Piecewise`, dérivées, intégrales, sommes,
  fonctions indéfinies), déclinées par l'hôte natif ; post-traitement
  Python (`_high_order_roots`, `_polynomial_roots`) qui appelle
  `solve(gen - t, x)` et repasse par la chaîne Python (sans effet sur le
  résultat) ; algeb_04 (ordre des solutions à plusieurs générateurs) ; branches
  complexes LambertW en plus (`x**3 - 3**x`, SUPERSET antérieur) ; formes
  `cov_not_sym`/`ratio_den` différentes (valeurs égales) ; `_take` lu par
  le seul `unrad` public, comme la référence
  ([détail](docs/algorithmic-parity.md#solve-2026-11-06-wave-42-native-_tsolve-fallback-and-unrad-before-_tsolve-for-mixed-equations)).
- prde : restent, en corps Python stock, les domaines de base à générateurs
  non symboles (`log(5)`, refusés à l'admission ; 4 appels enregistrés des
  tests, aucun dans les corpus), `prde_cancel_liouvillian` appelé
  directement sur un niveau tan/non linéaire (`UnboundLocalError` de
  SymPy, message dépendant de la version de Python), `real_roots` hors ZZ
  et les autres déclins listés dans le contrat (section prde3), dont le
  bras tan qui appelle `parametric_log_deriv_heu` avec le monôme courant
  dans le domaine des coefficients ; troncature `v == 0` à `W` non vide
  portée mais jamais rencontrée (la `PolyMatrix` stock lève `ValueError`
  là où SymPy tronque) ; écarts
  hérités restants : messages `TypeError` nommant `symox._core.Poly`,
  entrées brutes de `M` non simplifiées (formes publiques identiques), base
  `h` d'échelle différente (2 appels), messages `CoercionFailed` /
  `PolynomialError` des chemins stock ; le corps stock calcule encore le
  pgcd kkinv du cas 44 (seul le noyau natif lève la `PolynomialError` de
  SymPy) ; surcoût restant de l'admission splitfactor sur QQ (2,8 ms contre
  1,0 ms SymPy, voir l'écart ci-dessus).
- ratint/exroots : restent ratint_logpart sur `QQ<a>[p]` paramétrique et
  divisions non rétractables de `QQ<a>(t)` (repli de référence, erreur),
  `EX(a)` Python qui simplifie son élément, ordre des étapes d'expand hors
  produits de radicaux, constructeur RootSum en Python, bruit imaginaire de
  evalf (approché par l'évaluateur double dans l'ordre d'impression).
- Limites : 82/1293 asserts encore différents (vague 43, lot limits2 ;
  154 avant) ; 84 déclins Gruntz (87 avant), dont 67 oscillants, tous
  repris en Rust (décomposition bornée ou marcheur `AccumBounds` natif,
  vague 41) ou gardés. Restent :
  - `fps` (26 asserts, domaine Python non porté) ;
  - séries à exposants symboliques ou logarithmiques
    (NotImplementedError là où SymPy développe : aseries simple 2, 6, 7,
    bug2, issue 3463, log_power1, 14384:3, 6350) ;
  - puissances de log dans le terme O (3508, 5223:19), `exp(sqrt(x))`
    (exp_sqrt_1), `sin(x)**y` avec `logx` (5223:12), Taylor générique
    d'une sous-classe hôte de `Function` (3978:5), 21938 (ordre plus
    faible mais vrai) ;
  - 30 différences de forme seule (même valeur) ;
  - écarts à valider : 5172 (×2, SymPy rend 0, symox la valeur attendue
    par le test SymPy), 19770:2 (0 contre limite gardée), 7535 (×3,
    `AccumBounds(-oo, oo)` où SymPy garde la limite),
    gruntz_eval_special_slow:2 (SymPy expire), gruntz_eval_special_fail
    (PoleError contre NotImplementedError), order_failing_due_to_solveset
    (SymPy lève TypeError) ;
  - `Limit(x + y, x + y, a).doit()` reste gardé (SymPy rend `a`).
  `test_series_AccumBounds` (×2) garde une autre forme du même intervalle
  (`simplify` rend une autre forme de la somme) ; `heuristics()` est
  exposé (natif, vague 43).
  Restent en Python autour des bornes : `oscillating_integral` (aiguillage
  de `integrate` vers les limites de la primitive), la réévaluation des
  bornes par `subs`/`doit` (`_accum_eval._needs`, `reevaluate`), et les
  `Min`/`Max` à arguments bornes appelés par le natif via le crochet
  `lattice`. Écarts : `oo*sign(AccumBounds(...))` de SymPy rendu
  `AccumBounds(-oo, oo)` (divergence consignée), ordre d'impression de
  `AccumBounds*nan`, LaTeX `\left\langle a, b\right\rangle` absent.
- Matrices : perf lambdarr gamma/poly2.
- Tableaux : 19/201 tests SymPy de tableaux ne passent pas (vague 41 ;
  36 avant). Les causes sont :
  - `MatrixSymbol.args[0]` en `str` (le `Str` hôte active les chemins hôtes
    lents) : dérivées de traces, `matrixelement_diff`, `matexpr_subs` ;
  - l'ordre de `Sum` (`vectors_and_scalars`) ;
  - les postprocesseurs MatMul/MatAdd de `Mul`/`Add` (aiguillage Python de
    `_compat.py`, logique à porter en Rust) ;
  - `applyfunc` (2), `matadd_simplify`, `simplify_matrix_expressions`,
    `MatrixSet`, numpy, Hadamard `check=False`, `Trace(A).arg is A`,
    indexation symbolique des tableaux immuables ;
  - les écarts délibérés à valider : `combine_kronecker` sans la règle de
    somme ni les puissances non entières (fausses), donc `combine_mul` et
    `combine_pow`.

  À cela s'ajoutent le repli vendu de `convert_array_to_matrix` (6/527),
  l'ordre des arguments d'`ArrayDiagonal` (arrindex u_int_mul) et le coût
  de `matrix_derivative` par le chemin tableau.

Procédure de validation (leçons du 30 octobre et du 5 novembre) : lancer la
suite avec `XLA_PYTHON_CLIENT_PREALLOCATE=false` (sinon
`tests/extra/test_jax_codegen.py::test_custom_functions_and_standalone_source`
échoue dans le shard 23 par CUDA_ERROR_OUT_OF_MEMORY du sous-processus jax,
le pytest parent ayant préalloué la mémoire GPU ; échec d'environnement,
reproduit avec le .so de main) ; n'accepter que "SUITE CLEAN" (le compte de
lots "passed" incluait des lots en échec) ; le build doit sortir en succès
avant de copier `lib_core.so` (un build tué en silence laissait l'ancien
`.so`) ; doctests parallèles à relancer en mono-thread ("Couldn't
compile") ; tests Python en worktree avec `PYTHONPATH=python` ; jamais de
`pkill` large. Le compte de base « 60 020 réussis » de main comptait
15 tests sautés de plus (« Slang cuda device unavailable », 60 sautés contre
45) : les lots se comparent sur les tests collectés et les non-réussis
(main 60 085 collectés ; vague 39 : 51 non réussis, 6 xfail dont le xfail
strict prde).

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
| Corpus EDO, rejeu du 5 novembre (`feat-ode-20261105`) | 315 / 287 / 268 identiques sur 315 (main `15e75d3a` : 315 / 285 / 267) | [Preuves](audits/performance/2026-11-05-ode/README.md) ; fact_12 et liouville_08 corrigés, aucune régression |
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
- [ ] Limites (après le lot limits du 5 novembre) : porter en Rust le
  marcheur oscillant `AccumBounds` (68 déclins « series in mrv variable ») ;
  formes de valeur du chemin méromorphe (`test_rational`, 5740, 6364,
  14393, 16708, 17751), `series_AccumBounds`, 21785 ; variable non symbole
  (9205) ; `re`/`im` (21756) ; lever `Not sure of sign of 1/a` pour
  `Mod(x, a)` en `a` (25230) ; `x**n` et `2*(x + 1)**w` en `-oo` (SymPy
  `(-oo)**n`, `2*(-oo)**w`) ; `I**x` en `oo` (SymPy `nan`),
  `(3/x)**w` en `0+` (SymPy `oo**w`), `sin(x)**w` en `0+` et
  `x**a/(x**a + 1)` en `oo` (SymPy lève) ; `mrv_leadterm` rend le
  coefficient développé (`b/x - 1/x` pour `(b - 1)/x`) ; `test_aseries_trig`. Ces cas restent non
  évalués, jamais faux.
- [ ] Réconcilier les derniers indices et formes du corpus EDO avec les
  correctifs ultérieurs des solveurs et de `simplify`. Les solutions
  séparables/Riccati auparavant implicites sont maintenant résolues ;
  les différences de forme seules ne sont pas des solutions manquantes.

Preuves : [limrest2](docs/algorithmic-parity.md#limrest2-wave-14),
[limits](audits/performance/2026-11-05-limits/README.md),
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
