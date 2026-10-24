# TODO - Capacités mathématiques et parité SymPy

État documenté le **22 octobre 2026**, après la vague 3. Référence : **SymPy 1.14.0**.

Objectif : atteindre au moins les capacités mathématiques de la référence,
en conservant les domaines acceptés, l'arithmétique exacte, la complétude des
solutions, la précision demandée et le comportement en cas d'échec.
La présence des noms dans l'API et les performances sur quelques exemples ne
constituent pas une preuve de parité.

Ce fichier couvre les algorithmes mathématiques et leur validation. Le chantier
des extensions et sous-classes Python est suivi séparément et n'est pas inclus.

Contrat détaillé : [parité algorithmique et écarts ouverts](docs/algorithmic-parity.md).
Les priorités ci-dessous sont un ordre de travail proposé, pas une estimation
quantifiée de l'importance de chaque écart.

## Écarts reproduits dans l'environnement local

Comparaisons réalisées avec `.venv/bin/python`, Symox 0.1.1 et SymPy 1.14.0.
Les symboles `t` et `s` du test de Laplace sont déclarés positifs ; `f` désigne
une fonction inconnue de `x` dans les tests d'EDO. La colonne d'état décrit la
branche `integrate-wave1-20261021`.

| Calcul | SymPy 1.14.0 | Symox avant la vague 1 | État |
| --- | --- | --- | --- |
| `QQ.algebraic_field(2**Rational(1, 4))` | Construit le corps exact | `ValueError` : extension hors des modèles pris en charge | Fermé |
| `integrate(exp(-x**2), x, risch=True)` | Renvoie un `NonElementaryIntegral`, avec preuve de non-élémentarité | `NotImplementedError` pour l'option `risch` | Fermé |
| `laplace_transform(cos(sqrt(t)), t, s)` | Calcule une expression explicite | Conserve un `LaplaceTransform` non calculé | Fermé (calculé, forme différente) |
| `integrate(t*exp(-s*t), (t, 0, oo))`, `s` sans hypothèse | `Piecewise((s**(-2), Abs(arg(s)) < pi/2), (Integral(...), True))` | **Valeur fausse** : `(-oo/s - 1/s**2)*exp(-oo*s) + s**(-2)` | Fermé |
| `integrate(exp(-s*t), (t, 0, oo))`, `s` sans hypothèse | `Piecewise((1/s, Abs(arg(s)) < pi/2), (Integral(...), True))` | **Valeur fausse** : somme de `Piecewise` contenant `exp(-oo*s)` | Fermé |
| `integrate(exp(-t)*exp(-s*t), (t, 0, oo))` | `Piecewise((1/(s + 1), Abs(arg(s)) <= pi/2), ...)` | Conserve l'intégrale non calculée | Fermé |
| `dsolve(f.diff(x, 2) - x*f, f)` | Solution en `airyai` et `airybi` | `NotImplementedError` | Fermé |
| `dsolve(x**2*f.diff(x, 2) + x*f.diff(x) + (x**2 - 4)*f, f)` | Solution en `besselj` et `bessely` | `NotImplementedError` | Fermé |
| `pi.evalf(30, strict=True)` | Calcule la valeur | `NotImplementedError` pour l'option `strict` | Fermé |
| `pi.evalf(30, maxn=100)` | Calcule la valeur | `NotImplementedError` pour l'option `maxn` | Fermé |
| `meijerg([[], []], [[Rational(1, 3)], []], Rational(2, 3)).evalf(50)` | Donne une valeur à 50 chiffres | Conserve l'expression symbolique | Fermé |

Ces observations établissent des écarts concrets sur la version locale testée.
Elles ne remplacent pas un audit de tous les domaines et options, ni une nouvelle
validation après modification du moteur.

## Bilan de la vague 1 (21 octobre 2026)

Trois branches fusionnées dans `integrate-wave1-20261021` :
`feat-integration-convergence-20261020` (24 commits),
`feat-algebraic-fields-20261020` (4 commits) et
`feat-ode-precision-20261020` (11 commits).

- Intégration : corpus de référence passé de 988 identiques / 91 manquants à
  1 033 identiques / 54 manquants. Les deux seules valeurs numériques
  différentes sont des erreurs de la référence (`laplace_transform(1/(t+1))`,
  transformée de Fourier de `1/(1+x**2)` pour `k < 0`), vérifiées par
  quadrature mpmath.
- Corps algébriques : général, plongements, `CRootOf` certifié.
- EDO et précision : modules d'EDO et d'EDP de SymPy 1.14 portés, solveurs
  natifs conservés comme repli. Corpus de 315 EDO de la suite de SymPy :
  classification identique pour 295, indice du test pour 237, `dsolve` par
  défaut pour 245.

Reruns downstream sur la build release fusionnée, sans régression :

| Projet | Avant la vague 1 | Après |
| --- | --- | --- |
| Lcapy | 305 réussis / 12 échecs (référence 317) | 314 / 3 |
| Sage `sage/interfaces/sympy.py` | 311/319 | 317/319 |
| Torch `test_sympy_utils.py` | 34 échecs | 1 échec (`SingletonInt._eval_Eq`) |
| Optlang | 346 | 346 |
| EinsteinPy | 287 / 9 xfail | 287 / 9 xfail |
| Cadabra `sympy_cdb.cdb` | 12/12 | 12/12 |
| OctSymPy | 131/131 | 131/131 |

Journaux : `~/.cache/symox-validation/wave1/<projet>/`.

Changements de comportement alignés sur la référence, à connaître :

- `2nd_power_series_ordinary` refuse désormais les équations avec second membre
  (`ValueError`) et ignore `ics`, comme SymPy. Symox résolvait ces équations
  auparavant. **Décidé le 22 octobre 2026** : rester
  aligné sur SymPy.
- `sysode_linear_2eq_order1` lève `UnboundLocalError`, bug de SymPy reproduit.
- Racines indicielles de Frobenius répétées : Symox renvoie une base avec terme
  logarithmique, là où SymPy ne renvoie que la première série. C'est une
  extension, pas une correspondance.

## Bilan de la vague 2 (22 octobre 2026)

Dix-huit branches fusionnées dans `integrate-wave2-20261021` (au-dessus de la
vague 1, non poussée, main non avancé) : nested-rings, transforms,
solve-forms, frobenius, solve-gaps, misc-parity, lcapy-fixes, noncommutative,
gaussian-domains, domains-misc, piecewise-periodic, simplify-perf,
integrate-options, union-sets, solve-sets, inexact-composite, pde-perf,
ratint, plus des correctifs d'intégration (mpmath hors de l'import,
`AppliedUndef.fdiff` : un argument symbole n'est tenu directement que si aucun
autre argument n'en dépend, ce qui supprimait un débordement de pile natif dans
`pdsolve`).

Validation finale (`7bac15d`) : `cargo test -p symox` 569 + 5 réussis, suite
Python complète 36 947 réussis, 0 échec, clippy `correctness` propre.

Reruns downstream (build debug, `~/.cache/symox-validation/wave2/<projet>/`) :

| Projet | Après la vague 1 | Après la vague 2 |
| --- | --- | --- |
| Lcapy | 314 / 3 | 317 / 0 |
| Sage `sage/interfaces/sympy.py` | 317/319 | 319/319 |
| Torch `test_sympy_utils.py` | 1 échec | 1 échec (`SingletonInt._eval_Eq`, inchangé) |

Régressions trouvées par la suite et corrigées avant fusion finale : blocage de
`pdsolve` dans `heurisch` (`find_non_syms`), test `checkpdesol` passé de plus de
300 s à 41 s (`Integral.doit` et `Derivative.doit` comme SymPy, cache de
`count_ops`), `solve` sur une `Matrix` avec flottants, `solve` sur une
expression Sage non sympifiée.

**Décisions utilisateur du 22 octobre 2026 :**

- main avancé sur `integrate-wave2-20261021` (local, non poussé). Le `dsolve`
  de `(x + cos(f))*f' + f + 2*x` reste un écart de performance : 395 s en
  release contre 82 s pour SymPy (qui lève ensuite `NotImplementedError`).
- Corrections intentionnelles validées : résultats justes conservés et
  regroupés dans `docs/algorithmic-parity.md` (section "Intentional
  corrections validated on 2026-10-22").
- `2nd_power_series_ordinary` reste aligné sur SymPy (pas d'extension).

## Bilan de la vague 3 (22 octobre 2026)

Vingt branches fusionnées dans `integrate-wave3-20261022` (au-dessus de
`0a329a8`, non poussée). Premier lot : domains-inexact, sets-reals,
definite-params, domains-inexact2, trig-engine, evalf-cancel. Second lot :
integrate-rest, ode-float-corpus, piecewise-series, transforms-rest,
gaussian-fractions, sets-rest, solve-forms, noncomm-rest, sets-rest2,
sets-rest3, core-misc, forms-rest, heaviside-limit, extension-root.

Résultats faux trouvés et fermés : `limit(Heaviside(x), x, 0)` donnait `1/2`
dans toutes les directions, `Pow(2, Rational(10**20 + 1, 2))` donnait
`sqrt(2)`, `1/(cos(x) + 2)` sur `(0, 2*pi)` donnait `0`, l'intersection de
`ImageSet(n*I + 1)` avec `[0, 1]` donnait `EmptySet` (faux aussi dans SymPy),
régression flottante du corpus EDO introduite par `efa251c`.

Validation finale (`29d6bf3`, plus `394620b` qui met à jour un test Rust
écrit pour les anciennes limites de `as_real_imag`) : `cargo test -p symox`
569 + 5 réussis, clippy `correctness` propre, suite Python complète 38 669
réussis, 0 échec ; corpus EDO 296/263/259 (classification, indice, défaut) ;
corpus d'intégration 0 faux, 0 dépassement de temps.

Toutes les différences conservées où SymPy se trompe sont listées comme
corrections intentionnelles dans les items ci-dessous et dans
`docs/algorithmic-parity.md`.

**Décisions prises le 22 octobre 2026 :** implémenter `AccumBounds` comme
SymPy (`sin(oo)`, limites oscillantes, les 14 intégrales oscillantes) ; garder
le commit de documentation `1b42d9a` ; critère d'arrêt des items ouverts : le
TODO est fini quand le tri des 1 363 objets aux observations différentes de
l'audit versionné (`docs/sympy-parity-audit.md`) ne laisse que des différences
de forme, des corrections intentionnelles ou des lacunes documentées comme non
prises en charge.

main avancé par fast-forward sur `394620b` le 22 octobre 2026 et poussé sur
origin. Vague 4 lancée (branche d'intégration `integrate-wave4-20261022`) :
accumbounds, nroots-precision, domains-rrcc, audit-triage.

## Vague 4 (22 octobre 2026) : tri de l'audit versionné

Rejeu complet des exemples de documentation (2916 objets, empreinte
`dfb31b50...`) : 1749 identiques (1551 avant), 1163 différents, 2 incomplets
(délais natifs `Poly.which_all_roots`, `stats.coskewness`). Tri dans
`audits/sympy-1.14.0/triage-wave4.json` et `docs/sympy-parity-audit.md` :
312 forme, 137 équivalents, 557 fonctionnalités manquantes, 3 résultats
incomplets, 4 précisions insuffisantes, 45 erreurs différentes, **71 résultats
faux**, 15 cas où SymPy se trompe, 23 non déterministes.

Selon le critère d'arrêt, il reste à corriger les catégories g, e et d, puis
à implémenter ou documenter comme non prises en charge les 557 de la
catégorie c. Lots d'environ 1 h :

- [x] WP-G1 nombres et hypothèses (`zoo + zoo`, `comp`, `all_close`,
  `check_assumptions`, `unchanged`, `WildFunction.match`) et WP-G2 variables
  liées et `Factors` (`edc2cee`, `d32e433`, `0f4b474` : `nargs` propre à
  chaque WildFunction).
- [x] WP-G3 fonctions spéciales (`genocchi`, `rf`/`ff` de `Poly`, `Ei` polaire)
  et WP-G4 sommes, limites, intégrales, ensembles (`cbfe247`).
- [x] WP-G5 solveurs, WP-G6 matrices (`d78539d`), WP-G7 géométrie et
  vecteurs, WP-G8 combinatoire et holonome (`0631ef6`).
- [x] WP-G9 à G11 physique quantique et seconde quantification (`04ec521`,
  `d336bce`, `05990e8` ; sondes 161/235 → 234/235, secondquant 6/88 → 76/88),
  WP-G12 mécanique et WP-G13 contrôle et milieux continus (`16b4b68` ; 45/149
  → 135/149). Restes documentés : `qapply` d'un produit tensoriel
  d'opérateurs, coordonnées de repère en `Dummy`, `substitute_dummies`.
- [x] WP-G14 optique, Wigner, stats, impression (`59fc3a7`, mergé dans
  l'intégration). Restes : `P(Eq(A, 2) & Eq(B, 3))` de variables de Poisson
  indépendantes renvoie encore un `Piecewise` faux (à traiter avec C10) ;
  `Tuple(1, 2) == (1, 2)` faux dans symox, vrai dans SymPy ; `fcode` des
  nœuds codegen.
- [x] WP-E1 précision (`expm1`/`log1p` de `1e-99`, `lambdify` mpmath, `gbt`)
  et WP-D1 résultats incomplets (`substitution`, `get_contraction_structure`,
  `MusculotendonDeGroote2016.p`) : `089392a`, `e6a481d`, `aae930a`. Reste :
  `solveset` complexe natif incomplet sur `exp(x) - exp(-x)`,
  `exp(2*x) - 3*exp(x) + 2`, `2**x - 8`, `sin(x) + cos(x)` (Rust).
- [x] Valider les corrections intentionnelles (h) : les 15 de départ et
  celles des vagues 10 à 13 (symox juste, SymPy faux), vérifiées
  numériquement et inscrites dans `docs/algorithmic-parity.md` ; validées
  par l'utilisateur le 24 octobre 2026.
Étape validée le 22 octobre 2026 et main avancée à `91edb6a` : 39 388 tests,
corpus ODE 296/263/259, corpus d'intégration 1105/75/18/10/2 (0 faux).

- [x] Régression de `solve([x**2 + y**2 - 4, y - x**4 + 1], [x, y])`
  (`5cecd40`, mergé ; Rust) : plus de 200 s → 4-5 s (SymPy 12 s), mêmes 8
  solutions exactes. Cause : `expand(complex=True)` sans borne sur des
  `sqrt` de sommes radicales positives (signe décidé numériquement comme
  SymPy), vérification numérique superflue des racines (SymPy :
  `check=False`), `ln 2` recalculé. Reste : résolution par substitution au
  lieu de Groebner, résultats bien plus gros (1 Mo contre 16 ko).
- [x] WP-C4/C6 domaines et corps de nombres (`126f84a`, mergé) : `round_two`,
  `prime_decomp`, `prime_valuation`, modules, `integral_basis`,
  `galois_group`, `GF(p)`, méthodes `ZZ_I`/`QQ_I`. Restes : groupes de
  Galois de degré 5 et 6, alias `alpha`/`zeta` à l'affichage, `GF(n)` composé
  (SymPy l'accepte comme anneau Z/nZ, symox le refuse déjà sur main).
- [x] WP-C1/C2 fonctions spéciales (`a61b274`, mergé ; Rust) : sondes
  127/231 → 24/231 écarts. Restes : séries d'Airy et elliptiques,
  `dirichlet_eta(s, a)`, `f(nan)` évalué en `nan` au lieu de rester tel quel,
  quelques valeurs en `±I*oo` liées à `exp(I*oo)`.
- [x] WP-C3/C9 appariement, séries, sommes infinies, ensembles, matrices
  (`e7e1006`, mergé) : sondes 73/129 → 132/139. Restes :
  `solve_triangulated(extension=True)` (NotImplementedError), série d'`Abs`
  non factorisable, `expand()` ne sépare pas `(2*x)**y`, complément et unions
  mixtes de `ComplexRegion`.
- [x] WP-C5 API `DomainMatrix`/`DDM`/`SDM` (`1e72749`, mergé) : 85/85
  enregistrements de l'audit (28 avant), 984/984 sorties de docstrings, sans
  perte de performance. Restes : `DFM` (sans python-flint, comme la
  référence installée), module `eigen` (manque `FiniteExtension`).
- [x] WP-C8 géométrie (`6309214`, mergé) : 234/236 sondes (27 avant) ;
  moments d'aire, `reflect`, intersections exactes de coniques, plans.
  Restes : `Plane.parameter_value` à un paramètre peut différer de 2π
  (solveur de système incomplet), `reflect` d'ellipse sur droite oblique.
- [x] WP-C10 stats (`2960054` à `79e9ceb`, mergés) : probabilités jointes de
  variables indépendantes corrigées (le `Piecewise` faux), `piecewise_fold`
  ne sort plus une condition liée d'une `Sum`/`Integral` (faux aussi dans
  SymPy), `rewrite`/`expand`, `given` avec Tuple, solitons symboliques,
  `coskewness` 60 s → 0,14 s, chaînes de Markov symboliques. Restes : états
  `Str`, `canonical_form`, `communication_classes`, `decompose`,
  `fundamental_matrix`.
- [x] WP-C7 groupes (`cecc913`, mergé) : doctests SymPy 1145/1588 →
  1534/1588 exemples, comparaisons de valeurs 22/98 → 96/98 ; groupes de
  présentation finie, tables de cosets, homomorphismes, systèmes de
  réécriture, `pc_groups` (SymPy y est faux sur `induced_pcgs`). Restes :
  `_holder_formula`, énumération B(2,4) lente, repr récursif de `Collector`.
- [x] WP-C11 quantique et secondquant (`24f20db`, mergé) : sondes 42/126 →
  114/126 ; `represent` complet, commutateurs `doit`/`expand`, `couple`/
  `uncouple` à 3 spins et plus, module `pauli`, `hbar` numérique, mer de
  Fermi dans `apply_operator`, `substitute_dummies` complet. Restes :
  `couple` à spins symboliques, `XOp` dans la base `PxKet`, `KroneckerDelta`
  à travers le niveau de Fermi, `Sum` d'opérateurs.
- [x] WP-C13/C15 contrôle, quaternions, tenseurs (`88cf3ee`, mergé) :
  sondes contrôle 24/61 → 58/61 (`lti.py` de référence porté, MIMO,
  sous-espaces contrôlables/observables), quaternions 3/171 → 170/171,
  tenseurs 142/199 → 178/199. Corrections où SymPy est faux :
  `from_rotation_matrix` pour une rotation de pi, `replace_with_arrays` de
  `PartialDerivative`. Restes : `TransferFunctionMatrix` depuis un `Tuple`,
  `StateSpace.dsolve`, expressions de tableaux simplifiées, `Integral`
  d'un tableau, pprint des objets de contrôle.
- [x] WP-C12 mécanique (`960bbf0`, mergé ; Rust) : sondes 137/244 →
  215/244 ; `physics.vector`/`mechanics`, `beam`, `truss`, `arch` de
  référence portés (arguments complets des joints, `System.form_eoms`,
  `JointsMethod`, rotules et lignes d'influence de poutre, torsion de
  `Beam3D`) ; coordonnées de repère en `CoordinateSym` (avant, `A[0] + B[0]`
  fusionnait faussement en `2*A_x`) ; `lambdify` de `SingularityFunction`.
  `Cable` reste natif (équilibre faux dans SymPy). Restes : repr des tracés,
  `latex(Matrix)` général lève `TypeError` (déjà sur main).
- [x] WP-C14 génération de code et imprimeurs (`b31e33a`, mergé) :
  imprimeurs Fortran, Rust, GLSL, JS, Julia, Octave et R portés (Matrix,
  `assign_to`, Piecewise) ; `cfunctions` 1/14 → 14/14 en comportement ;
  nœuds `ast`, `fnodes`, `cnodes` complets. `rust_code(x + 1)` imprime
  désormais `x + 1` comme SymPy (avant `1.0`). Corrections : `do concurrent`
  Fortran 2008 valide, déclarations Fortran 77 simples. Restes : repr des
  nœuds codegen, arithmétique d'`Idx`, smtlib, `codegen.rewriting`.
- [x] WP-C1 à C15 : toutes les fonctionnalités manquantes traitées (restes
  listés par lot ci-dessus).
- [x] Écarts du noyau (`d5a9bbf`, `f584c0d`, mergés ; Rust) : `Tuple` égal
  et de même hachage qu'un tuple Python, règles d'assumptions de SymPy sur
  sommes, produits, puissances et `exp` avec `I` (`(2+I).is_complex` vrai,
  `(oo+I).is_extended_positive` corrigé : il valait faussement True),
  `GF(n)` composé comme anneau Z/nZ. Restes : `sqrt(2+I).is_extended_real`
  inconnu (SymPy False).
- [x] `Poly` sur Z/nZ composé (`66dd7cc`, mergé) : 205/228 sondes contre 8
  (polys denses de référence vendorisés dans `polys/_dense/`, graine fixe
  pour la factorisation aléatoire ; chemin premier 1-6 % plus lent).
- [x] `solveset` complexe natif (`b65ab71`, mergé ; Rust) : familles
  `ImageSet` complètes sur exp, trigo et hyperbolique (sondes 40/78 → 61/78
  identiques, 8 autres égales à la forme près), `ConditionSet` au lieu d'un
  ensemble fini incomplet (avant : faux sur `x**x - 1`, `exp(exp(x)) - 1`,
  `exp(x) - exp(-x)`). Corrections où SymPy est incomplet : `2**x - 3**x`
  (SymPy `{0}`), `log(x) + log(x - 1)` (SymPy garde une fausse racine).
- [x] `solve(sin(x) + cos(x))` (`4f80d00`, mergé ; Rust) : expansion
  d'angles, réécriture en tan et repli exp de SymPy portés ; 65/77 sondes
  contre 28, `solveset(sin(a*x) - 1)` en `ConditionSet(Ne(a, 0), ImageSet)`.
  Restes : `exp(exp(x)) - 1` (ImageSet sur `ProductSet`), `sec(x) - 3`,
  coefficient symbolique avec plusieurs fonctions trigo.
- [x] `limit_seq` imbriqué (`9320c8d`) : algorithme de référence porté
  (`dominant`, boucle `trials`, parité de l'indice) ; 131/159 sondes contre
  111, aucun résultat faux, les restes sont des `None` (nombres harmoniques
  et asymptotique de gamma à l'infini).
- [x] `limit(factorial(n)*exp(n)/(n**n*sqrt(n)), n, oo)` (`de932a6`,
  `484b6c7`, mergé ; Rust) : réécritures gamma de gruntz portées, séries de
  Stirling et polygamma à l'infini, `Limit` non évalué au lieu d'un 0 faux
  pour `Ei`, `erfi`, `besseli`, `fibonacci`... Le segfault ne se reproduit
  pas. Restes : 42 sondes non évaluées (jamais fausses).
- [x] `limit_seq((n**2+a)**(n+b)/(n+a)**(2*n), n)` donnait 0 (SymPy `None`) :
  corrigé dans gruntz (`211c5a5`, mergé).
- [x] `sqrt(-sqrt(2))` non auto-évalué en `2**(1/4)*I` (`737cf78`, mergé ;
  Rust) : puissance d'un produit à facteurs négatifs découpée comme
  `_eval_expand_power_base` ; corpus ODE revenu à 296/263/259, sondes
  1232 → 1588/1669 (le reste ne diffère que par la forme). Corrige aussi un
  blocage de `solve_poly_system([x**4+x+1])`.
- [x] Ordre canonique d'`Add` (`d339ff8`, mergé ; Rust) : l'égalité entre
  un nœud tenu (`evaluate=False`) et un nœud canonique dépendait de l'ordre
  de hachage, donc de l'ordre d'import (vrai par hasard pour 209/870 paires
  de symboles, SymPy 0) ; `AccumBounds` en coefficient de tête comme dans
  `Add.flatten`. Restes : `x + Ket('k')` différent de `Add(Ket('k'), x)`,
  `subs`/`expand` ne fusionnent pas les nombres dans un `AccumBounds`.
- [x] Re-triage de l'audit (`7d071c6`, `audits/sympy-1.14.0/triage-wave4-after.md`)
  sur 1167 objets : 587 résolus. Avant → après : g 71 → 7, e 4 → 1,
  d 3 → 0, c 557 → 155, f 45 → 41, h 15 → 19, b 137 → 128, a 312 → 199.
- [x] 7 g et 1 e restants : `Integral` à `t` libre, `Range.inf`/`sup`
  inversés, `solve` sur `IndexedBase` qui rend `[]`, `lambdify` qui passe
  les entiers en flottant (lot A) ; `Dagger` de non-commutatifs sans
  inversion, `A*X` de matrices aléatoires réordonné, `Expectation` imbriqué,
  `nonlinsolve` rationnel (lot B). Corrigés (`d4dea76`, `c7659d4`, mergés).
- [ ] 155 c restants, en 4 lots. c1 (stats, holonomic, sommes, séries,
  intégrales ; `211c5a5`, mergé) : 19/35 corrigés, 3 corrections
  documentées, 12 lacunes (processus stochastiques, `to_meijerg`, `fps`
  symbolique). c2 (functions, core, sets, solvers, simplify ; mergé) : 35/41,
  6 lacunes (forme simplifiée de `Derivative`, `nseries` de `x**y`,
  `trigsimp(old=True)`, `devise_plan_meijer`, `solve_triangulated`
  avec extension). c3 (matrices, polys, printing, codegen ; `8f97520`) : 29/41,
  `Indexed` corrigé au niveau de l'impression seulement, `combine_kronecker`
  sans la règle de somme (SymPy rend à tort `KroneckerProduct(A + C, B + D)`,
  candidat h). c4 (physics, tensor, vector, categories ; `d8092c1`) : 12/37,
  25 lacunes, surtout `pprint` d'objets définis en Python ; étape
  `DiracDelta` gardée dans l'intégration, `orient_new` sans contrôle
  d'orthogonalité comme SymPy.
- [x] 41 f (type d'erreur, `2e29090`) : 39/39 classes d'exception
  identiques, 32 messages identiques ; `solve_poly_system(strict=True)` lève
  `UnsolvableFactorError`.
- [x] Formes stats (47 b, `d68587a`) : 19 identiques, 8 en partie ; le
  reste vient d'une autre méthode (moments fermés + simplify au lieu de
  l'intégration) ou de divergences gardées (garde de support de LogCauchy,
  `MatrixStudentT` exact, ensembles circulaires évalués).
- [x] Impression des expressions matricielles (`a36f497`) : écarts `str`
  17 → 2 sur 82 sondes, doctests `schur`/LU de référence passent.
  Restes : `g(mat)` en julia/octave (argument matriciel non supporté par le
  noyau). latex/pretty matriciels (`4a73615`, mergé) : écarts 121/94/117 →
  2/2/2 sur 126 sondes.
- [ ] 167 a (impression) en 2 lots. a1 (polys, core, functions..., `da91415`,
  mergé) : 40/90 corrigés, le reste documenté (9 ordre seul, 40 lacunes :
  `DDM`/`SDM` qui affichent des MPQ, pretty des ensembles, `hyper`,
  `meijerg`, alias d'`AlgebraicField`, `Add` de `MatrixSymbol`). a2
  (physics, diffgeom, combinatorics, tensor, `60ce230`, mergé) : 57/77
  (contre 8), 14 ordre/forme, 1 erreur de référence (`wigner_9j`), 5
  lacunes (xypic, backends matplotlib, `Reshape`). print-a3 (`e58512f`, mergé) :
  pretty/latex des ensembles, `Order`, `hyper`, `meijerg`, `Integral` ;
  `DDM` en `MPQ` gardé (SymPy y garde des int bruts, d'où des flottants dans
  `rref`). print-order (`58e437e`, mergé ; Rust) : ordre des termes et
  facteurs à l'impression selon `as_ordered_terms`/`as_ordered_factors` et
  `default_sort_key` (tables `class_key`) ; sur 1973 sondes, écarts `str`
  182 → 138, "ordre seul" 57 → 17 ; latex/pretty d'un polynôme de 200 termes
  environ 12 % plus lents. Pretty de `Sum`/`Limit` porté en Rust
  (gaps-misc, `20e7154`, mergé). Agent en cours : small-rest (ordre
  `Sum`/`Limit` dans un `Add`, `re`/`im` non commutatifs, `hbar`, pretty
  d'`Integral`).
- [x] Formes équivalentes, lot b2 (`4b4fd1c`, mergé) : 41/81 identiques
  (contre 8), corpus inchangés.
- [x] Derniers chiffres d'evalf (`5d39e69`, mergé ; Rust) : 2796 sondes
  (60 fonctions, 15/30/50 chiffres), exactes 2641 → 2680, erreurs 53 → 26 ;
  `beta`, `rf`, `ff`, `lowergamma`, `laguerre`, `hermite`, `assoc_legendre`
  identiques à SymPy. Gardés : cas où symox est correctement arrondi et pas
  SymPy (vérifiés avec mpmath). Restes : `elliptic_pi` (7 sondes au-delà de
  10 s), `jacobi`/`assoc_laguerre` de degré négatif non entier.
- [x] 7 lacunes diverses (gaps-misc, `20e7154`, mergé ; Rust) :
  `reduce_inequalities` à un symbole d'intérêt (8/9), `wrt`/`sort` dans 21
  fonctions polys (134/134), `integrate` de vecteurs, `lambdify` qui garde
  les entiers numpy, impression `ANP`. Divergence : `Add(A, 1)` d'une
  `MatrixSymbol` lève `TypeError` (SymPy construit `1 + A`).
- [x] Restes c1/c2 (c-rest-a, `25b75f5`, mergé) : `devise_plan_meijer`,
  `nseries` de `x**y`, `trigsimp(old=True)`, `fps` symbolique, densités
  stats symboliques, `BernoulliProcess`. 3 erreurs de la référence gardées
  (candidats h). Regroupement `2**n*3**n` → `6**n` porté ensuite par
  numbase (voir plus bas).
- [x] Restes c3/c4 (c-rest-b, `7fec210`, mergé) : c4 12 → 23/37 (pprint des
  objets Python avec le `stringpict` de SymPy, catégories et dessin de
  diagrammes portés, `cancel(expand=False)`), c3 inchangé à 33/41 (reste :
  `galois_group` degré 6, `prime_decomp` symbolique...).
- [x] Segfault `erf(B*A - 2*x)` avec non-commutatifs (nc-segv, `2c7c8b2`,
  mergé ; Rust) : boucle de `could_extract_minus_sign` ; 1248 sondes, 88
  crashs → 0. `DiracDelta`/`Heaviside` lèvent la `ValueError` de SymPy pour
  un argument non réel.
- [x] Test `li` instable : `mp.dps = 30` au niveau module dans
  `tests/test_special_numeric.py` fuyait vers les autres tests (`742cbc0`).
- [x] stats-rest (`86cca9f`, mergé) : doctests stats 137 → 163/171
  (symboles aléatoires imprimés par leur nom, conditions tuples dépendantes,
  Laplace multivariée, `stationary_distribution`, chaînes de Markov à temps
  flottant ou symbolique, Wiener/Gamma avec `Contains`, forme close GSE).
  Divergences gardées (SymPy faux) : probabilité d'union Wiener = 1,
  conditionnelle Poisson `45*exp(-3)/(2*(-4+exp(3)))` (SymPy `5*exp(-3)`),
  valeurs de `E` conditionnelles. Reste : `MultivariateNormal` avec
  `MatrixSymbol` symbolique (NotImplementedError), impression de
  `(T**5)[1,2]`, ordre `X + x`, forme de la cdf de Moyal.
- [x] small-rest (`2005b4d`, mergé ; Rust) : ordre de Sum/Limit/Subs dans
  un Add 19 → 2, `re`/`im`/`Abs`/`arg`/`sign` d'arguments non commutatifs,
  `hbar.is_number`, pprint 2D des conteneurs Python, options `greedy` et
  `symmetric` de `div`, `ValueError` des polynômes orthogonaux de degré
  négatif. Divergence : degré positif non entier laissé non évalué (SymPy
  tronque à l'entier).
- [x] numbase (`f9981f5`, mergé ; Rust) : règle `num_exp` de SymPy portée
  (`2**n*3**n` → `6**n`) sans casser `meijerint` (étape gammasimp/combsimp) ;
  sondes 54/86 différences → 0.
- [x] holo-rest (`c733938`, mergé ; Rust) : holonomic 8 → 37/40, dérivées
  multiples 4 → 13/13, `solve_triangulated` avec extension 0 → 6/12,
  `sec`/`csc(k*x)`. Changement de politique à valider : `nonlinsolve([a*x], x)`
  rend `{(0,)}` générique comme SymPy au lieu de refuser. Divergence gardée :
  racines annulant un dénominateur retirées (SymPy liste des non-solutions).
- [x] ellippi (`7285fd0`, mergé) : `elliptic_pi(2, 7/2).evalf(50)` 120 s →
  0,18 s (intégrale initiale contrôlée dans l'`ellippi` de mpmath) ; 342
  évaluations, 0 timeout, jamais moins de chiffres justes que SymPy (SymPy
  en a 23 fausses).
- [x] Correctif `hbar` en pretty (`5d7ebd5`) ; main avancé à `59e4d2b`.
- [x] nc-rest (`0129b89`, mergé ; Rust) : sondes non commutatives 174 → 31
  (29 = RecursionError de SymPy, 2 = impression `(-I)` de csc/csch) ;
  commutativité des têtes que SymPy déclare commutatives, faits
  d'hypothèses impliquant `commutative` à False, `Eq(A*B, B*A)` non évalué
  comme SymPy, `(M**5)[1, 2]` parenthésé.
- [x] Re-triage final (`e0fe187`, `audits/sympy-1.14.0/triage-wave4-final.*`) :
  580 → 256 objets encore différents ; g 7 → 0, e 1 → 0, d 0, c 155 → 30
  (5 non documentés), f 41 → 15, h 19 → 21, b 68, i 33, a 89. Main avancé à
  `7589b91`.
- [x] Dernière vague : c-final (`4663d57` : `harmonic`/`loggamma`
  rationnels, série de `loggamma` en 0, équations radicales et
  `implicit=True` de `solve`, `reduce_inequalities` multi-symboles,
  `lambdify.__globals__`, `autowrap` à un helper), c-doc (`1a635fa` : 15/24
  objets c documentés alignés, tensor/array, couple/uncouple, categories,
  diffgeom, `simplify_gpgp`, `AlgebraicField` symbolique,
  `rational_parametrization`), f-final (`1fce7e8` : les 15 f lèvent le type
  d'exception de SymPy, 14 avec le message exact ; `Rational(pi)` lève
  `TypeError` ; tolérance et mesure de résidu de `nsolve` de SymPy).
  Validée le 23 octobre 2026 (8 shards, cargo test, clippy, corpus EDO
  296/263/259, corpus d'intégration 1105/75/18/10/2) ; main avancé à
  `793b98b`. Critère 4A atteint sur l'audit : plus de g/e/d, les c restants
  sont documentés, les f ont le type d'exception de SymPy ; restent des
  formes (a/b), des corrections h à valider et l'environnement (i).
- [x] Vague 5 (chantiers de fond des sections 1 à 4), mergée :
  algfield (`b0b799e` : sondes 158 → 166/177 et 53 → 57/59 ; `CRootOf` en
  coefficient, signe de `CRootOf` décidé nativement, racines entières
  extraites premier par premier, contenu flottant extrait, `apart` avec
  flottants ; écarts : racine double de `nroots`, `minimal_polynomial` des
  sommes de `CRootOf`), solve-rest (`3ef3b67` : sondes 249 → 273/303 ;
  **résultat faux corrigé** : `solveset(x + Abs(x - 1) - 1, x, Reals)` donnait
  `{1}` au lieu de `Interval(-oo, 1)` ; `periodicity` porté,
  `minsolve_linear_system`, options `exclude`/`manual`, `Abs`/`Min`/`Max`/
  `Piecewise` dans `solveset` ; choix : `Abs(sin(x)) < 1/2` sur `[0, 2*pi]`
  résolu là où SymPy laisse un `ConditionSet` ; reste :
  `solveset(x**(1/3) + 2, x, Reals)` donne `EmptySet`, SymPy `{-8}`),
  ode-forms (`f677927` : corpus EDO 296/263/259 → 298/272/267 avec
  algfield ; série `hyper` 2F1 de SymPy fausse, la nôtre gardée),
  integ-rest (`822a0c3` : règles de `Integer._eval_power` et regroupement
  de `Mul.flatten` portés, `1922**(2/3)` → `31*124**(1/3)` ; puissance
  énorme gardée tenue ; `ratint(1/(x**3 + x + 1))` toujours en xfail).
  Validée le 23 octobre 2026 (8 shards, cargo test, clippy, corpus EDO
  298/272/267, corpus d'intégration 1109/71/18/10/2) ; main avancé à
  `f55878a`.
- [x] Vague 6, mergée : prec (`92c6171` : 6225 évaluations comparées à
  SymPy et à mpmath en haute précision, aucun cas où SymPy est juste et symox
  faux, 188 où seul symox est juste ; **résultats faux corrigés** :
  dépassements silencieux des flottants 53 bits (`Float(2)**1024` donnait
  `inf`, `Float(1e300)*Float(1e300)`...), partie imaginaire parasite de
  `Pow(-2, 10**20 + 1).evalf()` ; correction intentionnelle :
  `gamma(10**20/3)` faux de 52 ordres de grandeur dans SymPy ;
  `factorial(10**6).evalf()` 169 s → 3 ms ; reste : exposants binaires bornés
  à `2**100`, symox laisse `exp(10**30).evalf()` non évalué), simp-mat
  (`c0292d8`, Python : sondes `sympy/simplify/tests` 1544 → 1573 identiques,
  `sympy/matrices/tests` 1649 → 1728 et exceptions 145 → 73 ; `nsimplify`
  PSLQ et `rational=True`, `Expr.equals` par `minimal_polynomial`, étapes de
  `simplify` de SymPy autour de la recherche native, définitude, valeurs
  propres de matrices flottantes par mpmath), integ-risch (`c646cf0`,
  Python : sonde de 732 cas 620 → 634 identiques, non évalués 19 → 14 ;
  deltas de polynômes, produits avec `Heaviside`/`sign`/`Abs`,
  `f(x)*DiracDelta(x)`, `Eq`/`Sum`/bornes `None`/`O` ; corrections
  intentionnelles : échelle des dérivées de delta d'argument affine, racine
  réelle de `DiracDelta(x**5 - x - 1)` (SymPy 0)), ode-solve (`8103e48` :
  équations `Mod` sur les entiers, équations à radicaux, repli d'inversion de
  `solve`, `roots(trig=True)`, asserts `solveset` de SymPy 25 → 60/65 ;
  **résultat faux corrigé** dans `exptrigsimp` (une racine cubique perdue,
  bogue aussi dans SymPy) ; `solveset(x**(1/3) + 2, x, Reals)` reste vide,
  SymPy `{-8}` est faux pour la racine principale). Validée le 23 octobre
  2026 (8 shards, cargo test, clippy, corpus EDO 298/273/268, corpus
  d'intégration 1109/71/18/10/2) ; main avancé à `5ed3090`.
- [x] Vague 7, mergée : mat-rest (`1ba86db` : sonde des matrices 3494 →
  3625 identiques, exceptions 434 → 387, aucune régression ; une matrice
  symbolique singulière lève `NonInvertibleMatrixError` au lieu de rendre
  `zoo`/`nan` ; valeurs propres `CRootOf` au-delà du degré 4, Jordan
  flottant, `adjoint(sin(A))`, erreurs d'indexation, `solve` avec une matrice
  de symboles), simp-rest (`6a52595` : sonde `sympy/simplify/tests` 1636 →
  1769 identiques avec le banc corrigé ; port complet de `cse`, `sech`/`csch`
  dans `trigsimp`, `inverse=True` par `inversecombine`, `kroneckersimp`),
  bigexp (`170e033` : exposants binaires bornés à `2**120` au lieu de
  `2**100`, `exp(10**30)`, `Float(2)**(2**101)`, `factorial(10**30)`
  s'évaluent ; bogue de précision de `x**y` corrigé, parties réelles ou
  imaginaires minuscules exactes ; corrections intentionnelles vérifiées par
  mpmath : `(Float(2)**(2**40))**(1/3)`, `atan(I*(1 + 10**-20))` fini là où
  SymPy rend `oo`), solint-rest (`2553dce` : méthode de Lambert bivariée,
  `x*exp(x**2) - 1`, `x**a*exp(x) - 1`, `x**x - 2` ;
  `sqrt(x)*sqrt(x - 7) - 12` donne `{16}`, forme `Complement` ; **résultats
  faux corrigés** : `solveset((x - a)/(x - b))` donnait `{a}` sans condition,
  `limit(Ei(-x) - Ei(-2*x))` donnait `nan` ; racines non réelles de radicaux
  filtrées sur `Reals` ; sonde 864 → 872 ; restent Frullani,
  `log(log(x))/x**n`, `sqrt(log(1/x))`, `SingularityFunction`,
  `ratint(1/(x**3 + x + 1))`, le délai de `1/(1 + exp(x - 1/2)/(1 + exp(x)))`
  ; résultat faux préexistant trouvé : `solveset(Eq((-2)**x, 4), x, Reals)`
  vide, SymPy `{2}`). Validée le 23 octobre 2026 (8 shards, 43 273 tests,
  cargo test, clippy, corpus EDO 298/273/268, corpus d'intégration
  1109/71/18/10/2) ; main avancé à `a2a0038`.
- [x] Vague 8, mergée : bigexp2 (`801b8ba` : exposants binaires de taille
  quelconque, entier jusqu'à `2**120` puis grand entier, limite à 65536 bits ;
  `exp(10**100)`, `exp(10**1000)`, `Float(2)**(2**121)`, littéraux `1e(10**50)`
  ; `floor`/`int` d'un flottant géant lève `OverflowError` comme SymPy, avant
  `floor(Float(2)**(2**100))` donnait 1 ; précision de `ln` corrigée
  (`log(2).evalf(2500)`) ; corrections intentionnelles vérifiées par mpmath :
  `gamma(10**40/3)`, `exp(10**40).evalf()**(1/3)` ; restent
  `exp(2**70000).evalf()` non évalué et l'affichage lent près de la limite),
  simp-rest2 (`96b2919` : **résultat faux corrigé** : `(-2)**x = 4` sur
  `Reals` donne `{2}`, port de `_invert_real` pour les bases négatives ;
  `I**x = 1` reste un `ConditionSet` au lieu de `{0}` faux ; `nthroot`,
  `_sqrt_numeric_denest`, aides `minpoly` ; sonde simplify 1776 → 1782 ;
  corrections intentionnelles : `(-2)**x = 1/4` → `{-2}` (SymPy lève
  `ValueError`), `(-2)**x = 2**x` reste un `ConditionSet` (SymPy `EmptySet`
  faux) ; restent `nc_simplify`, `hyperexpand`, `trigsimp_groebner`,
  `powsimp_polar`, `posify`), mat-rest2 (`0276995` : sonde des matrices 3641 →
  4174 identiques, exceptions 389 → 245, aucune régression ; `M**S.Half`,
  puissances flottantes, `pinv(method='ED')` (blocage → 0,1 s), dérivée par
  une matrice, modules hérités `matrices.common`/`_legacy_*` ; restent
  `ImmutableMatrix` comme vraie classe, `Eq` avec une matrice, `Dagger` dans
  l'adjoint, formes non simplifiées), replay (`c80ea3f` : rejeu de l'audit
  des docstrings, 256 → 228 objets différents, c 30 → 13, f 15 → 5 ;
  régressions corrigées : `FunctionMatrix` (52 s → 0,5 s), `str` de
  `BinaryRelation`, `gaunt` (93 s), `expand` qui évaluait `1/(0.99 + 0.2*I)` ;
  nouvel écart de précision (e) : somme de `CRootOf` qui s'annule sous `abs`
  dans `eval_sum`). Référence EDO de SymPy régénérée après la coupure de
  courant (18 cas où SymPy dépasse 600 s), base recalculée sur `a2a0038` :
  301/275/270. Validée le 23 octobre 2026 (8 shards, 43 346 tests, cargo test,
  clippy, corpus EDO 301/275/270, corpus d'intégration 1109/71/18/10/2) ;
  main avancé à `0f528d3`.
- [x] Vague 9, mergée : solint-rest2 (`e99c78f` : **résultat faux corrigé** :
  la somme de `CRootOf` qui s'annule sous `abs` donnait `6.6e-24`, elle donne
  un zéro sans chiffre significatif comme SymPy ; développements de `Ei`,
  `Ci`, `Chi`, `expint` en 0, intégrales de Frullani `log(2)` ;
  `SingularityFunction` sur un intervalle infini (correction intentionnelle :
  SymPy rend 0 pour deux intégrales divergentes, symox `oo`) ;
  `ratint(1/(x**3 + x + 1))` (forme radicale plus grosse que SymPy) ; formes
  `asec`/`acsc`/`acot` de `solveset` ; restent
  `integrate(1/(1 + exp(x - 1/2)/(1 + exp(x))), x)`, les limites de
  `Ei`/`Ci`/`Si`/`erfc` à l'infini), simp-rest3 (`7a427e7` : sonde simplify
  1695 → 1710, corpus d'intégration SAME 1109 → 1134 ; `(I**x).is_finite`,
  `solveset(I**x*(x - 2), x, Reals)` → `{2}`, `(r**2)**x` → `Abs(r)**(2*x)`,
  budget de Hensel de `factor_list` (`nthroot1`), `nc_simplify`,
  `trigsimp_groebner`, `powsimp` polaire ; correction : `trigsimp` groebner
  sur des angles fractionnaires (SymPy lève `ValueError`) ; restent la forme
  de `simplify((x + 5)/(16*x**2 + 28*x + 6))`, issues 17141 et 11004, `log`
  à deux arguments tenu, ordre symbolique de dérivée), audit9 (`dbb4143` :
  5/16 objets c/f identiques à SymPy : groupes de Galois de degré 5 et 6,
  `julia_code`/`octave_code` à argument matriciel, `shape`,
  `aesara_function` ; `nsolve` suit `findroot` de mpmath ; corrections :
  dérivée en paramètre d'`assoc_laguerre` fausse dans SymPy ; restent
  `DiscreteMarkovChain` à états `Str`, `MultivariateNormal` à paramètres
  `MatrixSymbol` (`NotImplementedError`)), mat-rest3 (`dc910d3` :
  `MutableDenseMatrix`/`ImmutableDenseMatrix` vraies classes, `Eq`/`Ne` de
  matrices évalués comme SymPy, `_csrtodok` ; sonde 4182 → 4212, assertions
  3468 → 3544 ; restent `Dagger` dans l'adjoint, `Trace(X*X).diff(X)`,
  `DFT`/`GenericIdentity`, entrées `Poly`). Validée le 24 octobre 2026
  (8 shards, 43 402 tests, cargo test, clippy, corpus EDO 301/275/270, corpus
  d'intégration 1134/46/18/10/2) ; main avancé et poussé à `fa9f814`.
- [x] Vague 10, mergée : liminf (`8468eef`) : limites des fonctions
  spéciales en ±oo (Si, Ci, Ei, li, erf, expint, fresnel...), sonde 15 → 61
  sur 71, 9 des 10 restes sont des corrections intentionnelles (SymPy faux,
  vérifié avec mpmath) ; `log(log(x))/x**2` sur `(1, oo)` = `-EulerGamma`,
  `sqrt(log(1/x))` sur `(0, 1)` = `sqrt(pi)/2`, `limit(li(x) - li(2*x), x, 0)`
  = 0 ; l'intégrale risch `1/(1 + exp(x - 1/2)/(1 + exp(x)))` répond en
  0,02 s. Reste : trois limites en `AccumBounds` (`x*Ci(x)`,
  `x*(Si(x) - pi/2)`, `x*(fresnels(x) - 1/2)`). icorpus (`10d1cd9`) : résultat
  faux corrigé (`k**2 > 1/(4*pi**2)` avec `k` positif incluait
  `0 < k < 1/(2*pi)`, d'où `fourier_transform(sin(x)/x)`), forme de Laplace de
  `t**a`, point de coupure `x = 0` de Meijer G ; corpus d'intégration SAME
  1134 → 1161, EQUIV 46 → 21, COND 2 → 0. Reste : 17 MISSING volontaires,
  `cos(s*x)` sur `(0, oo)` (limite de `sin(s*x)`), 10 DIFF où SymPy est faux
  (`laplace_transform(1/(t+1))`), 21 formes. simp-rest4 (`15911a5`) : les
  flottants de haute précision gardent leur précision dans `cancel`,
  `factor`, `gcd`, `div`, `apart`, `Poly`... (issue 11004, domaine `RR` à
  précision variable comme `RealField`) ; `Add._eval_is_zero` porté (sonde de
  `is_zero` des puissances complexes : 105 → 13 écarts). Reste : règle
  d'égalité de taille du choix de `simplify` (`(x + 5)/(16*x**2 + 28*x + 6)`,
  issue 17141), `log` à deux arguments retenu, ordre de dérivation symbolique,
  `MatrixExpr.doit(inv_expand)`, `I**I`. mat-rest4 (`34bcf7b`) : dérivées par
  un `MatrixSymbol` portées (`Trace(X*X).diff(X)` = `2*X.T`,
  `test_derivatives.py` 13 → 67 sur 107, sonde matricielle 4212 → 4266),
  résultat faux corrigé (`(x*x.T*x).diff(x)`), `simplify(SparseMatrix)` rend
  une `ImmutableSparseMatrix`, `Eq(M, N).lhs` est une `ImmutableDenseMatrix`.
  Reste : `Dagger(a)` dans l'adjoint, noms de classe des expressions
  matricielles, `DFT`/`GenericIdentity`, entrées `Poly`, dérivées de Hadamard
  et par élément. Validée le 24 octobre 2026 (8 shards, 43 526 tests, cargo
  test, clippy, corpus EDO 301/275/270, corpus d'intégration
  1161/21/18/10/0) ; main avancé et poussé à `6072b72`.
- [x] Vague 11, mergée : mat-rest5 (`7177fc0`) : adjoint d'une matrice
  d'`Operator` en `Dagger(a)`, `adjoint`/`transpose` de `Abs(A)` et
  `sqrt(A)`, les expressions matricielles sont de vraies classes (`type(X*Y)`
  est `MatMul`), avertissement de SymPy pour les entrées `Poly`,
  `test_smith_normal` ; sonde matricielle 4266 → 5320, aucune régression.
  Reste : `DFT`, `IDFT`, `GenericIdentity`, `GenericZeroMatrix`,
  `Matrix == RawMatrix`, dérivées de Hadamard et par élément. limacc
  (`3a07044`) : sonde `test_limits` + `test_gruntz` 462 → 493 sur 648, 3
  segfaults de gruntz supprimés ; limites `AccumBounds` de `x*Ci(x)`,
  `x*(Si(x) - pi/2)` et Fresnel ; `sin(s*x)` en `oo` ; résultats faux corrigés
  (`cos(m*x)/x` complexe, `sin(x)**(1/x)`, `oo/x`, limites unilatérales sur
  les coupures). Reste en résultats faux : issues 22334, 6364, 22220, 6366,
  6843, 18399, `frac(x)**(1/x)`, `acoth(exp(1/x))` en `0-`. intrest3
  (`72f7dca`) : résultats faux corrigés (garde `Ne(n, 0)` des puissances
  trigonométriques, `integrate(sin(n*x)**2, (x, 0, pi))` donnait `pi/2` ;
  intégrales rationnelles impropres à pôles non placés ; pente de
  `DiracDelta(a*x - b)`) ; pôles en `±I*a`. Corrections intentionnelles :
  `1/(x**2 + a)` sur `(0, oo)` reste non calculée quand `a` peut être négatif
  (SymPy `pi*a/(2*polar_lift(a)**(3/2))` donne `-I*pi/2` en `a = -1` pour une
  intégrale divergente), `DiracDelta(a*x)` sur `(-1, 1)` vaut `1/Abs(a)`
  (SymPy 0). simp-rest5 (`7191c74`) : faits réels de `I**I` et des puissances
  imaginaires, `log(exp(x), 2, evaluate=False)` retenu puis réévalué comme
  SymPy, `MatrixExpr.doit(inv_expand=False)`. La règle d'égalité du choix de
  `simplify` a été mesurée (sonde +10, EDO 303/279/274) mais retirée : elle
  cassait 7 tests. Validée le 24 octobre 2026 (8 shards, 43 674 tests, cargo
  test, clippy, corpus EDO 301/275/270, corpus d'intégration
  1158/21/21/10/0, les 4 MISSING en plus étant les corrections
  intentionnelles ci-dessus) ; main avancé et poussé à `1d44c59`.
- [x] Vague 12, mergée : audit12 (`368b1cb`) : `power_base` de `expand`
  porté, lot c2 30 → 31/40 (c1 et c4 déjà complets). mat-rest6 (`43bdd4a`) :
  `DFT`, `IDFT`, `GenericIdentity`, `GenericZeroMatrix` natifs,
  `Matrix == RawMatrix`, `Derivative(X, A)` retenu ; sonde matricielle
  5305 → 5331. limwrong (`9384a36`) : 7 des 8 limites fausses corrigées
  (coefficients dominants nuls cachés prouvés) ; sonde `test_limits` +
  `test_gruntz` 493 → 516 sur 648 ; l'issue 22334 reste une `Limit`. simp-tie
  (`6267956`) : `simplify` rend la séquence de SymPy quand elle s'applique,
  `count_ops` mesure comme `fraction` (`(x + 5)/(2*(8*x**2 + 14*x + 3))`,
  issue 17141), sélection native gardée dans les solveurs de systèmes
  polynomiaux, sommes directes des pmf discrètes ; sonde simplify 1721 →
  1730, aucune régression. Reste : issue 22210 (ordre de dérivation
  symbolique), elle touche `solve.rs`. Validée le 24 octobre 2026 (8 shards,
  43 713 tests, cargo test, clippy, corpus EDO 303/279/274, corpus
  d'intégration 1159/20/21/10/0) ; main avancé et poussé à `d7db1d5`.
- [x] Vague 13, mergée : limrest (`a087132`) : Bessel en `±oo`, séries de
  `LambertW`, `asec`/`acsc`/`asech`/`acsch`, `asin(exp(x))` ; sonde 516 → 560
  sur 648 ; correction intentionnelle
  `limit(besseli(1, x)*sqrt(x)/exp(x), x, oo) = 1/sqrt(2*pi)` (SymPy
  `sqrt(2*pi)/2`). intrest4 (`1f0648c`) : résultat faux corrigé,
  l'antiprimitive Meijer G de `Heaviside(x - a)` n'est valable que pour
  `x > 0` (SymPy donne 2 pour `integrate(Heaviside(x - 1/2), (x, -1, 1))` au
  lieu de 1/2) ; `_simpc` de `laplace` porté. mat-rest7 (`36ee81e`) :
  `(A*B)[i, j]` en `Sum` à indice muet, `KroneckerDelta` à intervalle,
  dérivées par un élément de matrice (`Inverse`, `Sum`, `Trace`) ; sonde
  matricielle 5331 → 5348. Deux tests laissaient la précision de mpmath
  modifiée (`ac29ef8`). Restes : issue 22334, floor/ceiling, `Mod`,
  `PoleError` de gamma/digamma, dérivées de Hadamard et `applyfunc`,
  `Derivative(X, X)`, forme Meijer G de `cos(x)/x`. Validée le 24 octobre
  2026 (8 shards, 43 795 tests, cargo test, clippy, corpus EDO 303/279/274,
  corpus d'intégration 1161/18/21/10/0) ; main avancé et poussé à `0fe7228`.
- [x] Vague 14, mergée : limrest2 (`8155e64`) : limites de
  gamma/digamma/loggamma à argument à croissance logarithmique, `besselj(1, x)`
  en `±I*oo`, `besselk(1, x)` en `-oo`, `x*Shi(x)*I` ; sonde 560 → 570 sur
  648. `limit(x**2*Ci(x), x, oo)` reste non évaluée (la limite n'existe pas,
  SymPy rend `Ci(1/_w)/_w**2`). mat-rest8 (`ea37e76`) : dérivées de
  `HadamardProduct`, `HadamardPower` et `applyfunc` par le chemin des
  expressions de tableaux de SymPy ; `Derivative(X, X)` retenu (le `doit` de
  SymPy y rend 0, faux) ; sonde matricielle 5350 → 5366. audit14
  (`b34e1ea`) : tri des variables de `Derivative` comme
  `_sort_variable_count`, `Derivative(f(g(x)), g(x))`, `Subs.doit` et dérivée
  d'un `Subs` portés (3 résultats faux corrigés) ; correction intentionnelle
  `diff(f(g(x), g(x)**2), x)` (la forme de SymPy donne `5*sin**2*cos` au lieu
  de `3*sin**2*cos`). Restes : issue 22334, floor/ceiling/`Mod`,
  `powdenest(polar=True)`, forme Meijer G de `cos(x)/x`. Validée le 24
  octobre 2026 (8 shards, 43 850 tests, cargo test, clippy, corpus EDO
  303/279/274, corpus d'intégration 1161/18/21/10/0) ; main avancé et poussé
  à `74ce45c`.
- [x] Vague 15, mergée : pw-nc (`4271436`) : commutativité d'un `Piecewise`
  tirée de ses branches (`C*p*A` gardait un mauvais ordre), `evalf` partiel
  gardant `-1` et `I`, `latex`/`pretty` de `adjoint` ; matrices à entrées non
  commutatives : `PolynomialError` de SymPy pour `charpoly` et ses dérivées
  (des valeurs fausses étaient rendues), `echelon_form`/QR non pris en charge ;
  sonde 261 → 319 sur 352. Correction intentionnelle : `inv(method='GE')` et
  `'ADJ'` de SymPy sont faux en non commutatif (écart à l'identité 2,05 et 480
  sur blocs 2x2 aléatoires, LU 6e-16) ; symox passe par LU. matnum
  (`28c5f2b`) : `exp`/`log`/`inv`/vecteurs propres de matrices flottantes par
  la route de SymPy, précision > 53 bits gardée, pivot faux corrigé
  (`cos(2*x) - cos(x)**2 + sin(x)**2` pris pour non nul dans `nullspace`) ;
  sonde 178 → 222 sur 230. simp-rest6 (`bba1825`) : `powdenest(polar=True)`
  (`mellin_transform(besselj(0, x))` comme SymPy), `expand(c**-1*a**-1)` non
  commutatif corrigé (résultat faux), `nc_simplify`, `MatPow`/`MatMul`,
  `besselsimp`, `zoo` ; sonde simplify 1797 → 1822. audit15 (`5e33e81`) :
  matrices dans `subs`/`xreplace`, intersection `FiniteSet` et `ImageSet`,
  `beta`, impression de tableaux ; lot a 127 → 138 sur 199. Correctif
  d'intégration (`73756c8`) : les gardes non commutatifs s'installent après
  les routes flottantes. Validée le 24 octobre 2026 (suite sur demeter en 48
  shards, 44 003 tests plus 72 tests de compilateurs en local, cargo test,
  clippy, corpus EDO 303/279/274, corpus d'intégration 1171/8/21/10/0) ; main
  avancé et poussé à `73756c8`.
- [x] Vague 16, mergée : algroot (réduction exacte des expressions contenant des
  objets racines, radicaux par facteur premier), rrcc (domaines `RR`, `CC` et
  gaussiens composites), branchcut (précision sur les coupures et près des
  singularités). matparam mergé (`f39245b`) : `to_DM` et domaines comme
  SymPy (`ZZ[sqrt(a)]`, `RR[a]`, `ZZ_I[a]`, EXRAW), `inv` sur `RR[a]`/`CC[a]`,
  coefficients de `charpoly` simplifiés (`eigenvals` et `singular_values`
  de `sin(x)**2 + cos(x)**2` comme SymPy), `cancel(log(1.0*I))` en
  `I*pi/2`, `exp` symbolique par la forme de Jordan ; sonde 273 → 301 sur
  308. Validé (44 104 tests, corpus inchangés) ; main poussé à `8878eab`.
  rrcc mergé (`e19b01e`) : `invert`/`gcdex` sur `RR` gardent la précision
  au-delà de 53 bits (des rationnels étaient rendus), `trunc` sur `RR`,
  `Poly(0.0, x)` sur `RR`, `numer`/`denom` de `RR(x)`, unification de
  `RR[y]` avec `ZZ[y]`/`QQ(y)` (levait `TypeError`), `content` à plusieurs
  générateurs sur `RR[y]` (rendait `3/2` au lieu de `1.0*y`) ; sonde 220 →
  233 sur 269. Correction intentionnelle : `Poly((I*x + 1)**2, x,
  domain=ZZ_I).sqf_list()` de SymPy donne un produit de signe opposé ; symox
  garde ses facteurs avec l'unité `-1`. Reste : `domain=CC` explicite non
  pris en charge, `invert`/`gcd`/`gcdex` sur `CC` inféré, PGCD sur `RR[y]`
  multivarié. Validé (44 052 tests, corpus inchangés) ; main poussé à
  `ea48129`. algroot mergé (`8cd91d7`) : recombinaison de Zassenhaus
  (test de norme l1, budget par candidat), donc `factor(x**6 + 2,
  extension=(-2)**(1/6))` comme SymPy ; intégrales rationnelles définies
  évaluées quand Sturm prouve l'absence de pôle réel sur l'intervalle ;
  `evalf`, `N`, `latex` et `subs` de `RootSum` ; `AlgebraicNumber.to_root`,
  `is_aliased`, `to_primitive_element(radicals=False)` ; sonde 257 → 262 sur
  273. Correction intentionnelle : `RootSum(x**5 - x - 1, Lambda(y, log(x -
  y))).subs(x, 2)` vaut 0 chez SymPy, `log(29)` chez symox (vérifié
  numériquement). Reste : `factor(x**12 + 2)` sur extension,
  `sqrt(CRootOf**2)`. Validé (44 058 tests, corpus inchangés) ; main poussé
  à `cedd024`. branchcut mergé (`acfbb82`, `51c4be3`) :
  argument flottant sur une singularité comme SymPy (`erfinv(1.0)`,
  `elliptic_k(1.0)`, `gamma(-3.0)` en `zoo`, `Float(0.0)**-1` en `zoo`),
  réévaluation à plus haute précision près d'un pôle, bits ajoutés au noyau
  `gamma` près d'un pôle, côté de la coupure de `meijerg` sur `exp_polar`
  exact ; les interruptions (`KeyboardInterrupt`, signaux) arrêtent
  l'évaluation au lieu d'être réessayées. Corrections intentionnelles
  (vérifiées avec mpmath à 60 chiffres) : `N(zeta(1 + 10**-30))`,
  `N(elliptic_k(1 - 10**-30))`, `N(erfinv(-1 + 10**-30))`,
  `N(gamma(-3 + 10**-30), 30)`, `N(loggamma(-3 + 10**-30))` et
  `N(log(1 + 10**-40), 30)` justes là où SymPy lève, donne `oo` ou perd des
  chiffres. Reste : `N(hyper([1, 1], [2], 1.0))` non évalué (SymPy `oo`),
  `uppergamma(-n, 0.0)`, `exp(2**70000).evalf()`. Vague validée le 24
  octobre 2026 (44 138 tests, cargo test, clippy, corpus EDO 303/279/274,
  corpus d'intégration 1171/8/21/10/0) ; main poussé à `099dd0e`.
- [x] Vague 17, mergée (dernière vague de la session du 24 octobre 2026).
  ccdom (`88bdd53`) : `domain=CC` explicite et inféré (`ComplexField(100)`,
  `CC[y]`, `CC.frac_field`), PGCD multivarié sur `RR[y]`, `sqf_list`
  multivarié sur `ZZ_I[y]`/`QQ_I(y)`, branche demi-entière de
  `Pow._eval_power` (`sqrt(CRootOf(x**3 - 2, 1)**2)` donne `-CRootOf(...)`) ;
  résultat faux corrigé : `Poly(1.0*x + I, x).gcd(Poly(x**2 + 1, x))` rendait
  `1.0` ; sonde 267 → 374 sur 385. Reste : échelle du PGCD/PPCM sur `RR[y]`,
  `cancel(..., domain=RR)`, `Poly(x + 0.5*I)` sans générateur.
  reslimits (`78abdca`) : inventaire des budgets dans
  `docs/algorithmic-parity.md` ; résultats faux corrigés à l'épuisement d'un
  budget : Wang déclarait irréductible après 60 points malchanceux, test de
  nullité des coefficients de série au-delà de 400 nœuds, partie entière de
  `Mod(exp(300), 7)` à 160 bits, `eigenvals` à grands termes constants,
  `checkodesol` au-delà de l'ordre 256, `decompogen` silencieux ; plafonds de
  factorisation relevés vers ceux que SymPy atteint (`factor(x**12 + 2,
  extension=(-2)**(1/12))` comme SymPy). Reste : `count_roots(x**5000 - x -
  1)`, `factor_list(x**300 + x + 1)`, grands `resultant`,
  `cyclotomic_poly(20001)`, `bernoulli(5002)` et grands nombres
  combinatoires non évalués là où SymPy calcule.
  defparam (`4335809`) : les 21 MISSING et 10 DIFF du corpus d'intégration
  sont des valeurs de SymPy fausses pour certains paramètres (vérifié par
  quadrature), conservées non évaluées ; distributions : primitive de
  `DiracDelta(a*x + b)`, deltas polynomiaux à coefficients symboliques,
  dérivées de delta fois fonction indéfinie, `Heaviside` d'argument pair par
  Meijer G, terme d'ordre des séries de `Heaviside(x - 1)` ; sonde 128 → 139
  sur 196. Reste : produits de distributions (`DiracDelta(x)**2`), delta en
  borne d'intervalle (demi-masse, SymPy 0).
  Validée le 24 octobre 2026 (44 338 tests, cargo test, clippy, corpus EDO
  303/279/274, corpus d'intégration 1171/8/21/10/0) ; main poussé à
  `f425cf0`.
- [x] Restes de solve (`6a4d3f9`, `a461a7d`, mergés ; Rust) : `sec`/`csc`
  en logarithmes, `exp(exp(x)) - 1` sur `ProductSet(Integers, Integers)`,
  coefficients symboliques (`sin(a*x) - b`...) ; 8 → 27/40 sondes. Le
  blocage de `nonlinsolve([x**3 + 1/y - 2, x - y])` est corrigé (0,2 s, les
  4 solutions de SymPy) ; sondes `nonlinsolve` 22 → 26/46 identiques, 0
  blocage. Restes : `sec(2*x) - 3`, systèmes à coefficients symboliques.
- [x] Test `test_limit_seq_sum_over_exponential` non déterministe
  (`da90fd7`, mergé ; Rust) : la règle du produit de `difference_delta`
  prenait le premier facteur stocké au lieu de l'ordre canonique, et le
  hachage des Dummy dépendait de l'horloge. 50/50 exécutions donnent 2.
  Reste : `limit_seq` imbriqué (SymPy 4, symox `None`), étape `dominant` non
  portée.

## 1. Corps algébriques généraux et racines complexes certifiées

- [x] Étendre les corps de coefficients au-delà des constructions actuellement
  prises en charge : racines quatrièmes et tours algébriques. Preuves :
  `factor(x**4-2*y**4, extension=2**(1/4))`, `gcd(x**4-2, x-2**(1/4),
  extension=True)` et `minpoly(sqrt(1+2**(1/4)))` (degré 8) identiques à SymPy.
- [x] Généraliser la construction de corps, les éléments primitifs, les plongements
  choisis et les tests exacts d'appartenance. Preuves :
  `QQ.alg_field_from_poly(x**5-x-1)`, `cyclotomic_field(7)`, `to_number_field`,
  `field_isomorphism`.
- [ ] Propager ces domaines dans les opérations polynomiales. Fait : PGCD
  multivariés, `sqf`, `factor`, `gcdex`, `invert`, `cancel`, `groebner`.
  Vague 2 : `resultant`, `discriminant`, `gcd`, `lcm` et `cancel` avec
  `extension=True` et des coefficients paramétrés (le `gcd` renvoyait `1` à
  tort). Vague 3 (`c9abbe9`) : `parallel_poly_from_expr` infère un domaine
  commun quand les domaines séparés diffèrent, comme SymPy.
- [ ] Compléter les anneaux de polynômes et corps de fractions. Fait : `Poly`
  sur `ZZ_I`, `QQ_I`, `EX` et `ZZ[sin(a)]`, unification `QQ(a)` avec `ZZ[pi]`.
  Vague 2 : anneaux imbriqués `ZZ[x][y]` et fabriques `ZZ_I[x]`,
  `QQ<sqrt(2)>[x]` (`tests/test_nested_rings_parity.py`) ; `Poly` sur
  `ZZ_I[y]`, `QQ_I[y]`, `ZZ_I(...)`, `QQ_I(...)` avec l'inférence de SymPy
  (`ZZ_I(N,pi)` pour `-2*I*pi*n/N`), `is_unit`/`canonical_unit` comme SymPy,
  141/148 cas identiques (`tests/test_gaussian_domains_parity.py`). Reste :
  éléments des corps de fractions sur `ZZ_I`/`QQ_I` (`construct_domain(...,
  field=True)` donne `EX`), `CC[y]` non inféré avec flottants,
  `quo_ground`/`exquo_ground` inexacts, normalisation du PGCD sur `QQ_I`,
  `sqf_list` multivarié, affichage de `.rep`, représentation plus lente que
  SymPy ; noms de classes et types d'erreurs différents. Doctests Sage à
  relancer. Vague 3 (`tests/test_gaussian_fractions_parity.py`, sonde 49 puis
  143 sur 151) : éléments des corps `ZZ_I(y)` et `QQ_I(y)` comme SymPy
  (`construct_domain(..., field=True)` ne tombe plus sur `EX`). Reste : `sqf_list`
  multivarié sur les domaines gaussiens composites (SymPy y est lui-même faux
  d'un signe dans plusieurs cas), petites différences de conversion.
- [x] Ajouter une isolation certifiée des racines complexes (rectangles
  certifiés, `CRootOf` avec ordre complexe exact).
- [ ] Étendre la réduction algébrique exacte des expressions contenant des
  objets racines et leur conversion en corps algébriques. Non commencé.
- [ ] Différences de domaines numériques et d'expressions. Vague 2
  (`tests/test_domains_misc_parity.py`) : conversions flottantes, `mpf` et
  complexes par `RR`/`CC` comme SymPy (`QQ.convert(2.5)` donne `5/2`), erreurs
  de refus identiques sur 400 cas, `EXRAW` refuse les objets non `Expr`,
  options de `construct_domain` (`extension=sqrt(2)` donne `ZZ`), affichage des
  éléments composites (`pi` et non `ZZ[pi](pi)`), `RR.invert` corrigé.
  Domaines composites inexacts à un générateur (vague 2,
  `tests/test_inexact_composite_parity.py`) : `RR[y]`, `CC[y]` et `RR(y)`
  inférés par `Poly` avec l'arrondi de `RealField(prec)`, division, PGCD,
  `cancel`, sqf, résultants, `gcdex` ; 61 lignes sondées sur 69 identiques.
  Vague 3 (`04d0855`, `c9abbe9`) : plusieurs générateurs (`Poly(1.5*x*z + y,
  x, z)` sur `RR[y]`), `construct_domain([1.5*x])` et `composite=True`
  (`RR[x]`, `RR[x,y]`), `ground_roots` sur `RR` par `factor_list` et
  `TypeError` de SymPy sur `RR[y]`, affichage des éléments (`2`, éléments de
  `CC[x]`), `CC.invert` lève le `TypeError` de SymPy. `cancel` sur `RR`
  garde les coefficients flottants (`rr_cancel`), `Poly.cancel` Python sur
  `RR` reste à aligner. Audit des suites `RR`/`CC` (rien corrigé) :
  produits et PGCD mêlant `RR[y]` et `ZZ[y]` lèvent `TypeError`,
  `factor_list`/`sqf_list`/`content` à plusieurs générateurs sur `RR`
  différents (facteurs affichés `1.0*x - 1.0*y`), `apart`/`invert`/`gcdex`
  sur `RR` rendent des rationnels, `domain=CC` explicite lève
  `NotImplementedError`, `RR.frac_field(x).field` et `RR.poly_ring(x).ring`
  lèvent `AttributeError`.
- [x] `(-2)**(1/6)` avec `extension=True` donne le corps `QQ<(-2)**(1/6)>` de
  SymPy (vague 3, `tests/test_extension_root_parity.py`, sonde 36 puis 79 sur
  86) : le noyau garde une base négative à exposant rationnel en une seule
  puissance, replie les puissances de `-1` dans les produits et distribue les
  puissances de produits de radicaux positifs comme `Mul.flatten`. Reste :
  extraction des racines par facteur premier (`(-12)**(1/6)` en
  `(-3)**(1/6)*2**(1/3)`), retirée car elle bloquait `simplify` dans le PGCD
  gaussien multivarié (borne à ajouter d'abord) ; `factor(x**6 + 2,
  extension=(-2)**(1/6))` atteint la limite de travail ; `nsimplify((-2)**(1/6))`
  non développé.
- [ ] Domaine `EX` : formes des racines et chiffres de `nroots` parfois
  différents de SymPy. Vague 3 (sonde 45 puis 47 sur 56) : corrections
  intentionnelles où symox est plus précis (racine triple de
  `(x - 1/3)**3` là où SymPy lève `NoConvergence`, petite racine `1e-20` de
  `x**2 - 1e20*x + 1`, racines doubles exactes). Reste : racines groupées de
  `(x - 1)*(x - 1.0001)*(x - 1.0002)` moins précises (2e-8 contre 1e-12),
  `roots` sur `EX` ne factorise pas d'abord (`x**2 - (1 + sqrt(2))*x +
  sqrt(2)` donne des radicaux imbriqués au lieu de `{1, sqrt(2)}`).

Validation : reproduire le test du corps contenant `2**(1/4)`, puis comparer les
domaines, plongements, identités exactes, multiplicités et raffinements de racines
sur des familles variées, y compris des racines complexes très proches.

Voir [polynômes](docs/polynomials.md), [domaines](docs/domains.md) et
[racines algébriques](docs/algebraic-roots-quaternions.md).

## 2. Intégration symbolique et transformées

- [x] **Résultat faux, fermé.** Les intégrales définies à borne infinie avec un
  paramètre sans hypothèse de signe produisent désormais les conditions de
  convergence de la référence (étapes Meijer G portées) ou restent non
  calculées.
- [x] Intégration définie Meijer G fidèle : étapes, conditions de convergence,
  `hyperexpand`, table, `meijerint_indefinite`, intégration des `Piecewise`
  dans `Integral.doit`.
- [ ] Compléter la couverture des stratégies Risch et heurisch. Fait : cas
  particuliers de `heurisch_wrapper`, forme Risch de poly*exp. Les branches de
  cas particuliers de heurisch diffèrent encore en forme (égales à une
  constante près).
- [x] Options de stratégie `risch`, `heurisch` et `manual` et
  `NonElementaryIntegral` (vague 2, `tests/test_integrate_methods_corpus.py`,
  1072 cas). Corpus de 294 exemples x 5 modes : 1339 identiques, 90 égaux par
  dérivation, 27 formes différentes, aucun résultat faux, 10
  `NonElementaryIntegral` prouvés des deux côtés. Reste : 9 dépassements de
  temps (`atan(x)*log(x)`, `heurisch=True` sur `x**(4/3)*log(x)`,
  `test_issue_10847_slow`), 3 intégrales non calculées où SymPy répond (Piecewise
  flottant de `test_issue_15498`), formes natives au lieu de `RootSum`
  (`atan(exp(x))`, `exp(x)/(1+exp(2*x))`, `1/(x*(1+log(x)**2))`, `risch=True`
  rend des fractions partielles natives), terme de correction `floor` absent
  pour `1/(cos(x) + 2)`. Vague 2 (`tests/test_ratint_parity.py`) : `RootSum`
  porté (construction, affichage, `diff`, `doit`), `ratint(1/(x**3+x+a))` donne
  le `RootSum` exact de SymPy, `manual=True` sur `atan(x)*log(x)` répond comme
  SymPy (2 fois plus lent). Reste : `ratint(1/(x**3 + x + 1))` ne lève plus mais
  sa forme fait 80 Ko contre 17 Ko et il est 5 fois plus lent (réduction des
  radicaux numériques : `1922**(2/3)` au lieu de `31*124**(1/3)` ; le portage
  essayé bloque la division `EX` de `log_to_atan`, correctif non appliqué).
  `checkpdesol` : `Integral.doit` et `Derivative.doit` refont comme SymPy, le
  `simplify` d'une EDP substituée donne 0 en 5 ms ; SymPy répartit encore
  l'intégrale retenue sur les termes de son intégrande, symox la garde entière.
  Vague 3 (`tests/test_integrate_rest_parity.py`) : terme de correction `floor`
  porté (corrige un `0` faux pour `1/(cos(x) + 2)` sur `(0, 2*pi)`), formes
  `RootSum` de SymPy, transformée de Fourier gaussienne sous condition,
  `test_issue_15498` résolu (`solve_lin_sys` sur `RR` comme SymPy, `cancel`
  sur `RR` garde les coefficients flottants) ; corpus 1374 identiques / 72
  égaux / 20 formes / 3 non calculés / 0 faux / 0 dépassement de temps.
  Correction intentionnelle : le `Piecewise` de correction `floor` est ajouté
  une fois (SymPy l'ajoute deux fois). Reste : forme de 80 Ko de
  `ratint(1/(x**3 + x + 1))`.
- [ ] Étendre les intégrales définies avec paramètres. 54 cas du corpus restent
  non calculés là où la référence répond : intégrales oscillantes sur
  `(0, oo)` rendues en `AccumBounds`, valeurs divergentes (`1/sin(x)` sur
  `(0, pi)`, `x**a` sur `(1, oo)`), `1/cosh(x)` sur la droite réelle,
  intervalles finis avec paramètre (`1/(x+a)`, `1/x` sur `(a, 1)`, `log(a*x)`,
  `Heaviside`, `DiracDelta`). Certaines valeurs de référence (`-log(a)` pour
  `1/x` sur `(a, 1)`) ne valent que si l'intervalle évite le pôle : ne pas les
  recopier sans condition. Vague 3 (`tests/test_definite_params_parity.py`) :
  42 puis 31 cas non calculés, aucun nouveau résultat faux ; 17 valeurs de
  SymPy volontairement non reprises (fausses pour certains paramètres), la
  valeur `-1 - 1/a` de SymPy corrigée en `oo`. Décision du 22 octobre : les 14
  intégrales oscillantes que SymPy rend en `AccumBounds` doivent être rendues
  comme SymPy (vague 4, branche accumbounds).
- [x] Compléter les transformées. Fait : Laplace, Fourier, Hankel et Mellin
  passent par l'intégrale, inverses de Mellin et de Laplace portées. Vague 2 :
  sinus, cosinus, Fourier, Hankel et leurs inverses passent d'abord par
  l'intégrale, tables en repli ; 94/111 cas identiques à SymPy (conditions,
  forme du tuple, `noconds`, erreurs), `tests/test_transforms_remaining_parity.py`.
  Reste : 15 cas de même valeur sous une autre forme (`simplify`), condition
  `True` au lieu de celle de SymPy pour `exp(-a*x)*sin(b*x)*Heaviside(x)`,
  `hankel_transform(1, r, k, 0)` lève `IntegralTransformError` là où SymPy
  renvoie un `Piecewise` de `nan`/`oo`, `Abs(w**2) < 1` non réécrit en
  intervalle réel. Détail dans `docs/algorithmic-parity.md`. Vague 3
  (`tests/test_transforms_wave3_parity.py`) : conditions réécrites par
  `piecewise_simplify` comme SymPy (`Abs(w**2) < 1` fermé). Reste : les 15
  formes, `hankel_transform(1, r, k, 0)`, la condition de Fourier de
  `exp(-a*x)*sin(b*x)*Heaviside(x)` (sans hypothèse, symox répond là où SymPy
  laisse la transformée non calculée).
- [x] `Piecewise` : inégalités périodiques, morceaux conservés, formes
  ITE/Min/Max (vague 2, `tests/test_piecewise_periodic_parity.py`, sonde
  42/90 -> 79/90). Ports de `piecewise_simplify`, `piecewise_exclusive`,
  `simplify_univariate`, puissance d'un `Piecewise` distribuée par branche.
  Reste : `sin(x) + cos(x) > 0` et `exp(sin(x)) > 1` restent en
  `ConditionSet`, `series` d'une somme ou d'un produit contenant un
  `Piecewise` lève `NotImplementedError`, terme d'ordre par branche
  (`O(x**3)` contre `O(x**5)`), simplification And/Or par motifs non portée
  (`Or(x < 0, x < 1)`, substitution `Eq`), `Heaviside(x, H0)` en `Piecewise`
  imbriqué (égal, forme différente).
  **Corrections intentionnelles validées le 22 octobre 2026** (SymPy 1.14 est faux) :
  `solveset(sin(x) > 0, x, Interval(2*pi, 3*pi))` donne `(2pi, 3pi)` (SymPy
  `EmptySet`) ; `sin(x) >= 0` sur `[0, 2pi]` inclut `2pi` ; l'intégrale de
  `Piecewise((1, sin(x) > 0), (0, True))` sur `(0, 4pi)` vaut `2pi` (SymPy
  `pi`) ; `simplify` garde `sin(x) > 0` au lieu de `0 < x < pi` sur la droite
  réelle ; `simplify_univariate(Or(x <= 0, x >= 1))` garde `x = 0`.
  Vague 3 (`tests/test_piecewise_series_parity.py`, sonde 48 puis 97 sur 112) :
  `series` des sommes et produits contenant un `Piecewise` (chemin
  `Expr.series`), terme d'ordre par branche, simplification And/Or par motifs
  de `boolalg`, `Order.subs` non affine, limites unilatérales d'un `Piecewise` ;
  `Heaviside(x, H0)` est un noeud natif à deux arguments
  (`tests/test_core_misc_wave3_parity.py`). Correction intentionnelle :
  `limit(Piecewise((sin(x), x < 0), (cos(x), True)) + x, x, 0, '-')` donne `0`
  (SymPy `1`, branche de droite). Reste : `sin(x) + cos(x) > 0` et
  `exp(sin(x)) > 1` en `ConditionSet`, `O(x*y).subs(y, x)` ne lève pas,
  `series` d'une fonction d'un `Piecewise` non décidé, motifs `Xor` non
  branchés.
- [x] `meijerint_indefinite` propage les erreurs que SymPy laisse sortir
  (vague 2, 60 intégrandes identiques en résultat, `None` ou erreur).
- [ ] Étendre l'intégration des distributions au-delà des deltas affines simples
  et des cas polynomiaux certifiés : paramètres, ordres dérivés et supports plus
  généraux, en conservant les conditions de validité.
  **Résultats faux, fermés (vague 3, `tests/test_heaviside_limit_parity.py`)** :
  `limit(Heaviside(x), x, 0)` donnait `1/2` dans toutes les directions (`1`
  par la droite, `0` par la gauche, `ValueError` pour `'+-'`), erreur héritée
  par toutes les limites passant par le saut ; `limit(DiracDelta(x), x, 0)`
  donne `0` ; limites, séries et termes dominants de `floor`, `ceiling`,
  `frac` aux sauts ; sonde 70 puis 194 sur 254. Une intégrale définie dont un
  delta de la primitive peut tomber sur une borne finie reste non calculée.
  Corrections intentionnelles : `Heaviside(sin(x))` en 0 (SymPy `1/2` des deux
  côtés), `Piecewise` par la gauche, `nseries` unilatéral,
  `integrate(DiracDelta(x - 1, 1), (x, 0, oo))` donne `0` (SymPy
  `DiracDelta(oo)`). Reste : `series(Heaviside(x)/x)` garde un `O(x**3)`,
  `DiracDelta(x).evalf(subs={x: 0})`.
- [x] Lcapy : 317 réussis, 0 échec (vague 2, `tests/test_lcapy_wave2_parity.py`).
  `cancel` fait le prétraitement de SymPy (`signsimp`, `factor_terms`), le
  `simplify` natif garde `factor_terms(sign=False)` et applique `exptrigsimp`
  comme SymPy, les nombres d'un `Mul` non évalué sont des coefficients. Reste :
  aucun (`test_DFT` construit maintenant `ZZ_I(N,pi)` comme SymPy).

Validation : vérifier les primitives par différentiation exacte là où elle
s'applique, et les intégrales définies avec leurs conditions de convergence.
Une continuation analytique ne doit pas être présentée comme la valeur d'une
intégrale divergente. Les conventions de branches et de poids aux bornes
doivent être explicites.

Voir [contrat d'intégration](docs/algorithmic-parity.md),
[compatibilité des transformées](docs/compatibility.md) et
[revue d'intégration](audits/sympy-1.14.0/integration-review.json).

## 3. Résolution, paramètres et complétude des solutions

- [ ] Étendre la résolution transcendante aux mélanges de puissances,
  exponentielles et logarithmes au-delà des familles Lambert W actuelles.
- [ ] Traiter les différents cas selon les paramètres, notamment les coefficients
  pouvant s'annuler et les composantes singulières des systèmes.
- [ ] Étendre les familles de solutions représentées par ensembles : branches
  complexes, paramètres entiers, restrictions réelles et intersections.
- [ ] Auditer toutes les options de `solve`, `solveset` et `roots` ; les options
  radicales déjà couvertes ne suffisent pas à établir la parité générale.
- [x] Forme des solutions de `solve` pour les polynômes de degré 3 et plus et
  pour les produits (vague 2) : `roots` composite comme `_solve`, `CRootOf` via
  `all_roots`, inversion par le générateur, chaque facteur d'un produit avec sa
  propre décision de simplification. Cubiques, quartiques, `x**6+x**3+1`,
  quintiques et septiques `CRootOf` identiques à SymPy
  (`tests/test_solve_polynomial_forms.py`). Audit : 21 différences de forme et
  5 différences d'ensemble sur ~120 cas, toutes racines vérifiées à 30 chiffres.
- [x] Solutions manquantes de `solve` (vague 2) : `x*log(x)` donne `[0, 1]`,
  `(x-1)*log(x-1)` `[1, 2]`, bases puissances parfaites (`4**x - 8`,
  `(1/4)**x - 8`, `8**x - 4`) avec leurs deux branches, `x**5 + a*x + 1` donne
  `[]` avec `incomplete=True` par défaut (lève avec `incomplete=False`),
  `tests/test_solve_completeness_parity.py`.
- [ ] Formes restantes de `solve`. Vague 2 (`tests/test_solve_sets_parity.py`) :
  repli `_tsolve` des bases puissances, rationalisation des flottants,
  hyperboliques par `exp`, exposants négatifs, bases mélangées, Lambert avec
  `log`, `f(x)**g(x) == c`, `nsimplify` des résultats Lambert ; sonde de 129
  équations : 89 puis 110 identiques, 0 erreur (9 avant). Reste : branche
  `-I*pi` de `2**(3*x)-3*2**x+1` (ordre), `exp(x)-x-2` garde la branche
  Lambert -1 que SymPy écarte, `2**x-4*x` et `x*2**x+1/4` lèvent
  `NotImplementedError` (SymPy résout), `exp(a-x)+1` donne `a-I*pi` contre
  `log(-exp(a))`, l'étape `shorter(expand_log, logcombine)` de `simplify` manque
  (`x**(2*x)-4` donne `log(4)/2` contre `log(2)`), ordre des clés des
  dictionnaires de solutions. Vague 3 (`tests/test_solve_forms_parity.py`,
  sonde 46 puis 65 sur 74, systèmes 21 sur 21, `rational=False` 11 sur 11) :
  `2**x - 4*x` et `x*2**x + 1/4` résolus, étape `expand_log`/`logcombine` de
  `simplify`, `exp(a - x) + 1` donne `log(-exp(a))`. Corrections
  intentionnelles : branches `W(-1)` supplémentaires conservées, `2**x - 4*x`
  ajoute `4`, `3**x - 3*x` ajoute `1`, `2**x - x - 2` ajoute `2` (solutions
  exactes que SymPy omet). Reste : formes seulement (base puissance, ordre des
  racines cubiques).
- [x] **Résultat faux, fermé (22 octobre 2026, `efa251c`).**
  `solve([1.5*x*y - 2, x + y - 1], [x])` donne `[]` comme SymPy : un reste
  d'élimination à numérateur développé non nul rend un système polynomial
  incohérent (paramètres génériques), dans les chemins linéaire et non
  linéaire ; les équations sans inconnue sont écartées comme dans SymPy
  (`tests/test_solve_overdetermined_parity.py`, 51 tests ; sonde de 158
  appels : 66 différences avant, 7 après, dont l'ordre des clés). Vague 3 :
  prétraitement `_invert` de SymPy porté, flottants gardés avec
  `rational=False`.
  `solve([x**2 - a, x - sqrt(a)], [x])` donne `[(sqrt(a),)]` (SymPy `[]`,
  alors que `x = sqrt(a)` convient) : résultat juste conservé.
- [ ] Développer la vérification exacte des candidats avec gestion des branches,
  des pôles et des identités algébriques. Un résidu indécidable conservé protège
  une solution possible, mais n'en prouve pas la validité.
- [ ] Étendre les domaines des inégalités et de l'analyse réelle lorsque la
  référence sait décider : signes paramétriques, continuité, points stationnaires
  et expressions non lisses.
- [ ] Relationnels et ensembles : inégalités périodiques et
  `as_expr_set_pairs` faits (vague 2). Reste :
  `solve_univariate_inequality(x**2 + 2 > 0, x, relational=False)` affiche
  `Interval(-oo, oo)` au lieu de `Reals` (égaux et de même hash depuis la
  vague 3, `9310a60` : différence d'affichage seulement).
- [x] Vacuité et inclusion des ensembles (vague 2,
  `tests/test_union_sets_parity.py`) : `is_subset` porté (cas triviaux,
  `_eval_is_subset`, règles `is_subset_sets`, repli par intersection),
  `Union.is_empty` et `ProductSet.is_empty`, `is_finite_set`, repli de la
  monotonie sur l'algorithme de SymPy (`is_increasing(sin(x))` donne `False`) ;
  sonde : 29 puis 49 cas identiques sur 50. Vague 3 (`9310a60`) : `Reals` égal
  à `Interval(-oo, oo)` avec le même hash, simplification native des unions
  sur les règles `union_sets` de SymPy, `Range.size` comme SymPy. Reste :
  `Range` symbolique dans `FiniteSet`.
- [ ] Ensembles, vague 3 (`tests/test_sets_rest_parity.py`,
  `tests/test_sets_intersection_parity.py`, `tests/test_sets_contains_parity.py`) :
  gestionnaires `imageset` et des fonctions d'ensembles, réseaux, `ProductSet`,
  ordre des unions (sonde 214 puis 281 sur 285) ; intersections `Range` avec
  `Range` et `Interval`, `ImageSet` par diophantienne, avec `Reals` et
  `Complexes` (147 puis 211 sur 218) ; inversion de `solveset` sans hypothèse
  sur le symbole, diophantienne parabolique, appartenance à `Range`,
  `ImageSet` et `ProductSet` (30 puis 96 sur 101). Corrections
  intentionnelles : `imageset(tan, Interval(0, 1))` donne
  `Interval(0, tan(1))` ; l'intersection d'un `ImageSet` complexe avec un
  intervalle passe par sa partie réelle (`ImageSet(n*I + 1)` sur `[0, 1]` donne
  `{1}`, SymPy `EmptySet`) ; `solveset(exp(I*x) - 1, x, Reals)` donne
  `ImageSet(2*n*pi)` et `exp(I*x) - I` sur `[0, 10]` donne `{pi/2, 5*pi/2}`.
  Reste : retrait de l'absorption `is_subset` du noyau Rust (bloqué), quelques
  formes ; détail dans `docs/algorithmic-parity.md`.
  **Correction intentionnelle validée le 22 octobre 2026** : `is_increasing(a*x, symbol=x)` avec
  `a` réel reste `None` ; SymPy 1.14 donne `False`, faux dès que `a >= 0`.

Validation : comparer les ensembles complets de solutions, les exclusions et
les multiplicités pertinentes. Vérifier seulement les solutions retournées ne
détecte pas les solutions manquantes. Conserver une représentation conditionnelle
ou une erreur explicite lorsqu'un calcul reste indécidable ou non pris en charge.

Voir [contrat des solveurs](docs/algorithmic-parity.md),
[revue des solveurs](audits/sympy-1.14.0/solver-review.json),
[inégalités](docs/inequalities.md) et [analyse réelle](docs/real-analysis.md).

## 4. Équations différentielles

- [x] Ajouter les familles d'EDO résolues par la référence : modules d'EDO de
  SymPy 1.14 portés, tous les indices de `allhints` disponibles, solveurs natifs
  utilisés seulement quand toutes les méthodes de référence échouent.
- [ ] Compléter les méthodes et options. Fait : `hint='all'`, `'all_Integral'`,
  `'best'`, `simplify=False`, `ics`, `xi`, `eta`, `n`, `x0`. Reste : formes de
  sortie quand `simplify`, `integrate` ou `factor` choisissent un autre candidat
  équivalent (corpus : 245/315 `dsolve` par défaut identiques), indice
  `factorable` dépendant de la forme (2 classifications), `rational_riccati` et
  cas dépendant de l'intégration. SymPy applique `factor_terms` sans condition ;
  Symox ne le garde que s'il ne fait pas grossir l'expression.
- [x] Séries de Frobenius (vague 2) : sortie par défaut de
  `2nd_power_series_regular` identique à SymPy sur les 37 équations auditées
  quand SymPy répond (17 identiques au caractère près, 5 mêmes erreurs),
  `tests/test_frobenius_parity.py`. Corrections communes : `cse` n'entre plus
  dans les `Order`, `terms_gcd` ne fusionne plus `x` et `x - 1`, `subs` découpe
  les produits par `fraction`. Extensions documentées : deuxième solution
  logarithmique (racines répétées ou écart entier, `sol.subs(C2, 0)` redonne
  SymPy), `p`/`q` non polynomiaux développés correctement, racines
  irrationnelles là où SymPy lève. Reste (comme SymPy) : racines indicielles
  complexes et ordre symbolique (`ValueError`), écart de racines `>= n` donne
  `O(1)`, `O(x**n)` en 0 quand `x0 != 0`, solution logarithmique non normalisée.
- [x] Limites de 64 coefficients et d'ordre ou de taille : ne s'appliquent plus
  qu'au repli natif. Les ordres au-delà de `10**6` lèvent `NotImplementedError`
  dans `classify_ode`, là où SymPy boucle ou épuise la mémoire.
- [x] Systèmes différentiels : `dsolve_system` de SymPy porté (types 1 à 6,
  ordre supérieur, décomposition en composantes, formes non linéaires à 2 et
  3 équations), méthodes natives en repli.
- [x] EDP du premier ordre : les 4 indices de SymPy, repli par
  caractéristiques conservé. EDP non linéaires et méthode de Charpit absentes,
  comme dans SymPy.
- [ ] **Performance** : la simplification native est beaucoup plus lente que
  SymPy sur les grandes expressions de l'heuristique de groupes de Lie
  `abaco2_unique_general` (build debug : 4 s / 56 s / plus de 900 s contre
  2,3 s / 9,6 s / 67 s). Le `dsolve` par défaut de
  `(x + cos(f))*f' + f + 2*x` dépasse 40 minutes en build release (le repli natif
  ne vient qu'après), contre 82 s pour SymPy avant son `NotImplementedError` et
  une réponse immédiate sur main (`Eq(x**2 + x*f(x) + sin(f(x)), C1)`). Le temps
  est passé dans `simplify` appelé par `lie_heuristic_abaco2_unique_general`
  (chemin `factorable`, puis `dsolve` interne, puis `lie_group`). Vague 2
  (`feat-simplify-perf`) : 395 s en release contre plus de 40 min, sorties
  identiques ; reste plus lent que SymPy (82 s) et que main (réponse
  immédiate). Main avancé le 22 octobre 2026. Vague 3 (`cf8b9a4`) : la suite
  de réécritures de SymPy est un candidat de plus du `simplify` natif ; les
  deux premiers appels rendent les formes de SymPy (142 et 150 opérations,
  0,3 s et 1,3 s contre 2,3 s et 10 s), le troisième reçoit l'entrée de SymPy
  et prend 208 s contre 71 s ; `dsolve` 210 s contre 85 s. Corpus EDO
  296/262/258 (avant 296/259/257). Reste : 159 s dans les candidats
  `trigsimp` natifs (sommes d'environ 12 000 termes reconstruits par
  `args()`, mémo à comparaison complète), `together` natif sans `deep`.
  Vague 3 (`e5dc7f7`, moteur de règles trigonométriques interné et indexé) :
  troisième appel de 230 s à 77 s (SymPy 71 s), sorties identiques ; le
  `dsolve` par défaut de `(x + cos(f))*f' + f + 2*x` rend
  `Eq(x**2 + x*f(x) + sin(f(x)), C1)` en 60 s en release, contre environ 85 s
  pour SymPy, qui lève ensuite `NotImplementedError`. Reste : `together`
  natif sans `deep`, troisième appel légèrement plus lent que SymPy.
  Régression flottante du corpus EDO introduite par `efa251c` corrigée en
  vague 3 (`4ecdacc`) : les flottants sont rationalisés quand les inconnues
  sont des fonctions ou des dérivées, comme SymPy, et `expand` n'a plus de
  limite de taille par défaut (`max_terms` reste une garde explicite) ;
  corpus de nouveau à 296/263/259. Reste : forme `tanh` de `separable_26`.

Validation : comparer les familles de solutions et les constantes arbitraires,
vérifier les résidus et les conditions initiales, et contrôler la complétude
des bases de solutions. Une solution particulière valide ne suffit pas.

Voir [équations différentielles](docs/differential-equations.md).

## 5. Précision et évaluation numérique

- [x] `evalf(strict=True)` et `maxn` : lèvent `PrecisionExhausted` avec le
  message de référence (tests/test_evalf_options.py).
- [x] Évaluation multiprécision des fonctions spéciales, Meijer G compris, par
  mpmath, branches `exp_polar` comme SymPy (tests/test_special_precision.py).
- [ ] Garanties de précision près des singularités, coupures et annulations.
  Fait pour les fonctions élémentaires (tests/test_elementary_precision.py).
  Reste : auditer la continuation analytique de chaque fonction multivaluée
  hors de ces chemins. Vague 3 (`ad71575`) : une somme qui s'annule exactement,
  par exemple `-I*exp(3*I/5)/2 - exp(3*I/10)*sin(3/10) + I/2`, donne les
  flottants de précision 1 de SymPy (`0.e-131 + 0.e-131*I`, SymPy
  `0.e-130 + 0.e-130*I`) ; `strict` lève toujours.
- [x] Compteurs d'exposants : exposants binaires de taille quelconque depuis la
  vague 8 (`801b8ba`), stockés jusqu'à 65536 bits ; au-delà le résultat reste
  non évalué, et la partie entière d'un flottant géant lève `OverflowError`
  comme SymPy. Reste : `exp(2**70000).evalf()` non évalué (mpmath n'a pas de
  limite), affichage lent près de la limite.
  **Résultat faux, fermé (vague 3, `927d768`)** : `Pow(2, Rational(10**20 + 1,
  2))` donnait `sqrt(2)` ; le résultat reste exact comme dans
  SymPy (`tests/test_core_misc_wave3_parity.py`). `expand(complex=True)` n'a
  plus de limite de taille par défaut (`as_real_imag_bounded` reste une garde
  réservée aux appelants internes).

Validation : comparer la précision effectivement obtenue, et pas seulement le
nombre de chiffres affichés. Une fonction sans évaluateur multiprécision doit
rester symbolique ou échouer explicitement, sans fabriquer des chiffres depuis
une valeur en double précision.

Voir [écarts numériques actuels](docs/algorithmic-parity.md) et
[fonctions spéciales](docs/special-functions.md).

## 6. Simplification et calcul matriciel

- [ ] Élargir la recherche de candidats de `simplify` et les stratégies avancées
  de simplification, dont la couverture reste plus petite que celle de SymPy.
  La simplification native est une recherche de candidats, pas la chaîne de
  réécritures complète de SymPy.
- [ ] Étendre le traitement des expressions non commutatives et des options
  avancées sur les opérations où la référence les accepte. Vague 2 : symboles
  et fonctions non commutatifs dans le noyau (`Mul.flatten` porté, `args_cnc`,
  affichage, `subs`, `expand`, `diff`, pickle), simplificateurs dans l'ordre de
  SymPy (`factor`, `together`, `cancel`, `apart`, `nc_simplify`, `simplify`),
  `Poly` lève `PolynomialError` ; 784 tests dans
  `tests/test_noncommutative_parity.py`. Reste : `adjoint`/`transpose` lèvent
  `NotImplementedError`, forme de `radsimp`, `nonlinsolve` lève `ValueError`
  (SymPy `{(B*A**-1,)}`), `InconsistentAssumptions`, `srepr` des fonctions non
  commutatives, `limit`/`nsimplify`/`evalf`/matrices non vérifiés.
  Vague 3 (`tests/test_noncommutative_rest_parity.py`, 70 puis 14 différences
  sur 143) : `adjoint`/`transpose`, fin de `radsimp`, `nonlinsolve`,
  `InconsistentAssumptions`, `srepr` ; `Mul(-1, somme non commutative)`
  conservé (`tests/test_core_misc_wave3_parity.py`). Correction
  intentionnelle : `nonlinsolve([A*B - C], [B])` donne `{(A**-1*C,)}`. Reste :
  contradiction de `solveset(A*x - B)`, `Matrix.inv`, `adjoint(sin(A))`,
  `latex`.
  **Corrections intentionnelles de `solve` validées le 22 octobre 2026** : `solve(B*A - C, B)`
  donne `[C*A**-1]` (SymPy `[A**-1*C]`, ordre faux), `solve(A*B*C - 1, B)`
  `[(C*A)**-1]` (SymPy `[(A*C)**-1]`), `solve(A*B - B, B)` `[0]` (SymPy `[]`),
  `solve(A*B - B*A, A)` lève `NotImplementedError` (SymPy `[]`, alors que
  `A = B` convient), `solve(x**2*A - A, x)` résolu (SymPy lève),
  `integrate(A*x, x)` calculé (SymPy lève).
- [ ] Améliorer la reconnaissance exacte des relations entre fonctions lors
  des tests de nullité, du calcul de rang et du choix des pivots.
- [ ] Compléter l'inférence des domaines et les distinctions selon les paramètres
  dans les calculs matriciels, au-delà des régressions déjà corrigées.
- [ ] Étendre les valeurs propres et fonctions de matrices à coefficients
  flottants au-delà des cas spectraux et des petites matrices antisymétriques
  actuellement couverts.
- [x] Divers relevés pendant la vague 1 (vague 2, `tests/test_misc_wave2_parity.py`) :
  tri des arguments de `Min`/`Max` à l'affichage et dérivée en `Heaviside`,
  `polarify` porté (paramètres de G et conditions inchangés), séries d'un
  `Piecewise` de premier niveau, `Order.subs` garde le point limite, `y < y`
  décidé seulement pour des côtés réels étendus, `as_expr_set_pairs` porté,
  `piecewise_fold` porté (Python). Reste : terme d'ordre par branche
  (`O(x**3)` au lieu de `O(x**5)`), `Piecewise` imbriqué dans une somme ou un
  produit, `Pow` d'un `Piecewise` (`_eval_power` absent), `Order.subs` limité
  aux substitutions affines, ordre du dictionnaire de `polarify(subs=True)`.

Validation : contrôler les identités exactes, les matrices singulières,
les spécialisations de paramètres et les erreurs attendues. Évaluer aussi la
croissance des expressions ; une accélération ne doit pas supprimer un repli.

Voir [simplification](docs/compatibility.md), [matrices](docs/matrices.md) et
[tests de nullité matriciels](docs/matrix-zero-tests.md).

## Validation transversale avant toute affirmation de parité

- [ ] Transformer les écarts reproduits en régressions différentielles et
  compléter les familles de cas autour de chacun, sans se limiter aux exemples.
- [ ] Rejouer les audits sur la version native effectivement livrée ; les preuves
  historiques liées à d'autres empreintes ne certifient pas un moteur modifié.
- [ ] Trier les différences d'observation : affichage, résultat équivalent,
  domaine absent, solution manquante, précision insuffisante ou erreur différente.
- [ ] Tenir séparées les corrections mathématiques intentionnelles de bugs de
  SymPy et les correspondances avec les résultats de référence.
- [ ] Auditer les limites de ressources et les complexités en fonction des degrés,
  tailles de coefficients, graphes d'expressions et tailles des sorties.
  Un budget épuisé ne prouve jamais l'absence de solution.
- [ ] Mettre à jour ce fichier et le contrat détaillé à chaque écart fermé,
  avec les preuves et la portée exacte de la validation.

L'[audit versionné](docs/sympy-parity-audit.md) enregistre **2 916 objets rejoués**,
dont **1 363 aux observations différentes à examiner** dans les exemples
ordinaires. Il contient plusieurs instantanés natifs : ce nombre n'est ni un
nombre de bugs actuels, ni un pourcentage de capacité manquante. La parité
algorithmique globale n'est pas certifiée.
