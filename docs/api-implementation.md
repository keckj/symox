# API internes ajoutées : comportements et limites

[Inventaire complet](sympy-api.md) · [État par section](api-connections.md) · [README](../README.md)

Ces implémentations utilisent les expressions natives de Symox. Elles n’importent pas SymPy pour effectuer leurs calculs. La présence d’un nom dans l’inventaire ne garantit pas toute la compatibilité de son protocole avec SymPy.

## Programmation linéaire et logique

`solvers.simplex` implémente un simplexe exact à deux phases avec pivots de Bland. `linprog` accepte les inégalités, égalités et bornes ; `lpmin` et `lpmax` convertissent les relations symboliques et traitent les variables libres. Les problèmes infaisables ou non bornés lèvent les exceptions correspondantes. Les coefficients doivent permettre les comparaisons exactes nécessaires aux pivots ; ce n’est pas un solveur de contraintes à paramètres symboliques indécidables.

`PropKB` conserve les clauses CNF et utilise le solveur SAT local pour les requêtes. `pl_true` gère aussi les affectations partielles. Ces opérations concernent la logique propositionnelle.

## Polynômes et matrices de domaines

Les cinq familles de polynômes d’Appell sont construites exactement, avec retour facultatif d’un `Poly`. DDM, SDM et les fonctions `ddm_*`/`sdm_*` effectuent des opérations sur leurs véritables représentations denses ou creuses et réutilisent `DomainMatrix` pour les algorithmes exacts.

**DFM utilise le moteur Symox, sans FLINT.** L’existence de ce conteneur ne signifie pas la disponibilité des optimisations ou de tous les protocoles du backend FLINT de SymPy. Les algorithmes de corps et d’anneaux restent distincts ; les méthodes de stockage non inventoriées ne sont pas toutes adaptées.

`polyconfig` fournit une configuration isolée par contexte, des valeurs par défaut et des surcharges temporaires. Les sélecteurs d’algorithmes ne contrôlent pas tous les chemins de calcul du moteur Rust.

Les morphismes de modules utilisent des matrices exactes pour leurs images, noyaux et compositions ; le calcul de polynôme minimal repose sur les puissances de Krylov. Les recherches finies modulo un entier sont bornées. Ces interfaces ne remplacent pas tous les algorithmes d’idéaux ou de modules de rang quelconque de SymPy.

## Géométrie, matrices et ensembles

Les lieux géométriques concernés héritent de `GeometrySet` et de `Set`. Les classes de lignes 2D/3D ont leurs intermédiaires dédiés. Les opérations géométriques utilisent les méthodes exactes existantes ; leur appartenance à `Set` n’ajoute pas une algèbre générale de tous les ensembles géométriques dans le moteur natif.

`CompanionMatrix` construit la matrice compagnon d’un polynôme monique univarié. `MatrixSet` vérifie les dimensions et l’appartenance de chaque coefficient. `MatrixKind` et `SetKind` fournissent des métadonnées internées ; les anciens discriminants internes du graphe matriciel natif gardent leur représentation actuelle.

Les helpers de matrices par blocs réutilisent les simplifications natives et des compléments de Schur. L’arithmétique des ensembles distingue les opérations élément par élément de l’union et de l’intersection. Elle énumère les ensembles finis, applique les règles d’intervalles prises en charge et conserve autrement une image symbolique.

## Génération de code et compilation

L’AST comprend types, variables, déclarations, affectations, boucles, appels, fonctions et retours, ainsi que des nœuds propres à C et Fortran. Les printers et `utilities.codegen` produisent du code pour C, Fortran, Python, Rust, Julia et Octave sur les constructions prises en charge. Les tests compilent et exécutent des routines C, Fortran et Rust, et exécutent le Python produit.

Les optimisations numériques incluent `expm1`, `log1p`, les puissances et des transformations de sommes d’exponentielles. Les intrinsics symboliques conservent leur intention pour le printer ; leurs dérivées et réécritures ne reproduisent pas tous les protocoles internes de SymPy. Les approximations polynomiales restent soumises aux bornes et vérifications prévues par leurs méthodes.

`autowrap(..., backend='f2py')` compile une extension avec NumPy f2py et un compilateur Fortran. Le backend Cython nécessite l’installation optionnelle de Cython. `ufuncify` compile une extension C fournissant une véritable `numpy.ufunc`, avec diffusion NumPy et plusieurs sorties, en double précision. Les signatures scalaires prises en charge ne constituent pas une adaptation générale des tableaux indexés et des types complexes. Les erreurs de compilation incluent les diagnostics du compilateur.

`implemented_function` peut combiner des fonctions Python avec les expressions numériques natives, y compris les appels imbriqués. `lambdastr` produit des lambdas exécutables avec arguments imbriqués. Les fonctions personnalisées doivent fournir une implémentation compatible avec le backend numérique choisi.

## Tableaux abstraits et dérivées tensorielles

`ArraySymbol`, les produits tensoriels, contractions, diagonales, permutations, sommes et changements de forme conservent leurs dimensions sans matérialiser les entrées. Les opérations explicites réutilisent les tableaux natifs et leur limite de 100 000 entrées. Les contractions de taille symbolique restent des sommes exactes.

Les entrées abstraites utilisent des applications symboliques `ArrayElement` pour conserver la dépendance envers les indices dans les calculs natifs. Leur représentation interne diffère des objets `Indexed` de SymPy. Les conversions couvrent les sommes et produits matriciels reconnus, les transpositions, produits de Kronecker et accès indexés directs ; elles ne reconnaissent pas toute expression d’Einstein arbitraire.

`PartialDerivative` calcule les indices libres et contractés et évalue les dérivées à partir des composantes fournies. Les composantes de dérivation doivent être des coordonnées différentiables ; l’algèbre symbolique de produits contenant ces nouveaux objets n’est pas entièrement adaptée.

## Espaces de probabilité

Les domaines sont séparés de leurs mesures. Les domaines simples, produits et conditionnels finis ou continus fournissent des prédicats et des intégrations ou énumérations exactes. `pspace` retourne les espaces spécialisés pour les variables scalaires prises en charge.

Les espaces finis généraux conservent leurs poids joints pour les probabilités, espérances et tirages. Les espaces continus généraux utilisent leur densité jointe et peuvent garder des intégrales non évaluées. L’échantillonnage conditionnel par inversion peut être coûteux ; il demande des supports en intervalles et des paramètres numériques. Les restrictions générales décrites dans le [guide Statistics](statistics.md) restent applicables.

## Validation

Les tests dédiés se trouvent dans `tests/test_remaining_api_algorithms.py`, `tests/test_remaining_codegen.py` et `tests/test_remaining_structural_api.py`. Ils vérifient des résultats mathématiques, les erreurs, les formes des tableaux, les dépendances probabilistes et l’exécution du code compilé. Le corpus différentiel de 273 expressions demeure un indicateur distinct de la présence des API.

La suite complète de cette passe compte 26 689 tests réussis. Deux tests ajoutés ensuite pour les noyaux entiers passent également ; les vérifications ciblées des modules modifiés restent vertes. Les tests conditionnés par des dépendances indisponibles restent sautés.
