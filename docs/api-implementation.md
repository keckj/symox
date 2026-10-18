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

## Noyau, hypothèses et évaluation numérique

Les classes de nombres, types de fonctions, métadonnées de `Kind`, singletons et conteneurs identifient ou construisent les objets natifs correspondants. `StdFactKB` déduit les hypothèses via le moteur natif ; les relations binaires personnalisées peuvent enregistrer des handlers pour `ask`. Les helpers de facteurs, de parcours, de symboles et de cache opèrent sur les expressions locales.

Les helpers de `core.evalf` exposent le protocole de tuples mpmath à précision exprimée en bits. L’évaluation, les sommes et les quadratures utilisent les moteurs numériques disponibles. Les champs de précision sont des estimations de travail, sans certification générale d’erreur par intervalles. L’égalité exacte des `Float` tient compte de leur valeur binaire ; `all_close` fournit une comparaison à tolérance distincte.

## Combinaisons, groupes et tenseurs

Les groupes transitifs standards de degrés 5 et 6, les générateurs du Rubik’s Cube, les codes de Gray et les rangs RGS disposent de représentations effectives. Les tables de permutations et de glyphes sont des données statiques ; SymPy n’est pas chargé pour les utiliser. `find_transitive_subgroups_of_S6` vérifie les groupes du catalogue, sans reproduire une recherche aléatoire exhaustive.

`tensor_can` explore les orbites signées des symétries d’emplacements et des indices muets. Il détecte les annulations dues aux signes opposés et choisit un représentant canonique. L’exploration est limitée à 100 000 états ; ce n’est pas une implantation de l’algorithme optimisé Butler–Portugal. `canonicalize_naive` utilise le même moteur d’orbites. Le certificat de graphe utilise une recherche exhaustive bornée et son format diffère du certificat tensoriel de SymPy.

Les critères arithmétiques pour compter les groupes cycliques, abéliens ou nilpotents sont exacts. Le nombre général de groupes reste limité aux ordres couverts par les formules et tables disponibles ; les autres cas sont refusés.

## Séries d’ODE et simplification

`dsolve` accepte les séries de Taylor de premier ordre, les séries ordinaires linéaires du second ordre et les séries de Frobenius homogènes dont les deux branches se calculent par récurrence. Les calculs sont limités à 64 coefficients. Les racines indicielles répétées et les résonances nécessitant un logarithme lèvent `NotImplementedError`. Les conditions initiales prises en charge sont substituées dans les solutions et résolues pour les constantes. `classify_ode(..., n=...)` demande explicitement la classification des séries ; la classification par défaut conserve les familles de solutions exactes.

Les systèmes non linéaires triangulaires de deux ou trois fonctions sont réduits à des ODE scalaires, avec substitution des solutions et constantes distinctes. Les systèmes non linéaires couplés généraux restent hors de ce solveur. Les inspecteurs internes du second ordre ne reconnaissent aucune famille non linéaire, comme leurs homologues de référence.

Les helpers CSE ordonnent les dépendances, suivent les arguments communs et libèrent les temporaires après leur dernière utilisation. Les combinaisons de sommes et produits traitent les limites identiques ou adjacentes. `nthroot` accepte une simplification algébrique candidate après vérification exacte de sa puissance. `inversecombine` est une transformation formelle qui ignore volontairement les restrictions de branches.

La réduction trigonométrique par Gröbner utilise les relations polynomiales de sinus/cosinus ou sinh/cosh et compare deux ordres de générateurs, avec au plus huit arguments distincts. Les angles commensurables sont ramenés à une variable commune avant réduction lorsque leurs multiplicateurs entiers ne dépassent pas 16. Les options avancées `hints` et `polynomial` ne pilotent pas encore une stratégie équivalente à celle de SymPy.

Les opérateurs contigus hypergéométriques et Meijer sont de vrais polynômes en l’opérateur d’Euler. Les décalages inverses proviennent de la division par l’opérateur affine dans l’équation différentielle ; les annulations singulières sont refusées. Les plans hypergéométriques recherchent un chemin monotone évitant ces singularités, dans une borne de 256 pas et 10 000 états. Les paramètres de chaque groupe sont limités à huit pour leur appariement.

`Formula` construit une base différentielle et sa matrice compagnon. Les classes de formules sont des adaptateurs Python, sans toute la hiérarchie `Expr` de SymPy. Le catalogue élémentaire et la reconnaissance Lerch sont restreints ; la recherche utilise aussi les développements déjà connus du moteur natif. Les séries terminées et certains décalages de sommes sont calculés exactement. La présence de ces classes ne signifie pas que tout le catalogue hypergéométrique de SymPy a été porté.

## Intégration polynomiale et analyse de programmes

`polytope_integrate` calcule les moments par intégrales de frontière et triangulation de facettes, avec arithmétique exacte en dimensions deux et trois. Il accepte les polygones orientés, les polyèdres à facettes et les demi-espaces bornés du plan. Le degré maximal est 64. L’orientation détermine le signe ; les entrées doivent décrire un domaine et des facettes cohérents.

Les adaptateurs `meijerint` appellent l’intégration native. Un résultat non résolu retourne `None` ; les conditions connues restent explicites. Ils ne constituent pas un second moteur général d’intégration par Meijer G.

`SymPyExpression` analyse un sous-ensemble numérique de C et Fortran : déclarations, affectations, expressions, fonctions et boucles reconnues. L’AST peut produire du C, du Fortran ou du Python. Ce parseur n’implémente pas les grammaires complètes, le préprocesseur, l’édition de liens ni toutes les règles de types et de portée de ces langages. Les programmes hors du sous-ensemble doivent être adaptés ou traités par un compilateur spécialisé.

## Affichage, sessions et runners

Les conventions de noms, tables de précédence, dispositions pretty avec ligne de base et mappings de langages sont reliés aux printers. Les expressions symboliques gardent les limites de rendu propres à chaque backend. Le printer Rust traite les déclarations et calculs scalaires entiers, les plages de boucles et les conversions d’entiers vers flottants dans les fonctions testées ; il ne remplace pas un vérificateur de types Rust général.

Les sessions Python peuvent transformer les littéraux entiers et créer les symboles chargés inconnus. Les sessions IPython utilisent leurs hooks de transformation lorsque la dépendance est installée. `AesaraPrinter` construit des graphes Aesara si le paquet compatible est disponible ; il n’installe pas ce paquet et échoue explicitement sinon.

`testing` réutilise pytest et les doctests Python. Les interfaces exécutent réellement les tests, gèrent les assertions, avertissements et sous-processus pris en charge, sans embarquer le moteur de collecte et de réécriture propre à SymPy. Les sélecteurs et événements internes non adaptés ne donnent pas accès à toutes les fonctions du runner SymPy.

## Validation

Les fichiers `tests/test_remaining_*` vérifient les résultats exacts, les dépendances et dimensions, les erreurs, les valeurs numériques indépendantes des opérateurs spéciaux et l’exécution du code généré. Les tests de code compilent des programmes C, Fortran et Rust ; certains dépendent des outils optionnels installés.

Le corpus différentiel de 273 expressions demeure un indicateur distinct de la présence des API. Le tableau généré ne mesure ni toutes les signatures ni tous les comportements de SymPy.

La suite complète de cette passe donne **26 874 tests réussis, aucun échec**, sur 26 927 cas collectés et répartis en quatre lots disjoints. Les cas conditionnels sautés et les échecs attendus restent signalés par pytest. Un processus bloquant tout import de SymPy charge également les 310 chemins natifs importables et exécute des calculs des nouvelles interfaces.
