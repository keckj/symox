# Implémentation de l’API SymPy, section par section

[Inventaire mesuré](sympy-api.md) · [Comportements et limites des nouveaux modules](api-implementation.md) · [README](../README.md)

Mesure du 2026-10-18 contre SymPy 1.14.0. L’audit initial reste disponible dans le commit `e2a1816`.

311 chemins : 301 importables avec des exports mesurables, 0 manquants, 9 sans exports mesurables et 1 erreur dans la référence. 0 lignes ont encore des noms absents.

Les nombres mesurent la présence des noms. Ils ne vérifient ni les signatures ni les comportements. Les lignes se recouvrent : leurs nombres ne s’additionnent pas en une couverture globale. Les chemins sans exports et les erreurs de référence restent visibles.

| Section | Chemins | Chemins absents | Lignes avec noms absents |
|---|---:|---:|---:|
| [root](#root) | 1 | 0 | 0 |
| [abc](#abc) | 1 | 0 | 0 |
| [algebras](#algebras) | 1 | 0 | 0 |
| [assumptions](#assumptions) | 5 | 0 | 0 |
| [calculus](#calculus) | 5 | 0 | 0 |
| [categories](#categories) | 2 | 0 | 0 |
| [codegen](#codegen) | 13 | 0 | 0 |
| [combinatorics](#combinatorics) | 16 | 0 | 0 |
| [core](#core) | 23 | 0 | 0 |
| [crypto](#crypto) | 2 | 0 | 0 |
| [diffgeom](#diffgeom) | 1 | 0 | 0 |
| [discrete](#discrete) | 3 | 0 | 0 |
| [functions](#functions) | 14 | 0 | 0 |
| [geometry](#geometry) | 9 | 0 | 0 |
| [holonomic](#holonomic) | 1 | 0 | 0 |
| [integrals](#integrals) | 5 | 0 | 0 |
| [interactive](#interactive) | 3 | 0 | 0 |
| [liealgebras](#liealgebras) | 1 | 0 | 0 |
| [logic](#logic) | 2 | 0 | 0 |
| [matrices](#matrices) | 8 | 0 | 0 |
| [ntheory](#ntheory) | 14 | 0 | 0 |
| [parsing](#parsing) | 3 | 0 | 0 |
| [physics](#physics) | 68 | 0 | 0 |
| [plotting](#plotting) | 3 | 0 | 0 |
| [polys](#polys) | 17 | 0 | 0 |
| [printing](#printing) | 31 | 0 | 0 |
| [series](#series) | 2 | 0 | 0 |
| [sets](#sets) | 5 | 0 | 0 |
| [simplify](#simplify) | 13 | 0 | 0 |
| [solvers](#solvers) | 9 | 0 | 0 |
| [stats](#stats) | 6 | 0 | 0 |
| [tensor](#tensor) | 7 | 0 | 0 |
| [testing](#testing) | 4 | 0 | 0 |
| [utilities](#utilities) | 12 | 0 | 0 |
| [vector](#vector) | 1 | 0 | 0 |

<a id="root"></a>
<details>
<summary>root</summary>

Les 921 exports de premier niveau sont présents. Ce résultat ne garantit pas leurs signatures ni leurs comportements.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="abc"></a>
<details>
<summary>abc</summary>

Le module fournit les 76 symboles latins/grecs préconstruits et les dictionnaires _clash. Le compteur 0/0 ne mesure pas ces constantes : le chemin existe désormais, mais N/A reste N/A.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="algebras"></a>
<details>
<summary>algebras</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="assumptions"></a>
<details>
<summary>assumptions</summary>

Prédicats et relations binaires extensibles, faits CNF et helpers refine sont reliés au moteur natif. Les règles disponibles restent celles du moteur et des handlers enregistrés.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="calculus"></a>
<details>
<summary>calculus</summary>

Les 17 exports mesurés sont raccordés, dont AccumBounds et euler_equations. Les sous-modules mesurés sont également complets en noms.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="categories"></a>
<details>
<summary>categories</summary>

Aucun raccordement manquant dans les deux chemins mesurés.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="codegen"></a>
<details>
<summary>codegen</summary>

AST de programmes, types, déclarations, boucles, fonctions et nœuds C/Fortran ajoutés. Réécritures numériques, méthode de Newton et générateurs C, Fortran, Python, Rust, Julia et Octave présents. Voir le guide des comportements et limites.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="combinatorics"></a>
<details>
<summary>combinatorics</summary>

Constructeurs de groupes transitifs de degrés 5 et 6, nombres de groupes, Gray/RGS et helpers de stabilisateurs ajoutés. La canonisation tensorielle effectue une exploration exacte bornée des orbites signées.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="core"></a>
<details>
<summary>core</summary>

Types et métadonnées, hypothèses, singletons, facteurs, parcours et helpers numériques sont ajoutés. Les tuples de précision evalf adaptent le moteur natif ; ils ne constituent pas une certification par intervalles.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="crypto"></a>
<details>
<summary>crypto</summary>

Aucun raccordement manquant dans les deux chemins mesurés.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="diffgeom"></a>
<details>
<summary>diffgeom</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="discrete"></a>
<details>
<summary>discrete</summary>

Aucun raccordement manquant dans les trois chemins mesurés.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="functions"></a>
<details>
<summary>functions</summary>

Classes de familles de fonctions reliées aux expressions natives, décomposition réelle des fonctions prises en charge et corrections numériques de bêta incomplète.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="geometry"></a>
<details>
<summary>geometry</summary>

GeometrySet et les classes intermédiaires LinearEntity2D/3D sont intégrés à la hiérarchie. Matrices de rotation, translation et mise à l’échelle, recherche de symbole et coplanarité sont implémentées.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="holonomic"></a>
<details>
<summary>holonomic</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="integrals"></a>
<details>
<summary>integrals</summary>

Intégration exacte des polynômes sur polygones et polyèdres, moments et helpers de facettes ajoutés. Les adaptateurs Meijer utilisent l’intégration native et conservent les échecs ou conditions non résolus.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="interactive"></a>
<details>
<summary>interactive</summary>

Transformation des entiers et création automatique de symboles reliées aux sessions Python/IPython. Les hooks IPython nécessitent son installation.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="liealgebras"></a>
<details>
<summary>liealgebras</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="logic"></a>
<details>
<summary>logic</summary>

Les bases de connaissances propositionnelles utilisent la mise en forme CNF et le solveur SAT local. gateinputcount et les évaluations logiques partielles sont ajoutés.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="matrices"></a>
<details>
<summary>matrices</summary>

MatrixKind, MatrixSet, CompanionMatrix et les réécritures de matrices par blocs sont ajoutés. Les classes de représentation adaptent les matrices natives ; les nouveaux objets ne remplacent pas le moteur natif.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="ntheory"></a>
<details>
<summary>ntheory</summary>

digits est raccordé à ntheory.factor_. Aucun nom manquant dans les chemins mesurés.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="parsing"></a>
<details>
<summary>parsing</summary>

SymPyExpression analyse un sous-ensemble explicite des programmes numériques C et Fortran et génère un AST exécutable. Les constructions hors grammaire sont rejetées.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="physics"></a>
<details>
<summary>physics</summary>

Aucun nom manquant dans les chemins importables mesurés. Le chemin amont continuum_mechanics.Arch avec une majuscule échoue dans SymPy lui-même ; conserver N/A, ne pas inventer un faux export. La présence des noms ne supprime pas les limites du guide Physics.

- `symox.physics.continuum_mechanics.Arch` : N/A, import en échec dans la référence.

</details>

<a id="plotting"></a>
<details>
<summary>plotting</summary>

Validation des arguments et fabrique de tracés reliées aux objets Plot existants. Les dépendances graphiques restent optionnelles.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="polys"></a>
<details>
<summary>polys</summary>

Polynômes d’Appell exacts, configuration locale de contexte, stockage DDM/SDM/DFM, noyaux de calcul dense/creux, systèmes linéaires polynomiaux et morphismes de modules ajoutés. DFM utilise le moteur Symox : aucun backend FLINT n’est revendiqué.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="printing"></a>
<details>
<summary>printing</summary>

Printers Python/NumPy/mpmath, C89/C99/C11, Fortran et Rust ajoutés avec rendu des nœuds de programme. Protocoles pretty, conventions, précédences et mappings des autres langages sont reliés aux printers. Le backend Aesara demande son installation optionnelle.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="series"></a>
<details>
<summary>series</summary>

Aucun raccordement manquant dans les deux chemins inventoriés. Le rapport ne couvre pas tous les sous-modules internes.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="sets"></a>
<details>
<summary>sets</summary>

SetKind, les régions complexes cartésiennes/polaires et la normalisation angulaire sont ajoutés. L’arithmétique élément par élément utilise des règles d’intervalles, l’énumération finie ou des images symboliques.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="simplify"></a>
<details>
<summary>simplify</summary>

Helpers CSE, combinaisons de sommes/produits, radicaux et réduction trigonométrique par Gröbner ajoutés. Les opérateurs hypergéométriques et Meijer utilisent de véritables polynômes différentiels ; le catalogue de formules reste limité.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="solvers"></a>
<details>
<summary>solvers</summary>

Solveur linéaire exact à deux phases, avec contraintes d’égalité, bornes et variables libres, ajouté dans simplex. Helpers de solveset, constantes et conditions initiales, récurrences de séries ODE et résolution de systèmes non linéaires triangulaires sont reliés aux API. Les limites des séries de Frobenius et des systèmes couplés sont documentées.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="stats"></a>
<details>
<summary>stats</summary>

Domaines et espaces de probabilité finis/continus ajoutés, reliés à pspace et aux moteurs d’espérance/probabilité. Les lois jointes finies gardent leur dépendance pendant les calculs et l’échantillonnage. Les limites de Statistics restent applicables.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="tensor"></a>
<details>
<summary>tensor</summary>

Expressions abstraites de tableaux, conversions matricielles/indexées, opérations d’axes et dérivation par composantes ajoutées. Les tailles symboliques restent paresseuses ; la matérialisation est bornée. Les conversions indexées générales restent limitées.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="testing"></a>
<details>
<summary>testing</summary>

Adaptateurs pytest, doctest et tests aléatoires exécutables. Ils utilisent les runners Python disponibles et ne reproduisent pas toute l’infrastructure interne de tests SymPy.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="utilities"></a>
<details>
<summary>utilities</summary>

Mémoïsation, avertissements, utilitaires de source et de chaînes ajoutés. lambdastr et implemented_function sont reliés aux calculs numériques. codegen, autowrap/f2py et ufuncify NumPy compilent de véritables routines ; Cython reste une dépendance optionnelle.

Aucun nom absent dans les chemins mesurés.

</details>

<a id="vector"></a>
<details>
<summary>vector</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom absent dans les chemins mesurés.

</details>
