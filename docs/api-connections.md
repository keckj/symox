# Raccordements de l’API SymPy, section par section

[Inventaire mesuré](sympy-api.md) · [README](../README.md)

État après raccordements du 2026-10-18, contre SymPy 1.14.0. L’audit initial est conservé dans le commit `e2a1816`.

Validation : 43 nouveaux tests API, dont les imports dans les deux ordres sans SymPy au runtime, et 26 609 tests Python réussis sur la suite complète. Le corpus différentiel reste à 273/273.

311 chemins examinés : 257 importables avec des exports mesurables, 49 chemins manquants, 4 modules sans exports mesurables et 1 erreur d’import de référence. 129 lignes ont encore des noms manquants. Les lignes se recouvrent : ne pas additionner leurs noms pour calculer une couverture globale.

Un **raccordement** désigne un objet déjà exposé par le moteur ou un module local, à réexporter sous le chemin compatible. Les objets raccordés ont été retrouvés dans les implémentations locales, puis sélectionnés par famille ; leur existence ne prouve pas la compatibilité des signatures ou du comportement. La rubrique « adaptation / développement » inclut les objets sans source directe identifiée : certains peuvent être construits à partir de méthodes existantes. Les homonymes sans rapport sont exclus.

Les raccordements retenus dans l’audit sont appliqués. Les 54 chemins rétablis réutilisent les implémentations existantes ; les chemins et noms encore absents sont détaillés ci-dessous. Les signatures, hiérarchies et comportements ne sont pas entièrement compatibles pour autant.

Pour chaque façade, vérifier l’identité des objets réexportés et les deux ordres d’import. Les noms `simplify`, `series`, `lambdify`, `ask` ou `refine` peuvent désigner à la fois une fonction et un module : un nouvel import ne doit pas remplacer accidentellement une fonction publique par un module. Ne pas importer tous les sous-modules depuis la racine pour faire monter le compteur.

| Section | Chemins | Chemins absents | Lignes avec noms manquants |
|---|---:|---:|---:|
| [root](#root) | 1 | 0 | 0 |
| [abc](#abc) | 1 | 0 | 0 |
| [algebras](#algebras) | 1 | 0 | 0 |
| [assumptions](#assumptions) | 5 | 1 | 4 |
| [calculus](#calculus) | 5 | 0 | 0 |
| [categories](#categories) | 2 | 0 | 0 |
| [codegen](#codegen) | 13 | 11 | 13 |
| [combinatorics](#combinatorics) | 16 | 4 | 11 |
| [core](#core) | 23 | 1 | 18 |
| [crypto](#crypto) | 2 | 0 | 0 |
| [diffgeom](#diffgeom) | 1 | 0 | 0 |
| [discrete](#discrete) | 3 | 0 | 0 |
| [functions](#functions) | 14 | 0 | 3 |
| [geometry](#geometry) | 9 | 0 | 3 |
| [holonomic](#holonomic) | 1 | 0 | 0 |
| [integrals](#integrals) | 5 | 3 | 2 |
| [interactive](#interactive) | 3 | 0 | 1 |
| [liealgebras](#liealgebras) | 1 | 0 | 0 |
| [logic](#logic) | 2 | 0 | 2 |
| [matrices](#matrices) | 8 | 1 | 7 |
| [ntheory](#ntheory) | 14 | 0 | 0 |
| [parsing](#parsing) | 3 | 1 | 1 |
| [physics](#physics) | 68 | 0 | 0 |
| [plotting](#plotting) | 3 | 0 | 1 |
| [polys](#polys) | 17 | 7 | 10 |
| [printing](#printing) | 31 | 5 | 22 |
| [series](#series) | 2 | 0 | 0 |
| [sets](#sets) | 5 | 0 | 3 |
| [simplify](#simplify) | 13 | 1 | 6 |
| [solvers](#solvers) | 9 | 1 | 4 |
| [stats](#stats) | 6 | 2 | 3 |
| [tensor](#tensor) | 7 | 2 | 3 |
| [testing](#testing) | 4 | 4 | 3 |
| [utilities](#utilities) | 12 | 5 | 9 |
| [vector](#vector) | 1 | 0 | 0 |

<a id="root"></a>
<details>
<summary>root</summary>

Les 921 exports de premier niveau sont présents. Ce résultat ne garantit pas leurs signatures ni leurs comportements.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="abc"></a>
<details>
<summary>abc</summary>

Le module fournit les 76 symboles latins/grecs préconstruits et les dictionnaires _clash. Le compteur 0/0 ne mesure pas ces constantes : le chemin existe désormais, mais N/A reste N/A.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="algebras"></a>
<details>
<summary>algebras</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="assumptions"></a>
<details>
<summary>assumptions</summary>

register_handler et remove_handler sont raccordés dans le paquet et ask. Les relations binaires, métaclasses de prédicats, tables de faits et fonctions refine_* nécessitent encore des interfaces supplémentaires.

**`symox.assumptions` — 10/12 actuellement ; `available`.**

- Adaptation / développement à examiner : `AppliedBinaryRelation`, `BinaryRelation`.

**`symox.assumptions.ask` — 4/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `get_all_known_facts`, `get_known_facts_dict`.

**`symox.assumptions.assume` — 4/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `PredicateMeta`, `UndefinedPredicate`.

**`symox.assumptions.predicates` — 0/0 actuellement ; `missing_module`.**

- Aucun nom compté ; vérifier le chemin et les objets non mesurés (constantes ou paquet de documentation).

**`symox.assumptions.refine` — 1/9 actuellement ; `available`.**

- Adaptation / développement à examiner : `refine_Pow`, `refine_abs`, `refine_arg`, `refine_atan2`, `refine_im`, `refine_matrixelement`, `refine_re`, `refine_sign`.

</details>

<a id="calculus"></a>
<details>
<summary>calculus</summary>

Les 17 exports mesurés sont raccordés, dont AccumBounds et euler_equations. Les sous-modules mesurés sont également complets en noms.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="categories"></a>
<details>
<summary>categories</summary>

Aucun raccordement manquant dans les deux chemins mesurés.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="codegen"></a>
<details>
<summary>codegen</summary>

Travail de développement majoritaire : AST de programme, types, déclarations, boucles, appels, retours, nœuds C/Fortran et réécritures. Les générateurs CUDA/PyTorch/JAX existants ne fournissent pas à eux seuls ce modèle AST. Les homonymes Element, reshape, shape, use ou log2 ne constituent pas des raccordements valides.

**`symox.codegen` — 2/13 actuellement ; `available`.**

- Adaptation / développement à examiner : `Attribute`, `Declaration`, `For`, `FunctionCall`, `FunctionDefinition`, `FunctionPrototype`, `Print`, `Scope`, `Variable`, `While`, `aug_assign`.

**`symox.codegen.algorithms` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `newtons_method`, `newtons_method_function`.

**`symox.codegen.approximations` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `SeriesApprox`, `SumApprox`.

**`symox.codegen.ast` — 3/43 actuellement ; `available`.**

- Adaptation / développement à examiner : `AddAugmentedAssignment`, `AssignmentBase`, `Attribute`, `AugmentedAssignment`, `BreakToken`, `CodegenAST`, `ComplexBaseType`, `ComplexType`, `ContinueToken`, `Declaration`, `DivAugmentedAssignment`, `Element`, `FloatBaseType`, `FloatType`, `For`, `FunctionCall`, `FunctionDefinition`, `FunctionPrototype`, `IntBaseType`, `ModAugmentedAssignment`, `MulAugmentedAssignment`, `Node`, `NoneToken`, `Pointer`, `Print`, `QuotedString`, `Raise`, `Return`, `RuntimeError_`, `Scope`, `SignedIntType`, `Stream`, `String`, `SubAugmentedAssignment`, `Token`, `Type`, `UnsignedIntType`, `Variable`, `While`, `aug_assign`.

**`symox.codegen.cfunctions` — 0/11 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `Cbrt`, `Sqrt`, `exp2`, `expm1`, `fma`, `hypot`, `isinf`, `isnan`, `log10`, `log1p`, `log2`.

**`symox.codegen.cnodes` — 0/11 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `CommaOperator`, `Label`, `PostDecrement`, `PostIncrement`, `PreDecrement`, `PreIncrement`, `alignof`, `goto`, `sizeof`, `struct`, `union`.

**`symox.codegen.cutils` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `render_as_source_file`.

**`symox.codegen.cxxnodes` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `using`.

**`symox.codegen.fnodes` — 0/32 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `ArrayConstructor`, `Do`, `Extent`, `F95Function`, `FFunction`, `FortranReturn`, `GoTo`, `ImpliedDoLoop`, `Module`, `Program`, `Subroutine`, `SubroutineCall`, `allocated`, `array`, `bind_C`, `cmplx`, `dimension`, `dsign`, `isign`, `kind`, `lbound`, `literal_dp`, `literal_sp`, `merge`, `product_`, `reshape`, `shape`, `size`, `sum_`, `ubound`, `use`, `use_rename`.

**`symox.codegen.futils` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `render_as_module`.

**`symox.codegen.matrix_nodes` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `MatrixSolve`.

**`symox.codegen.pyutils` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `render_as_module`.

**`symox.codegen.rewriting` — 0/18 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `FuncMinusOneOptim`, `Optimization`, `ReplaceOptim`, `cosm1_opt`, `create_expand_pow_optimization`, `exp2_opt`, `expm1_opt`, `log1p_opt`, `log2_opt`, `log2const_opt`, `logaddexp2_opt`, `logaddexp_opt`, `logsumexp_2terms_opt`, `matinv_opt`, `optimize`, `powm1_opt`, `sinc_opt1`, `sinc_opt2`.

</details>

<a id="combinatorics"></a>
<details>
<summary>combinatorics</summary>

DirectProduct, Coset et SymmetricPermutationGroup sont raccordés aux modules voisins. Les groupes nommés supplémentaires, utilitaires Gray/RGS et algorithmes de canonicalisation tensorielle demandent encore des ajouts. Les alias Perm/PermGroup nécessitent un examen de leur protocole.

**`symox.combinatorics.galois` — 4/21 actuellement ; `available`.**

- Adaptation / développement à examiner : `A4_in_S6`, `A4xC2`, `G18`, `G36m`, `G36p`, `G72`, `M20`, `PGL2F5`, `PSL2F5`, `S3_in_S6`, `S4m`, `S4p`, `S4xC2`, `S5TransitiveSubgroups`, `S6TransitiveSubgroups`, `find_transitive_subgroups_of_S6`, `four_group`.

**`symox.combinatorics.generators` — 4/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `rubik`, `rubik_cube_generators`.

**`symox.combinatorics.graycode` — 1/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `bin_to_gray`, `get_subset_from_bitstring`, `gray_to_bin`, `graycode_subsets`, `random_bitstring`.

**`symox.combinatorics.group_numbers` — 0/4 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `groups_count`, `is_abelian_number`, `is_cyclic_number`, `is_nilpotent_number`.

**`symox.combinatorics.named_groups` — 6/8 actuellement ; `available`.**

- Adaptation / développement à examiner : `set_alternating_group_properties`, `set_symmetric_group_properties`.

**`symox.combinatorics.partitions` — 5/7 actuellement ; `available`.**

- Adaptation / développement à examiner : `RGS_generalized`, `random_integer_partition`.

**`symox.combinatorics.perm_groups` — 3/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `PermGroup`.

**`symox.combinatorics.permutations` — 2/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `AppliedPermutation`, `Perm`.

**`symox.combinatorics.subsets` — 1/2 actuellement ; `available`.**

- Adaptation / développement à examiner : `ksubsets`.

**`symox.combinatorics.tensor_can` — 0/12 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `bsgs_direct_product`, `canonical_free`, `canonicalize`, `double_coset_can_rep`, `dummy_sgs`, `gens_products`, `get_minimal_bsgs`, `get_symmetric_group_sgs`, `get_transversals`, `perm_af_direct_product`, `tensor_gens`, `transversal2coset`.

**`symox.combinatorics.testutil` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `canonicalize_naive`, `graph_certificate`.

**`symox.combinatorics.util` — 0/0 actuellement ; `missing_module`.**

- Aucun nom compté ; vérifier le chemin et les objets non mesurés (constantes ou paquet de documentation).

</details>

<a id="core"></a>
<details>
<summary>core</summary>

Les 101 exports du paquet sont présents. Les façades add, mul, power, symbol, relational, traversal, cache et expr sont raccordées aux objets existants. assumptions/common_assumptions interrogent les propriétés natives ; evalf reste un module et sympify une fonction, comme dans la référence. Les classes de nombres singleton, métaclasses et routines internes evalf ci-dessous restent à compléter.

**`symox.core.assumptions` — 4/8 actuellement ; `available`.**

- Adaptation / développement à examiner : `ManagedProperties`, `StdFactKB`, `as_property`, `make_property`.

**`symox.core.basic` — 2/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `as_Basic`.

**`symox.core.cache` — 3/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `lazy_function`, `print_cache`.

**`symox.core.containers` — 2/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `OrderedSet`, `TupleKind`, `tuple_wrapper`.

**`symox.core.evalf` — 3/44 actuellement ; `available`.**

- Adaptation / développement à examiner : `EvalfMixin`, `add_terms`, `as_mpmath`, `bitcount`, `check_convergence`, `check_target`, `chop_parts`, `complex_accuracy`, `do_integral`, `evalf_abs`, `evalf_add`, `evalf_alg_num`, `evalf_atan`, `evalf_ceiling`, `evalf_exp`, `evalf_float`, `evalf_floor`, `evalf_im`, `evalf_integer`, `evalf_integral`, `evalf_log`, `evalf_mul`, `evalf_piecewise`, `evalf_pow`, `evalf_prod`, `evalf_rational`, `evalf_re`, `evalf_subs`, `evalf_sum`, `evalf_symbol`, `evalf_trig`, `fastlog`, `finalize_complex`, `get_abs`, `get_complex_part`, `get_integer_part`, `hypsum`, `iszero`, `pure_complex`, `quad_to_mpmath`, `scaled_zero`.

**`symox.core.expr` — 4/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `ExprBuilder`, `unchanged`.

**`symox.core.exprtools` — 3/7 actuellement ; `available`.**

- Adaptation / développement à examiner : `Factors`, `Term`, `decompose_power`, `decompose_power_rat`.

**`symox.core.function` — 22/28 actuellement ; `available`.**

- Adaptation / développement à examiner : `Application`, `AppliedUndef`, `ArgumentIndexError`, `DefinedFunction`, `UndefSageHelper`, `UndefinedFunction`.

**`symox.core.intfunc` — 9/10 actuellement ; `available`.**

- Adaptation / développement à examiner : `igcd_lehmer`.

**`symox.core.kind` — 1/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `KindDispatcher`, `KindMeta`.

**`symox.core.mul` — 2/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `NC_Marker`, `expand_2arg`.

**`symox.core.multidimensional` — 1/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `apply_on_element`, `iter_copy`, `structure_copy`.

**`symox.core.numbers` — 14/34 actuellement ; `available`.**

- Adaptation / développement à examiner : `ComplexInfinity`, `Exp1`, `Half`, `ImaginaryUnit`, `Infinity`, `IntegerConstant`, `NaN`, `NegativeInfinity`, `NegativeOne`, `One`, `Pi`, `RationalConstant`, `Zero`, `all_close`, `equal_valued`, `int_valued`, `mpf_norm`, `sympify_complex`, `sympify_fractions`, `sympify_mpmath`.

**`symox.core.random` — 0/4 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `random_complex_number`, `seed`, `test_derivative_numerically`, `verify_numerically`.

**`symox.core.singleton` — 1/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `Singleton`, `SingletonRegistry`.

**`symox.core.symbol` — 5/8 actuellement ; `available`.**

- Adaptation / développement à examiner : `Str`, `disambiguate`, `uniquely_named_symbol`.

**`symox.core.sympify` — 2/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `CantSympify`, `kernS`.

**`symox.core.traversal` — 4/7 actuellement ; `available`.**

- Adaptation / développement à examiner : `iterargs`, `iterfreeargs`, `walk`.

</details>

<a id="crypto"></a>
<details>
<summary>crypto</summary>

Aucun raccordement manquant dans les deux chemins mesurés.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="diffgeom"></a>
<details>
<summary>diffgeom</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="discrete"></a>
<details>
<summary>discrete</summary>

Aucun raccordement manquant dans les trois chemins mesurés.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="functions"></a>
<details>
<summary>functions</summary>

Les 178 exports du paquet sont raccordés à la racine, sans figer les constructeurs avant leur installation finale. Les façades elliptic_integrals et zeta_functions sont complètes ; error_functions et polynomials gardent les classes de base et auxiliaires manquants ci-dessous.

**`symox.functions.special.beta_functions` — 3/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `betainc_mpmath_fix`.

**`symox.functions.special.error_functions` — 18/21 actuellement ; `available`.**

- Adaptation / développement à examiner : `FresnelIntegral`, `TrigonometricIntegral`, `real_to_real_as_real_imag`.

**`symox.functions.special.polynomials` — 13/14 actuellement ; `available`.**

- Adaptation / développement à examiner : `OrthogonalPolynomial`.

</details>

<a id="geometry"></a>
<details>
<summary>geometry</summary>

deg/rad sont raccordés dans polygon. Point3D.are_coplanar existe pour les points ; la fonction libre exige un adaptateur couvrant aussi les autres entités géométriques. GeometrySet, LinearEntity2D/3D et les fonctions libres de transformation demandent un examen de hiérarchie et de signatures. Ne pas réutiliser find des expressions pour find de geometry.util.

**`symox.geometry.entity` — 1/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `GeometrySet`, `rotate`, `scale`, `translate`.

**`symox.geometry.line` — 10/12 actuellement ; `available`.**

- Adaptation / développement à examiner : `LinearEntity2D`, `LinearEntity3D`.

**`symox.geometry.util` — 7/9 actuellement ; `available`.**

- Adaptation / développement à examiner : `are_coplanar`, `find`.

</details>

<a id="holonomic"></a>
<details>
<summary>holonomic</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="integrals"></a>
<details>
<summary>integrals</summary>

Le paquet (30/30) et transforms (32/32) sont raccordés. Les chemins integrals, laplace et singularityfunctions exposent également les fonctions existantes. intpoly et meijerint restent à développer ; leurs homonymes dans d’autres familles ne sont pas des remplacements.

**`symox.integrals.intpoly` — 0/22 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `best_origin`, `cross_product`, `decompose`, `distance_to_side`, `gradient_terms`, `hyperplane_parameters`, `integration_reduction`, `integration_reduction_dynamic`, `intersection`, `is_vertex`, `left_integral2D`, `left_integral3D`, `lineseg_integrate`, `main_integrate`, `main_integrate3d`, `norm`, `plot_polynomial`, `plot_polytope`, `point_sort`, `polygon_integrate`, `polytope_integrate`, `strip`.

**`symox.integrals.meijerint` — 0/3 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `meijerint_definite`, `meijerint_indefinite`, `meijerint_inversion`.

**`symox.integrals.meijerint_doc` — 0/0 actuellement ; `missing_module`.**

- Aucun nom compté ; vérifier le chemin et les objets non mesurés (constantes ou paquet de documentation).

</details>

<a id="interactive"></a>
<details>
<summary>interactive</summary>

interactive_traversal, init_printing et init_session sont raccordés. Les transformations automatiques de saisie et initialisations spécifiques IPython/Python demandent encore des adaptateurs de session.

**`symox.interactive.session` — 1/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `enable_automatic_int_sympification`, `enable_automatic_symbols`, `init_ipython_session`, `init_python_session`, `int_to_Integer`.

</details>

<a id="liealgebras"></a>
<details>
<summary>liealgebras</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="logic"></a>
<details>
<summary>logic</summary>

Les opérateurs booléens, formes normales et satisfiable sont raccordés dans logic, boolalg et inference. gateinputcount et les interfaces de base de connaissances KB/PropKB restent à développer.

**`symox.logic` — 19/20 actuellement ; `available`.**

- Adaptation / développement à examiner : `gateinputcount`.

**`symox.logic.inference` — 1/7 actuellement ; `available`.**

- Adaptation / développement à examiner : `KB`, `PropKB`, `entails`, `literal_symbol`, `pl_true`, `valid`.

</details>

<a id="matrices"></a>
<details>
<summary>matrices</summary>

Les réexports du paquet et les façades immutable, sparse, matrixbase et sparsetools sont raccordés. Le paquet atteint 70/73. MatrixKind/MatrixSet, CompanionMatrix, les classes intermédiaires et les réécritures blockmatrix restent à adapter ou implémenter.

**`symox.matrices` — 70/73 actuellement ; `available`.**

- Adaptation / développement à examiner : `MatrixKind`, `MatrixSet`, `dotprodsimp`.

**`symox.matrices.expressions` — 37/40 actuellement ; `available`.**

- Adaptation / développement à examiner : `CompanionMatrix`, `MatrixSet`, `combine_kronecker`.

**`symox.matrices.expressions.blockmatrix` — 4/16 actuellement ; `available`.**

- Adaptation / développement à examiner : `bc_block_plus_ident`, `bc_dist`, `bc_inverse`, `bc_matadd`, `bc_matmul`, `bc_transpose`, `bc_unpack`, `blockinverse_1x1`, `blockinverse_2x2`, `bounds`, `deblock`, `reblock_2x2`.

**`symox.matrices.immutable` — 3/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `ImmutableRepMatrix`, `sympify_matrix`, `sympify_mpmath_matrix`.

**`symox.matrices.kind` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `MatrixKind`.

**`symox.matrices.matrixbase` — 5/7 actuellement ; `available`.**

- Adaptation / développement à examiner : `a2idx`, `classof`.

**`symox.matrices.sparse` — 2/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `SparseRepMatrix`.

</details>

<a id="ntheory"></a>
<details>
<summary>ntheory</summary>

digits est raccordé à ntheory.factor_. Aucun nom manquant dans les chemins mesurés.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="parsing"></a>
<details>
<summary>parsing</summary>

SymPyExpression manque. Il faut un adaptateur de conversion/parsing de programmes, pas un alias de parse_expr.

**`symox.parsing.sym_expr` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `SymPyExpression`.

</details>

<a id="physics"></a>
<details>
<summary>physics</summary>

Aucun nom manquant dans les chemins importables mesurés. Le chemin amont continuum_mechanics.Arch avec une majuscule échoue dans SymPy lui-même ; conserver N/A, ne pas inventer un faux export. La présence des noms ne supprime pas les limites du guide Physics.

- `symox.physics.continuum_mechanics.Arch` : **N/A** ; erreur dans le module de référence, aucun raccordement natif à déduire.

</details>

<a id="plotting"></a>
<details>
<summary>plotting</summary>

La façade pygletplot expose PygletPlot existant avec son backend actuel ; cela ne crée pas un nouveau backend Pyglet. plot_factory et check_arguments nécessitent encore des adaptateurs.

**`symox.plotting.plot` — 10/12 actuellement ; `available`.**

- Adaptation / développement à examiner : `check_arguments`, `plot_factory`.

</details>

<a id="polys"></a>
<details>
<summary>polys</summary>

Les 81 exports disponibles à la racine sont raccordés au paquet (171/176). DM est raccordé à domainmatrix. Les cinq constructeurs de polynômes spéciaux, DDM/SDM/DFM, routines ddm_*/sdm_*, endomorphismes et configuration des algorithmes restent à compléter.

**`symox.polys` — 171/176 actuellement ; `available`.**

- Adaptation / développement à examiner : `andre_poly`, `bernoulli_c_poly`, `bernoulli_poly`, `euler_poly`, `genocchi_poly`.

**`symox.polys.matrices._dfm` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `DFM`.

**`symox.polys.matrices._typing` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `RingElement`.

**`symox.polys.matrices.ddm` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `DDM`.

**`symox.polys.matrices.dense` — 0/17 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `R`, `T`, `ddm_berk`, `ddm_iadd`, `ddm_idet`, `ddm_iinv`, `ddm_ilu`, `ddm_ilu_solve`, `ddm_ilu_split`, `ddm_imatmul`, `ddm_imul`, `ddm_ineg`, `ddm_irmul`, `ddm_irref`, `ddm_irref_den`, `ddm_isub`, `ddm_transpose`.

**`symox.polys.matrices.sdm` — 0/13 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `SDM`, `binop_dict`, `sdm_berk`, `sdm_dotvec`, `sdm_irref`, `sdm_matmul`, `sdm_matmul_exraw`, `sdm_matvecmul`, `sdm_nullspace_from_rref`, `sdm_particular_from_rref`, `sdm_rref_den`, `sdm_transpose`, `unop_dict`.

**`symox.polys.numberfields.modules` — 3/13 actuellement ; `available`.**

- Adaptation / développement à examiner : `EndomorphismRing`, `InnerEndomorphism`, `Module`, `ModuleElement`, `ModuleEndomorphism`, `ModuleHomomorphism`, `find_min_poly`, `is_sq_maxrank_HNF`, `make_mod_elt`, `to_col`.

**`symox.polys.polyconfig` — 0/4 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `configure`, `query`, `setup`, `using`.

**`symox.polys.polyoptions` — 1/2 actuellement ; `available`.**

- Adaptation / développement à examiner : `build_options`.

**`symox.polys.solvers` — 0/5 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `PolyNonlinearError`, `RawMatrix`, `eqs_to_matrix`, `solve_lin_sys`, `sympy_eqs_to_ring`.

</details>

<a id="printing"></a>
<details>
<summary>printing</summary>

Les fonctions print_* déjà disponibles, les exports pretty et la façade codeprinter sont raccordés. Les classes de printers, exceptions et protocoles communs, fonctions de symbologie et stringPict ci-dessous demandent encore une adaptation ou une implémentation.

**`symox.printing.aesaracode` — 0/4 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `AesaraPrinter`, `aesara_code`, `aesara_function`, `dim_handling`.

**`symox.printing.c` — 2/10 actuellement ; `available`.**

- Adaptation / développement à examiner : `C11CodePrinter`, `C89CodePrinter`, `C99CodePrinter`, `CodePrinter`, `get_math_macros`, `known_functions_C89`, `known_functions_C99`, `requires`.

**`symox.printing.codeprinter` — 6/11 actuellement ; `available`.**

- Adaptation / développement à examiner : `AssignmentError`, `CodePrinter`, `PrintMethodNotImplementedError`, `print_rust_code`, `requires`.

**`symox.printing.conventions` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `requires_partial`, `split_super_sub`.

**`symox.printing.fortran` — 2/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `FCodePrinter`.

**`symox.printing.jscode` — 3/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `known_functions`.

**`symox.printing.julia` — 2/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `known_fcns_src1`, `known_fcns_src2`, `print_julia_code`.

**`symox.printing.lambdarepr` — 0/6 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `IntervalPrinter`, `LambdaPrinter`, `MpmathPrinter`, `NumPyPrinter`, `PythonCodePrinter`, `lambdarepr`.

**`symox.printing.latex` — 4/7 actuellement ; `available`.**

- Adaptation / développement à examiner : `accepted_latex_functions`, `latex_escape`, `translate`.

**`symox.printing.mathematica` — 2/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `known_functions`.

**`symox.printing.mathml` — 5/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `MathMLPrinter`.

**`symox.printing.octave` — 2/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `known_fcns_src1`, `known_fcns_src2`, `print_octave_code`.

**`symox.printing.precedence` — 2/12 actuellement ; `available`.**

- Adaptation / développement à examiner : `PRECEDENCE_FUNCTIONS`, `PRECEDENCE_VALUES`, `precedence_Float`, `precedence_FracElement`, `precedence_Integer`, `precedence_Mul`, `precedence_PolyElement`, `precedence_Rational`, `precedence_UnevaluatedExpr`, `precedence_traditional`.

**`symox.printing.pretty.pretty` — 7/19 actuellement ; `available`.**

- Adaptation / développement à examiner : `U`, `annotated`, `center_pad`, `hobj`, `is_subscriptable_in_unicode`, `pretty_atom`, `pretty_symbol`, `pretty_try_use_unicode`, `pretty_use_unicode`, `vobj`, `xobj`, `xsym`.

**`symox.printing.pretty.pretty_symbology` — 0/25 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `G`, `U`, `VF`, `annotated`, `atoms_table`, `center`, `center_pad`, `digit_2txt`, `frac`, `g`, `greek_letters`, `greek_unicode`, `hobj`, `pretty_atom`, `pretty_symbol`, `pretty_try_use_unicode`, `pretty_use_unicode`, `root`, `sub`, `sup`, `symb_2txt`, `vobj`, `xobj`, `xstr`, `xsym`.

**`symox.printing.pretty.stringpict` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `prettyForm`, `stringPict`.

**`symox.printing.printer` — 1/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `print_function`, `printer_context`.

**`symox.printing.pycode` — 1/7 actuellement ; `available`.**

- Adaptation / développement à examiner : `AbstractPythonCodePrinter`, `ArrayPrinter`, `CmathPrinter`, `MpmathPrinter`, `PythonCodePrinter`, `SymPyPrinter`.

**`symox.printing.python` — 2/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `PythonPrinter`.

**`symox.printing.rcode` — 3/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `known_functions`.

**`symox.printing.rust` — 0/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `RustCodePrinter`, `TypeCast`, `float_ceiling`, `float_floor`, `known_functions`.

**`symox.printing.str` — 3/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `StrReprPrinter`.

</details>

<a id="series"></a>
<details>
<summary>series</summary>

Aucun raccordement manquant dans les deux chemins inventoriés. Le rapport ne couvre pas tous les sous-modules internes.

Aucun nom manquant dans les chemins mesurés.

</details>

<a id="sets"></a>
<details>
<summary>sets</summary>

Les 26 exports du paquet et les façades conditionset, fancysets, powerset et sets sont raccordés aux objets existants. SetKind, les variantes CartesianComplexRegion/PolarComplexRegion et les auxiliaires d’arithmétique des ensembles restent à compléter.

**`symox.sets.conditionset` — 2/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `SetKind`.

**`symox.sets.fancysets` — 9/12 actuellement ; `available`.**

- Adaptation / développement à examiner : `CartesianComplexRegion`, `PolarComplexRegion`, `normalize_theta_set`.

**`symox.sets.sets` — 12/22 actuellement ; `available`.**

- Adaptation / développement à examiner : `SetKind`, `is_function_invertible_in_set`, `set_add`, `set_div`, `set_function`, `set_mul`, `set_pow`, `set_sub`, `simplify_intersection`, `simplify_union`.

</details>

<a id="simplify"></a>
<details>
<summary>simplify</summary>

Les 32 exports du paquet et les 31 de fu, règles TR* comprises, sont raccordés. Les façades réutilisent _simptools, _epath et les fonctions natives ; elles suivent les adaptateurs installés ultérieurement. Les internes CSE, certaines réécritures et collections/formules hypergéométriques restent à développer. Aucun opérateur quantique homonyme n’a été réexporté.

**`symox.simplify.cse_main` — 1/11 actuellement ; `available`.**

- Adaptation / développement à examiner : `FuncArgTracker`, `Unevaluated`, `cse_release_variables`, `cse_separate`, `match_common_args`, `opt_cse`, `postprocess_for_cse`, `preprocess_for_cse`, `reps_toposort`, `tree_cse`.

**`symox.simplify.hyperexpand` — 1/37 actuellement ; `available`.**

- Adaptation / développement à examiner : `Formula`, `FormulaCollection`, `G_Function`, `Hyper_Function`, `MeijerFormula`, `MeijerFormulaCollection`, `MeijerShiftA`, `MeijerShiftB`, `MeijerShiftC`, `MeijerShiftD`, `MeijerUnShiftA`, `MeijerUnShiftB`, `MeijerUnShiftC`, `MeijerUnShiftD`, `MultOperator`, `Operator`, `ReduceOrder`, `ShiftA`, `ShiftB`, `UnShiftA`, `UnShiftB`, `add_formulae`, `add_meijerg_formulae`, `apply_operators`, `build_hypergeometric_formula`, `debug`, `devise_plan`, `devise_plan_meijer`, `hyperexpand_special`, `make_derivative_operator`, `make_simp`, `reduce_order`, `reduce_order_meijer`, `try_lerchphi`, `try_polynomial`, `try_shifted_sum`.

**`symox.simplify.hyperexpand_doc` — 0/0 actuellement ; `missing_module`.**

- Aucun nom compté ; vérifier le chemin et les objets non mesurés (constantes ou paquet de documentation).

**`symox.simplify.radsimp` — 7/17 actuellement ; `available`.**

- Adaptation / développement à examiner : `collect_abs`, `collect_sqrt`, `denom_expand`, `expand_denom`, `expand_fraction`, `expand_numer`, `fraction_expand`, `numer_expand`, `rad_rationalize`, `split_surds`.

**`symox.simplify.simplify` — 10/21 actuellement ; `available`.**

- Adaptation / développement à examiner : `clear_coefficients`, `dotprodsimp`, `factor_sum`, `inversecombine`, `nc_simplify`, `nthroot`, `product_mul`, `product_simplify`, `sum_add`, `sum_combine`, `sum_simplify`.

**`symox.simplify.sqrtdenest` — 1/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `SqrtdenestStopIteration`, `is_algebraic`, `is_sqrt`, `sqrt_biquadratic_denest`, `sqrt_depth`.

**`symox.simplify.trigsimp` — 3/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `trigsimp_groebner`, `trigsimp_old`.

</details>

<a id="solvers"></a>
<details>
<summary>solvers</summary>

Les fonctions disponibles de récurrence, systèmes polynomiaux, décomposition et inégalités sont raccordées (43/46 au niveau du paquet). simplex et ses erreurs nécessitent un véritable solveur LP. Les internes des ODE/systèmes et les inversions de solveset restent à compléter.

**`symox.solvers` — 43/46 actuellement ; `available`.**

- Adaptation / développement à examiner : `linprog`, `lpmax`, `lpmin`.

**`symox.solvers.ode.ode` — 5/23 actuellement ; `available`.**

- Adaptation / développement à examiner : `check_linear_2eq_order1`, `check_nonlinear_2eq_order1`, `check_nonlinear_2eq_order2`, `check_nonlinear_3eq_order1`, `check_nonlinear_3eq_order2`, `classify_sysode`, `constant_renumber`, `get_numbered_constants`, `iter_numbered_constants`, `ode_1st_power_series`, `ode_2nd_power_series_ordinary`, `ode_2nd_power_series_regular`, `ode_sol_simplicity`, `odesimp`, `solve_ics`, `sysode_linear_2eq_order1`, `sysode_nonlinear_2eq_order1`, `sysode_nonlinear_3eq_order1`.

**`symox.solvers.simplex` — 0/6 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `InfeasibleLPError`, `UnboundedLPError`, `linprog`, `lpmax`, `lpmin`, `show_linprog`.

**`symox.solvers.solveset` — 8/14 actuellement ; `available`.**

- Adaptation / développement à examiner : `domain_check`, `invert_complex`, `invert_real`, `linear_coeffs`, `solve_decomposition`, `solvify`.

</details>

<a id="stats"></a>
<details>
<summary>stats</summary>

Les distributions et fonctions publiques existent, mais les abstractions ContinuousDomain/PSpace, FiniteDomain/PSpace et leurs produits/conditions ne sont pas implémentées sous ce modèle. Il faut une couche objet compatible reliée au moteur existant ; aliaser toutes ces classes vers Distribution serait incorrect.

**`symox.stats.crv` — 0/9 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `ConditionalContinuousDomain`, `ContinuousDistribution`, `ContinuousDomain`, `ContinuousPSpace`, `ProductContinuousDomain`, `SingleContinuousDistribution`, `SingleContinuousDomain`, `SingleContinuousPSpace`, `reduce_rational_inequalities_wrap`.

**`symox.stats.frv` — 0/9 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `ConditionalFiniteDomain`, `FiniteDensity`, `FiniteDomain`, `FinitePSpace`, `ProductFiniteDomain`, `ProductFinitePSpace`, `SingleFiniteDistribution`, `SingleFiniteDomain`, `SingleFinitePSpace`.

**`symox.stats.rv` — 27/41 actuellement ; `available`.**

- Adaptation / développement à examiner : `ConditionalDomain`, `Density`, `IndependentProductPSpace`, `MatrixDomain`, `NamedArgsMixin`, `ProductDomain`, `ProductPSpace`, `RandomDomain`, `RandomIndexedSymbol`, `RandomMatrixSymbol`, `SingleDomain`, `SinglePSpace`, `sample_iter_lambdify`, `sample_iter_subs`.

</details>

<a id="tensor"></a>
<details>
<summary>tensor</summary>

Les 19 exports du paquet, les objets indexés et les fonctions d’indices sont raccordés. Les expressions abstraites ArraySymbol/ArrayContraction/ArrayTensorProduct et leurs conversions sont distinctes des tableaux concrets existants. PartialDerivative tensoriel reste aussi à développer.

**`symox.tensor.array.expressions` — 0/16 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `ArrayAdd`, `ArrayContraction`, `ArrayDiagonal`, `ArrayElement`, `ArrayElementwiseApplyFunc`, `ArraySymbol`, `ArrayTensorProduct`, `OneArray`, `PermuteDims`, `Reshape`, `ZeroArray`, `array_derive`, `convert_array_to_indexed`, `convert_array_to_matrix`, `convert_indexed_to_array`, `convert_matrix_to_array`.

**`symox.tensor.indexed` — 3/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `IndexException`.

**`symox.tensor.toperators` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `PartialDerivative`.

</details>

<a id="testing"></a>
<details>
<summary>testing</summary>

Créer éventuellement les façades de lancement et les adaptateurs pytest. Ce bloc concerne l’outillage de tests SymPy, pas les calculs utilisateurs. Les helpers internes du runner demandent une adaptation au dépôt ; priorité faible.

**`symox.testing` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `doctest`, `test`.

**`symox.testing.pytest` — 0/16 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `ExceptionInfo`, `Failed`, `RaisesContext`, `SKIP`, `Skipped`, `XFAIL`, `XFail`, `XPass`, `nocache_fail`, `raises`, `skip`, `skip_under_pyodide`, `slow`, `tooslow`, `warns`, `warns_deprecated_sympy`.

**`symox.testing.randtest` — 0/0 actuellement ; `missing_module`.**

- Aucun nom compté ; vérifier le chemin et les objets non mesurés (constantes ou paquet de documentation).

**`symox.testing.runtests` — 0/21 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `DependencyError`, `PyTestReporter`, `Reporter`, `Skipped`, `SymPyDocTestFinder`, `SymPyDocTestRunner`, `SymPyDocTests`, `SymPyOutputChecker`, `SymPyTestResults`, `SymPyTests`, `TimeOutError`, `convert_to_native_paths`, `doctest`, `get_sympy_dir`, `raise_on_deprecated`, `run_all_tests`, `run_in_subprocess_with_hash_randomization`, `setup_pprint`, `split_list`, `sympytestfile`, `test`.

</details>

<a id="utilities"></a>
<details>
<summary>utilities</summary>

Les 25 exports du paquet et les six objets de utilities.enumerative sont raccordés. lambdify est disponible sous son chemin compatible, et as_int dans misc. lambdastr/implemented_function, autowrap/ufuncify et les objets CodeGen/Routine nécessitent encore des adaptateurs ou implémentations. symox.codegen reste un module, jamais utilisé comme substitut d’une fonction codegen.

**`symox.utilities.autowrap` — 0/9 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `CodeWrapError`, `CodeWrapper`, `CythonCodeWrapper`, `DummyWrapper`, `F2PyCodeWrapper`, `UfuncifyCodeWrapper`, `autowrap`, `binary_function`, `ufuncify`.

**`symox.utilities.codegen` — 0/16 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `Argument`, `CCodeGen`, `CodeGen`, `DataType`, `FCodeGen`, `InputArgument`, `JuliaCodeGen`, `OctaveCodeGen`, `OutputArgument`, `Result`, `Routine`, `RustCodeGen`, `codegen`, `default_datatypes`, `get_default_datatype`, `make_routine`.

**`symox.utilities.decorator` — 4/10 actuellement ; `available`.**

- Adaptation / développement à examiner : `T`, `conserve_mpmath_dps`, `deprecated`, `doctest_depends_on`, `no_attrs_in_subclass`, `threaded_factory`.

**`symox.utilities.exceptions` — 0/3 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `SymPyDeprecationWarning`, `ignore_warnings`, `sympy_deprecation_warning`.

**`symox.utilities.lambdify` — 1/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `implemented_function`, `lambdastr`.

**`symox.utilities.memoization` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `assoc_recurrence_memo`, `recurrence_memo`.

**`symox.utilities.misc` — 2/13 actuellement ; `available`.**

- Adaptation / développement à examiner : `Undecidable`, `debug`, `debug_decorator`, `debugf`, `find_executable`, `func_name`, `ordinal`, `rawlines`, `replace`, `strlines`, `translate`.

**`symox.utilities.source` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `get_class`, `get_mod_func`.

**`symox.utilities.timeutils` — 1/2 actuellement ; `available`.**

- Adaptation / développement à examiner : `timethis`.

</details>

<a id="vector"></a>
<details>
<summary>vector</summary>

Aucun raccordement manquant dans la surface inventoriée.

Aucun nom manquant dans les chemins mesurés.

</details>
