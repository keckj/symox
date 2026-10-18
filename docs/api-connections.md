# Implémentation de l’API SymPy, section par section

[Inventaire mesuré](sympy-api.md) · [Comportements et limites des nouveaux modules](api-implementation.md) · [README](../README.md)

Mesure du 2026-10-18 contre SymPy 1.14.0. L’audit initial reste disponible dans le commit `e2a1816`.

311 chemins : 287 importables avec des exports mesurables, 19 manquants, 4 sans exports mesurables et 1 erreur dans la référence. 69 lignes ont encore des noms absents.

Validation de cette passe : 26 689 tests réussis dans la suite complète, puis 2 tests supplémentaires réussis pour les noyaux entiers ; corpus différentiel 273/273. Les manifestes des quatre fractions de la suite ont été rapprochés de la collecte complète.

Les nombres mesurent la présence des noms. Ils ne vérifient ni les signatures ni les comportements. Les lignes se recouvrent : leurs nombres ne s’additionnent pas en une couverture globale. Les chemins sans exports et les erreurs de référence restent visibles.

| Section | Chemins | Chemins absents | Lignes avec noms absents |
|---|---:|---:|---:|
| [root](#root) | 1 | 0 | 0 |
| [abc](#abc) | 1 | 0 | 0 |
| [algebras](#algebras) | 1 | 0 | 0 |
| [assumptions](#assumptions) | 5 | 1 | 4 |
| [calculus](#calculus) | 5 | 0 | 0 |
| [categories](#categories) | 2 | 0 | 0 |
| [codegen](#codegen) | 13 | 0 | 0 |
| [combinatorics](#combinatorics) | 16 | 4 | 11 |
| [core](#core) | 23 | 1 | 18 |
| [crypto](#crypto) | 2 | 0 | 0 |
| [diffgeom](#diffgeom) | 1 | 0 | 0 |
| [discrete](#discrete) | 3 | 0 | 0 |
| [functions](#functions) | 14 | 0 | 3 |
| [geometry](#geometry) | 9 | 0 | 0 |
| [holonomic](#holonomic) | 1 | 0 | 0 |
| [integrals](#integrals) | 5 | 3 | 2 |
| [interactive](#interactive) | 3 | 0 | 1 |
| [liealgebras](#liealgebras) | 1 | 0 | 0 |
| [logic](#logic) | 2 | 0 | 0 |
| [matrices](#matrices) | 8 | 0 | 0 |
| [ntheory](#ntheory) | 14 | 0 | 0 |
| [parsing](#parsing) | 3 | 1 | 1 |
| [physics](#physics) | 68 | 0 | 0 |
| [plotting](#plotting) | 3 | 0 | 1 |
| [polys](#polys) | 17 | 0 | 0 |
| [printing](#printing) | 31 | 4 | 16 |
| [series](#series) | 2 | 0 | 0 |
| [sets](#sets) | 5 | 0 | 0 |
| [simplify](#simplify) | 13 | 1 | 6 |
| [solvers](#solvers) | 9 | 0 | 2 |
| [stats](#stats) | 6 | 0 | 1 |
| [tensor](#tensor) | 7 | 0 | 0 |
| [testing](#testing) | 4 | 4 | 3 |
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

register_handler et remove_handler sont raccordés dans le paquet et ask. Les relations binaires, métaclasses de prédicats, tables de faits et fonctions refine_* nécessitent encore des interfaces supplémentaires.

**`symox.assumptions` — 10/12 ; `available`.**

Noms absents : `AppliedBinaryRelation`, `BinaryRelation`.

**`symox.assumptions.ask` — 4/6 ; `available`.**

Noms absents : `get_all_known_facts`, `get_known_facts_dict`.

**`symox.assumptions.assume` — 4/6 ; `available`.**

Noms absents : `PredicateMeta`, `UndefinedPredicate`.

**`symox.assumptions.predicates` — 0/0 ; `missing_module`.**

Aucun export mesuré ; le chemin reste à ajouter.

**`symox.assumptions.refine` — 1/9 ; `available`.**

Noms absents : `refine_Pow`, `refine_abs`, `refine_arg`, `refine_atan2`, `refine_im`, `refine_matrixelement`, `refine_re`, `refine_sign`.

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

DirectProduct, Coset et SymmetricPermutationGroup sont raccordés aux modules voisins. Les groupes nommés supplémentaires, utilitaires Gray/RGS et algorithmes de canonicalisation tensorielle demandent encore des ajouts. Les alias Perm/PermGroup nécessitent un examen de leur protocole.

**`symox.combinatorics.galois` — 4/21 ; `available`.**

Noms absents : `A4_in_S6`, `A4xC2`, `G18`, `G36m`, `G36p`, `G72`, `M20`, `PGL2F5`, `PSL2F5`, `S3_in_S6`, `S4m`, `S4p`, `S4xC2`, `S5TransitiveSubgroups`, `S6TransitiveSubgroups`, `find_transitive_subgroups_of_S6`, `four_group`.

**`symox.combinatorics.generators` — 4/6 ; `available`.**

Noms absents : `rubik`, `rubik_cube_generators`.

**`symox.combinatorics.graycode` — 1/6 ; `available`.**

Noms absents : `bin_to_gray`, `get_subset_from_bitstring`, `gray_to_bin`, `graycode_subsets`, `random_bitstring`.

**`symox.combinatorics.group_numbers` — 0/4 ; `missing_module`.**

Noms absents : `groups_count`, `is_abelian_number`, `is_cyclic_number`, `is_nilpotent_number`.

**`symox.combinatorics.named_groups` — 6/8 ; `available`.**

Noms absents : `set_alternating_group_properties`, `set_symmetric_group_properties`.

**`symox.combinatorics.partitions` — 5/7 ; `available`.**

Noms absents : `RGS_generalized`, `random_integer_partition`.

**`symox.combinatorics.perm_groups` — 3/4 ; `available`.**

Noms absents : `PermGroup`.

**`symox.combinatorics.permutations` — 2/4 ; `available`.**

Noms absents : `AppliedPermutation`, `Perm`.

**`symox.combinatorics.subsets` — 1/2 ; `available`.**

Noms absents : `ksubsets`.

**`symox.combinatorics.tensor_can` — 0/12 ; `missing_module`.**

Noms absents : `bsgs_direct_product`, `canonical_free`, `canonicalize`, `double_coset_can_rep`, `dummy_sgs`, `gens_products`, `get_minimal_bsgs`, `get_symmetric_group_sgs`, `get_transversals`, `perm_af_direct_product`, `tensor_gens`, `transversal2coset`.

**`symox.combinatorics.testutil` — 0/2 ; `missing_module`.**

Noms absents : `canonicalize_naive`, `graph_certificate`.

**`symox.combinatorics.util` — 0/0 ; `missing_module`.**

Aucun export mesuré ; le chemin reste à ajouter.

</details>

<a id="core"></a>
<details>
<summary>core</summary>

Les 101 exports du paquet sont présents. Les façades add, mul, power, symbol, relational, traversal, cache et expr sont raccordées aux objets existants. assumptions/common_assumptions interrogent les propriétés natives ; evalf reste un module et sympify une fonction, comme dans la référence. Les classes de nombres singleton, métaclasses et routines internes evalf ci-dessous restent à compléter.

**`symox.core.assumptions` — 4/8 ; `available`.**

Noms absents : `ManagedProperties`, `StdFactKB`, `as_property`, `make_property`.

**`symox.core.basic` — 2/3 ; `available`.**

Noms absents : `as_Basic`.

**`symox.core.cache` — 3/5 ; `available`.**

Noms absents : `lazy_function`, `print_cache`.

**`symox.core.containers` — 2/5 ; `available`.**

Noms absents : `OrderedSet`, `TupleKind`, `tuple_wrapper`.

**`symox.core.evalf` — 3/44 ; `available`.**

Noms absents : `EvalfMixin`, `add_terms`, `as_mpmath`, `bitcount`, `check_convergence`, `check_target`, `chop_parts`, `complex_accuracy`, `do_integral`, `evalf_abs`, `evalf_add`, `evalf_alg_num`, `evalf_atan`, `evalf_ceiling`, `evalf_exp`, `evalf_float`, `evalf_floor`, `evalf_im`, `evalf_integer`, `evalf_integral`, `evalf_log`, `evalf_mul`, `evalf_piecewise`, `evalf_pow`, `evalf_prod`, `evalf_rational`, `evalf_re`, `evalf_subs`, `evalf_sum`, `evalf_symbol`, `evalf_trig`, `fastlog`, `finalize_complex`, `get_abs`, `get_complex_part`, `get_integer_part`, `hypsum`, `iszero`, `pure_complex`, `quad_to_mpmath`, `scaled_zero`.

**`symox.core.expr` — 4/6 ; `available`.**

Noms absents : `ExprBuilder`, `unchanged`.

**`symox.core.exprtools` — 3/7 ; `available`.**

Noms absents : `Factors`, `Term`, `decompose_power`, `decompose_power_rat`.

**`symox.core.function` — 22/28 ; `available`.**

Noms absents : `Application`, `AppliedUndef`, `ArgumentIndexError`, `DefinedFunction`, `UndefSageHelper`, `UndefinedFunction`.

**`symox.core.intfunc` — 9/10 ; `available`.**

Noms absents : `igcd_lehmer`.

**`symox.core.kind` — 1/3 ; `available`.**

Noms absents : `KindDispatcher`, `KindMeta`.

**`symox.core.mul` — 2/4 ; `available`.**

Noms absents : `NC_Marker`, `expand_2arg`.

**`symox.core.multidimensional` — 1/4 ; `available`.**

Noms absents : `apply_on_element`, `iter_copy`, `structure_copy`.

**`symox.core.numbers` — 14/34 ; `available`.**

Noms absents : `ComplexInfinity`, `Exp1`, `Half`, `ImaginaryUnit`, `Infinity`, `IntegerConstant`, `NaN`, `NegativeInfinity`, `NegativeOne`, `One`, `Pi`, `RationalConstant`, `Zero`, `all_close`, `equal_valued`, `int_valued`, `mpf_norm`, `sympify_complex`, `sympify_fractions`, `sympify_mpmath`.

**`symox.core.random` — 0/4 ; `missing_module`.**

Noms absents : `random_complex_number`, `seed`, `test_derivative_numerically`, `verify_numerically`.

**`symox.core.singleton` — 1/3 ; `available`.**

Noms absents : `Singleton`, `SingletonRegistry`.

**`symox.core.symbol` — 5/8 ; `available`.**

Noms absents : `Str`, `disambiguate`, `uniquely_named_symbol`.

**`symox.core.sympify` — 2/4 ; `available`.**

Noms absents : `CantSympify`, `kernS`.

**`symox.core.traversal` — 4/7 ; `available`.**

Noms absents : `iterargs`, `iterfreeargs`, `walk`.

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

Les 178 exports du paquet sont raccordés à la racine, sans figer les constructeurs avant leur installation finale. Les façades elliptic_integrals et zeta_functions sont complètes ; error_functions et polynomials gardent les classes de base et auxiliaires manquants ci-dessous.

**`symox.functions.special.beta_functions` — 3/4 ; `available`.**

Noms absents : `betainc_mpmath_fix`.

**`symox.functions.special.error_functions` — 18/21 ; `available`.**

Noms absents : `FresnelIntegral`, `TrigonometricIntegral`, `real_to_real_as_real_imag`.

**`symox.functions.special.polynomials` — 13/14 ; `available`.**

Noms absents : `OrthogonalPolynomial`.

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

Le paquet (30/30) et transforms (32/32) sont raccordés. Les chemins integrals, laplace et singularityfunctions exposent également les fonctions existantes. intpoly et meijerint restent à développer ; leurs homonymes dans d’autres familles ne sont pas des remplacements.

**`symox.integrals.intpoly` — 0/22 ; `missing_module`.**

Noms absents : `best_origin`, `cross_product`, `decompose`, `distance_to_side`, `gradient_terms`, `hyperplane_parameters`, `integration_reduction`, `integration_reduction_dynamic`, `intersection`, `is_vertex`, `left_integral2D`, `left_integral3D`, `lineseg_integrate`, `main_integrate`, `main_integrate3d`, `norm`, `plot_polynomial`, `plot_polytope`, `point_sort`, `polygon_integrate`, `polytope_integrate`, `strip`.

**`symox.integrals.meijerint` — 0/3 ; `missing_module`.**

Noms absents : `meijerint_definite`, `meijerint_indefinite`, `meijerint_inversion`.

**`symox.integrals.meijerint_doc` — 0/0 ; `missing_module`.**

Aucun export mesuré ; le chemin reste à ajouter.

</details>

<a id="interactive"></a>
<details>
<summary>interactive</summary>

interactive_traversal, init_printing et init_session sont raccordés. Les transformations automatiques de saisie et initialisations spécifiques IPython/Python demandent encore des adaptateurs de session.

**`symox.interactive.session` — 1/6 ; `available`.**

Noms absents : `enable_automatic_int_sympification`, `enable_automatic_symbols`, `init_ipython_session`, `init_python_session`, `int_to_Integer`.

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

SymPyExpression manque. Il faut un adaptateur de conversion/parsing de programmes, pas un alias de parse_expr.

**`symox.parsing.sym_expr` — 0/1 ; `missing_module`.**

Noms absents : `SymPyExpression`.

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

La façade pygletplot expose PygletPlot existant avec son backend actuel ; cela ne crée pas un nouveau backend Pyglet. plot_factory et check_arguments nécessitent encore des adaptateurs.

**`symox.plotting.plot` — 10/12 ; `available`.**

Noms absents : `check_arguments`, `plot_factory`.

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

Printers Python/NumPy/mpmath, C89/C99/C11, Fortran et Rust ajoutés avec rendu des nœuds de programme. Les autres protocoles d’affichage encore absents sont listés ci-dessous.

**`symox.printing.aesaracode` — 0/4 ; `missing_module`.**

Noms absents : `AesaraPrinter`, `aesara_code`, `aesara_function`, `dim_handling`.

**`symox.printing.conventions` — 0/2 ; `missing_module`.**

Noms absents : `requires_partial`, `split_super_sub`.

**`symox.printing.jscode` — 3/4 ; `available`.**

Noms absents : `known_functions`.

**`symox.printing.julia` — 2/5 ; `available`.**

Noms absents : `known_fcns_src1`, `known_fcns_src2`, `print_julia_code`.

**`symox.printing.latex` — 4/7 ; `available`.**

Noms absents : `accepted_latex_functions`, `latex_escape`, `translate`.

**`symox.printing.mathematica` — 2/3 ; `available`.**

Noms absents : `known_functions`.

**`symox.printing.mathml` — 5/6 ; `available`.**

Noms absents : `MathMLPrinter`.

**`symox.printing.octave` — 2/5 ; `available`.**

Noms absents : `known_fcns_src1`, `known_fcns_src2`, `print_octave_code`.

**`symox.printing.precedence` — 2/12 ; `available`.**

Noms absents : `PRECEDENCE_FUNCTIONS`, `PRECEDENCE_VALUES`, `precedence_Float`, `precedence_FracElement`, `precedence_Integer`, `precedence_Mul`, `precedence_PolyElement`, `precedence_Rational`, `precedence_UnevaluatedExpr`, `precedence_traditional`.

**`symox.printing.pretty.pretty` — 7/19 ; `available`.**

Noms absents : `U`, `annotated`, `center_pad`, `hobj`, `is_subscriptable_in_unicode`, `pretty_atom`, `pretty_symbol`, `pretty_try_use_unicode`, `pretty_use_unicode`, `vobj`, `xobj`, `xsym`.

**`symox.printing.pretty.pretty_symbology` — 0/25 ; `missing_module`.**

Noms absents : `G`, `U`, `VF`, `annotated`, `atoms_table`, `center`, `center_pad`, `digit_2txt`, `frac`, `g`, `greek_letters`, `greek_unicode`, `hobj`, `pretty_atom`, `pretty_symbol`, `pretty_try_use_unicode`, `pretty_use_unicode`, `root`, `sub`, `sup`, `symb_2txt`, `vobj`, `xobj`, `xstr`, `xsym`.

**`symox.printing.pretty.stringpict` — 0/2 ; `missing_module`.**

Noms absents : `prettyForm`, `stringPict`.

**`symox.printing.printer` — 1/3 ; `available`.**

Noms absents : `print_function`, `printer_context`.

**`symox.printing.python` — 2/3 ; `available`.**

Noms absents : `PythonPrinter`.

**`symox.printing.rcode` — 3/4 ; `available`.**

Noms absents : `known_functions`.

**`symox.printing.str` — 3/4 ; `available`.**

Noms absents : `StrReprPrinter`.

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

Les 32 exports du paquet et les 31 de fu, règles TR* comprises, sont raccordés. Les façades réutilisent _simptools, _epath et les fonctions natives ; elles suivent les adaptateurs installés ultérieurement. Les internes CSE, certaines réécritures et collections/formules hypergéométriques restent à développer. Aucun opérateur quantique homonyme n’a été réexporté.

**`symox.simplify.cse_main` — 1/11 ; `available`.**

Noms absents : `FuncArgTracker`, `Unevaluated`, `cse_release_variables`, `cse_separate`, `match_common_args`, `opt_cse`, `postprocess_for_cse`, `preprocess_for_cse`, `reps_toposort`, `tree_cse`.

**`symox.simplify.hyperexpand` — 1/37 ; `available`.**

Noms absents : `Formula`, `FormulaCollection`, `G_Function`, `Hyper_Function`, `MeijerFormula`, `MeijerFormulaCollection`, `MeijerShiftA`, `MeijerShiftB`, `MeijerShiftC`, `MeijerShiftD`, `MeijerUnShiftA`, `MeijerUnShiftB`, `MeijerUnShiftC`, `MeijerUnShiftD`, `MultOperator`, `Operator`, `ReduceOrder`, `ShiftA`, `ShiftB`, `UnShiftA`, `UnShiftB`, `add_formulae`, `add_meijerg_formulae`, `apply_operators`, `build_hypergeometric_formula`, `debug`, `devise_plan`, `devise_plan_meijer`, `hyperexpand_special`, `make_derivative_operator`, `make_simp`, `reduce_order`, `reduce_order_meijer`, `try_lerchphi`, `try_polynomial`, `try_shifted_sum`.

**`symox.simplify.hyperexpand_doc` — 0/0 ; `missing_module`.**

Aucun export mesuré ; le chemin reste à ajouter.

**`symox.simplify.radsimp` — 7/17 ; `available`.**

Noms absents : `collect_abs`, `collect_sqrt`, `denom_expand`, `expand_denom`, `expand_fraction`, `expand_numer`, `fraction_expand`, `numer_expand`, `rad_rationalize`, `split_surds`.

**`symox.simplify.simplify` — 10/21 ; `available`.**

Noms absents : `clear_coefficients`, `dotprodsimp`, `factor_sum`, `inversecombine`, `nc_simplify`, `nthroot`, `product_mul`, `product_simplify`, `sum_add`, `sum_combine`, `sum_simplify`.

**`symox.simplify.sqrtdenest` — 1/6 ; `available`.**

Noms absents : `SqrtdenestStopIteration`, `is_algebraic`, `is_sqrt`, `sqrt_biquadratic_denest`, `sqrt_depth`.

**`symox.simplify.trigsimp` — 3/5 ; `available`.**

Noms absents : `trigsimp_groebner`, `trigsimp_old`.

</details>

<a id="solvers"></a>
<details>
<summary>solvers</summary>

Solveur linéaire exact à deux phases, avec contraintes d’égalité, bornes et variables libres, ajouté dans simplex. Les interfaces internes ODE et solveset encore absentes restent listées.

**`symox.solvers.ode.ode` — 5/23 ; `available`.**

Noms absents : `check_linear_2eq_order1`, `check_nonlinear_2eq_order1`, `check_nonlinear_2eq_order2`, `check_nonlinear_3eq_order1`, `check_nonlinear_3eq_order2`, `classify_sysode`, `constant_renumber`, `get_numbered_constants`, `iter_numbered_constants`, `ode_1st_power_series`, `ode_2nd_power_series_ordinary`, `ode_2nd_power_series_regular`, `ode_sol_simplicity`, `odesimp`, `solve_ics`, `sysode_linear_2eq_order1`, `sysode_nonlinear_2eq_order1`, `sysode_nonlinear_3eq_order1`.

**`symox.solvers.solveset` — 8/14 ; `available`.**

Noms absents : `domain_check`, `invert_complex`, `invert_real`, `linear_coeffs`, `solve_decomposition`, `solvify`.

</details>

<a id="stats"></a>
<details>
<summary>stats</summary>

Domaines et espaces de probabilité finis/continus ajoutés, reliés à pspace et aux moteurs d’espérance/probabilité. Les lois jointes finies gardent leur dépendance pendant les calculs et l’échantillonnage. Les limites de Statistics restent applicables.

**`symox.stats.rv` — 39/41 ; `available`.**

Noms absents : `RandomIndexedSymbol`, `RandomMatrixSymbol`.

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

Créer éventuellement les façades de lancement et les adaptateurs pytest. Ce bloc concerne l’outillage de tests SymPy, pas les calculs utilisateurs. Les helpers internes du runner demandent une adaptation au dépôt ; priorité faible.

**`symox.testing` — 0/2 ; `missing_module`.**

Noms absents : `doctest`, `test`.

**`symox.testing.pytest` — 0/16 ; `missing_module`.**

Noms absents : `ExceptionInfo`, `Failed`, `RaisesContext`, `SKIP`, `Skipped`, `XFAIL`, `XFail`, `XPass`, `nocache_fail`, `raises`, `skip`, `skip_under_pyodide`, `slow`, `tooslow`, `warns`, `warns_deprecated_sympy`.

**`symox.testing.randtest` — 0/0 ; `missing_module`.**

Aucun export mesuré ; le chemin reste à ajouter.

**`symox.testing.runtests` — 0/21 ; `missing_module`.**

Noms absents : `DependencyError`, `PyTestReporter`, `Reporter`, `Skipped`, `SymPyDocTestFinder`, `SymPyDocTestRunner`, `SymPyDocTests`, `SymPyOutputChecker`, `SymPyTestResults`, `SymPyTests`, `TimeOutError`, `convert_to_native_paths`, `doctest`, `get_sympy_dir`, `raise_on_deprecated`, `run_all_tests`, `run_in_subprocess_with_hash_randomization`, `setup_pprint`, `split_list`, `sympytestfile`, `test`.

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
