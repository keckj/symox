# Raccordements de l’API SymPy, section par section

[Inventaire mesuré](sympy-api.md) · [README](../README.md)

Audit du 2026-10-18, SymPy 1.14.0, après le commit Physics `a845521`. Aucune API modifiée par cet audit.

311 chemins examinés : 204 importables avec des exports mesurables, 103 chemins manquants, 3 modules sans exports mesurables et 1 erreur d’import de référence. 163 lignes ont des noms manquants. Les lignes se recouvrent : ne pas additionner leurs noms pour calculer une couverture globale.

Un **raccordement** désigne un objet déjà exposé par le moteur ou un module local, à réexporter sous le chemin compatible. Les sources ci-dessous sont des objets locaux retrouvés, sélectionnés par famille ; leur existence ne prouve pas la compatibilité des signatures ou du comportement. La rubrique « adaptation / développement » inclut les objets sans source directe identifiée : certains peuvent être construits à partir de méthodes existantes. Les homonymes sans rapport sont exclus.

Ordre proposé : `functions` → `integrals` → `core` → `simplify` → `polys` / `matrices` / `sets` → `solvers` / `tensor` / `utilities` → compléments et nouveaux algorithmes.

Pour chaque façade, vérifier l’identité des objets réexportés et les deux ordres d’import. Les noms `simplify`, `series`, `lambdify`, `ask` ou `refine` peuvent désigner à la fois une fonction et un module : un nouvel import ne doit pas remplacer accidentellement une fonction publique par un module. Ne pas importer tous les sous-modules depuis la racine pour faire monter le compteur.

| Section | Chemins | Chemins absents | Lignes avec noms manquants |
|---|---:|---:|---:|
| [root](#root) | 1 | 0 | 0 |
| [abc](#abc) | 1 | 1 | 0 |
| [algebras](#algebras) | 1 | 0 | 0 |
| [assumptions](#assumptions) | 5 | 1 | 4 |
| [calculus](#calculus) | 5 | 0 | 1 |
| [categories](#categories) | 2 | 0 | 0 |
| [codegen](#codegen) | 13 | 11 | 13 |
| [combinatorics](#combinatorics) | 16 | 5 | 12 |
| [core](#core) | 23 | 16 | 23 |
| [crypto](#crypto) | 2 | 0 | 0 |
| [diffgeom](#diffgeom) | 1 | 0 | 0 |
| [discrete](#discrete) | 3 | 0 | 0 |
| [functions](#functions) | 14 | 4 | 6 |
| [geometry](#geometry) | 9 | 0 | 4 |
| [holonomic](#holonomic) | 1 | 0 | 0 |
| [integrals](#integrals) | 5 | 5 | 4 |
| [interactive](#interactive) | 3 | 2 | 3 |
| [liealgebras](#liealgebras) | 1 | 0 | 0 |
| [logic](#logic) | 2 | 2 | 2 |
| [matrices](#matrices) | 8 | 5 | 8 |
| [ntheory](#ntheory) | 14 | 0 | 1 |
| [parsing](#parsing) | 3 | 1 | 1 |
| [physics](#physics) | 68 | 0 | 0 |
| [plotting](#plotting) | 3 | 1 | 2 |
| [polys](#polys) | 17 | 7 | 11 |
| [printing](#printing) | 31 | 6 | 24 |
| [series](#series) | 2 | 0 | 0 |
| [sets](#sets) | 5 | 4 | 5 |
| [simplify](#simplify) | 13 | 13 | 12 |
| [solvers](#solvers) | 9 | 2 | 5 |
| [stats](#stats) | 6 | 2 | 3 |
| [tensor](#tensor) | 7 | 4 | 5 |
| [testing](#testing) | 4 | 4 | 3 |
| [utilities](#utilities) | 12 | 7 | 11 |
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

Créer le module de symboles préconstruits à partir de Symbol et des lettres grecques/latines. Le compteur 0/0 ne mesure pas ces constantes : le chemin manque réellement.

**`symox.abc` — 0/0 actuellement ; `missing_module`.**

- Aucun nom compté ; vérifier le chemin et les objets non mesurés (constantes ou paquet de documentation).

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

Réexporter register_handler et remove_handler depuis _public_tools. Les relations binaires, métaclasses de prédicats et fonctions refine_* demandent des interfaces supplémentaires ; le moteur ask/refine existe déjà.

**`symox.assumptions` — 8/12 actuellement ; `available`.**

- Source retrouvée dans `symox` : `register_handler`, `remove_handler`.
- Adaptation / développement à examiner : `AppliedBinaryRelation`, `BinaryRelation`.

**`symox.assumptions.ask` — 2/6 actuellement ; `available`.**

- Source retrouvée dans `symox` : `register_handler`, `remove_handler`.
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

Exporter AccumBounds et euler_equations au niveau du paquet. Les sous-modules mesurés sont déjà raccordés.

**`symox.calculus` — 15/17 actuellement ; `available`.**

- Source retrouvée dans `symox` : `AccumBounds`, `euler_equations`.

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

DirectProduct, Coset et SymmetricPermutationGroup existent dans les modules voisins. Les groupes nommés supplémentaires, utilitaires Gray/RGS et algorithmes de canonicalisation tensorielle demandent des ajouts. Les alias Perm/PermGroup nécessitent un examen de leur protocole.

**`symox.combinatorics.galois` — 4/21 actuellement ; `available`.**

- Adaptation / développement à examiner : `A4_in_S6`, `A4xC2`, `G18`, `G36m`, `G36p`, `G72`, `M20`, `PGL2F5`, `PSL2F5`, `S3_in_S6`, `S4m`, `S4p`, `S4xC2`, `S5TransitiveSubgroups`, `S6TransitiveSubgroups`, `find_transitive_subgroups_of_S6`, `four_group`.

**`symox.combinatorics.generators` — 4/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `rubik`, `rubik_cube_generators`.

**`symox.combinatorics.graycode` — 1/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `bin_to_gray`, `get_subset_from_bitstring`, `gray_to_bin`, `graycode_subsets`, `random_bitstring`.

**`symox.combinatorics.group_constructs` — 0/1 actuellement ; `missing_module`.**

- Source retrouvée dans `symox.combinatorics.named_groups` : `DirectProduct`.

**`symox.combinatorics.group_numbers` — 0/4 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `groups_count`, `is_abelian_number`, `is_cyclic_number`, `is_nilpotent_number`.

**`symox.combinatorics.named_groups` — 6/8 actuellement ; `available`.**

- Adaptation / développement à examiner : `set_alternating_group_properties`, `set_symmetric_group_properties`.

**`symox.combinatorics.partitions` — 5/7 actuellement ; `available`.**

- Adaptation / développement à examiner : `RGS_generalized`, `random_integer_partition`.

**`symox.combinatorics.perm_groups` — 1/4 actuellement ; `available`.**

- Source retrouvée dans `symox.combinatorics.coset` : `Coset`.
- Source retrouvée dans `symox.combinatorics.named_groups` : `SymmetricPermutationGroup`.
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

Priorité élevée : créer les façades add, mul, power, symbol, relational, traversal, cache, expr et compléter core.__init__. Les 70 noms déjà présents à la racine sur les 71 manquants du paquet sont des candidats de réexport. Les classes de nombres singleton, métaclasses et routines internes evalf ne se réduisent pas à des alias des constantes.

**`symox.core` — 30/101 actuellement ; `available`.**

- Source retrouvée dans `symox` : `AlgebraicNumber`, `Atom`, `AtomicExpr`, `BooleanKind`, `Catalan`, `E`, `Eq`, `Equality`, `EulerGamma`, `Ge`, `GoldenRatio`, `GreaterThan`, `Gt`, `I`, `Le`, `LessThan`, `Lt`, `Mod`, `N`, `Ne`, `NumberKind`, `PoleError`, `PrecisionExhausted`, `Rel`, `S`, `StrictGreaterThan`, `StrictLessThan`, `TribonacciConstant`, `UndefinedKind`, `Unequality`, `UnevaluatedExpr`, `bottom_up`, `cacheit`, `check_assumptions`, `count_ops`, `default_sort_key`, `diff`, `evalf`, `expand_complex`, `expand_func`, `expand_log`, `expand_mul`, `expand_power_base`, `expand_power_exp`, `expand_trig`, `factor_nc`, `factor_terms`, `failing_assumptions`, `gcd_terms`, `igcd`, `ilcm`, `integer_log`, `integer_nthroot`, `mod_inverse`, `nan`, `nfloat`, `num_digits`, `oo`, `ordered`, `pi`, `postorder_traversal`, `preorder_traversal`, `prod`, `symbols`, `trailing`, `use`, `var`, `vectorize`, `zoo`.
- Adaptation / développement à examiner : `assumptions`, `common_assumptions`.

**`symox.core.add` — 0/1 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Add`.

**`symox.core.assumptions` — 0/8 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `check_assumptions`, `failing_assumptions`.
- Adaptation / développement à examiner : `ManagedProperties`, `StdFactKB`, `as_property`, `assumptions`, `common_assumptions`, `make_property`.

**`symox.core.basic` — 0/3 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Atom`, `Basic`.
- Adaptation / développement à examiner : `as_Basic`.

**`symox.core.cache` — 0/5 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `cacheit`.
- Source retrouvée dans `symox.utilities.cache` : `clear_cache`.
- Adaptation / développement à examiner : `cached_property`, `lazy_function`, `print_cache`.

**`symox.core.containers` — 2/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `OrderedSet`, `TupleKind`, `tuple_wrapper`.

**`symox.core.evalf` — 3/44 actuellement ; `available`.**

- Adaptation / développement à examiner : `EvalfMixin`, `add_terms`, `as_mpmath`, `bitcount`, `check_convergence`, `check_target`, `chop_parts`, `complex_accuracy`, `do_integral`, `evalf_abs`, `evalf_add`, `evalf_alg_num`, `evalf_atan`, `evalf_ceiling`, `evalf_exp`, `evalf_float`, `evalf_floor`, `evalf_im`, `evalf_integer`, `evalf_integral`, `evalf_log`, `evalf_mul`, `evalf_piecewise`, `evalf_pow`, `evalf_prod`, `evalf_rational`, `evalf_re`, `evalf_subs`, `evalf_sum`, `evalf_symbol`, `evalf_trig`, `fastlog`, `finalize_complex`, `get_abs`, `get_complex_part`, `get_integer_part`, `hypsum`, `iszero`, `pure_complex`, `quad_to_mpmath`, `scaled_zero`.

**`symox.core.expr` — 0/6 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `AtomicExpr`, `Expr`, `UnevaluatedExpr`, `factor_terms`.
- Adaptation / développement à examiner : `ExprBuilder`, `unchanged`.

**`symox.core.exprtools` — 0/7 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `factor_nc`, `factor_terms`, `gcd_terms`.
- Adaptation / développement à examiner : `Factors`, `Term`, `decompose_power`, `decompose_power_rat`.

**`symox.core.function` — 8/28 actuellement ; `available`.**

- Source retrouvée dans `symox` : `BadArgumentsError`, `BadSignatureError`, `PoleError`, `count_ops`, `diff`, `expand`, `expand_complex`, `expand_func`, `expand_log`, `expand_mul`, `expand_power_base`, `expand_power_exp`, `expand_trig`, `nfloat`.
- Adaptation / développement à examiner : `Application`, `AppliedUndef`, `ArgumentIndexError`, `DefinedFunction`, `UndefSageHelper`, `UndefinedFunction`.

**`symox.core.intfunc` — 0/10 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `igcd`, `ilcm`, `integer_log`, `integer_nthroot`, `mod_inverse`, `num_digits`, `trailing`.
- Source retrouvée dans `symox.ntheory._intfunc` : `igcdex`, `isqrt`.
- Adaptation / développement à examiner : `igcd_lehmer`.

**`symox.core.kind` — 0/3 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Kind`.
- Adaptation / développement à examiner : `KindDispatcher`, `KindMeta`.

**`symox.core.mod` — 0/1 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Mod`.

**`symox.core.mul` — 0/4 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Mul`, `prod`.
- Adaptation / développement à examiner : `NC_Marker`, `expand_2arg`.

**`symox.core.multidimensional` — 1/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `apply_on_element`, `iter_copy`, `structure_copy`.

**`symox.core.numbers` — 8/34 actuellement ; `available`.**

- Source retrouvée dans `symox` : `AlgebraicNumber`, `Catalan`, `EulerGamma`, `GoldenRatio`, `TribonacciConstant`, `mod_inverse`.
- Adaptation / développement à examiner : `ComplexInfinity`, `Exp1`, `Half`, `ImaginaryUnit`, `Infinity`, `IntegerConstant`, `NaN`, `NegativeInfinity`, `NegativeOne`, `One`, `Pi`, `RationalConstant`, `Zero`, `all_close`, `equal_valued`, `int_valued`, `mpf_norm`, `sympify_complex`, `sympify_fractions`, `sympify_mpmath`.

**`symox.core.power` — 0/1 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Pow`.

**`symox.core.random` — 0/4 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `random_complex_number`, `seed`, `test_derivative_numerically`, `verify_numerically`.

**`symox.core.relational` — 0/14 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Eq`, `Equality`, `Ge`, `GreaterThan`, `Gt`, `Le`, `LessThan`, `Lt`, `Ne`, `Rel`, `Relational`, `StrictGreaterThan`, `StrictLessThan`, `Unequality`.

**`symox.core.singleton` — 0/3 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `S`.
- Adaptation / développement à examiner : `Singleton`, `SingletonRegistry`.

**`symox.core.symbol` — 0/8 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Dummy`, `Symbol`, `Wild`, `symbols`, `var`.
- Adaptation / développement à examiner : `Str`, `disambiguate`, `uniquely_named_symbol`.

**`symox.core.sympify` — 2/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `CantSympify`, `kernS`.

**`symox.core.traversal` — 0/7 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `bottom_up`, `postorder_traversal`, `preorder_traversal`, `use`.
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

Priorité élevée : les 178 noms manquants de functions existent tous à la racine. Créer les façades elliptic_integrals, error_functions, polynomials et zeta_functions. Les classes de base FresnelIntegral, TrigonometricIntegral et OrthogonalPolynomial ainsi que deux auxiliaires restent à adapter ou développer.

**`symox.functions` — 0/178 actuellement ; `available`.**

- Source retrouvée dans `symox` : `Abs`, `Chi`, `Ci`, `DiracDelta`, `E1`, `Ei`, `Eijk`, `FallingFactorial`, `Heaviside`, `Id`, `KroneckerDelta`, `LambertW`, `LeviCivita`, `Li`, `Max`, `Min`, `Piecewise`, `Rem`, `RisingFactorial`, `Shi`, `Si`, `SingularityFunction`, `Ynm`, `Ynm_c`, `Znm`, `acos`, `acosh`, `acot`, `acoth`, `acsc`, `acsch`, `adjoint`, `airyai`, `airyaiprime`, `airybi`, `airybiprime`, `andre`, `appellf1`, `arg`, `asec`, `asech`, `asin`, `asinh`, `assoc_laguerre`, `assoc_legendre`, `atan`, `atan2`, `atanh`, `bell`, `bernoulli`, `besseli`, `besselj`, `besselk`, `bessely`, `beta`, `betainc`, `betainc_regularized`, `binomial`, `bspline_basis`, `bspline_basis_set`, `carmichael`, `catalan`, `cbrt`, `ceiling`, `chebyshevt`, `chebyshevt_root`, `chebyshevu`, `chebyshevu_root`, `conjugate`, `cos`, `cosh`, `cot`, `coth`, `csc`, `csch`, `digamma`, `dirichlet_eta`, `divisor_sigma`, `elliptic_e`, `elliptic_f`, `elliptic_k`, `elliptic_pi`, `erf`, `erf2`, `erf2inv`, `erfc`, `erfcinv`, `erfi`, `erfinv`, `euler`, `exp`, `exp_polar`, `expint`, `factorial`, `factorial2`, `ff`, `fibonacci`, `floor`, `frac`, `fresnelc`, `fresnels`, `gamma`, `gegenbauer`, `genocchi`, `hankel1`, `hankel2`, `harmonic`, `hermite`, `hermite_prob`, `hn1`, `hn2`, `hyper`, `im`, `interpolating_spline`, `jacobi`, `jacobi_normalized`, `jacobi_symbol`, `jn`, `jn_zeros`, `kronecker_symbol`, `laguerre`, `legendre`, `legendre_symbol`, `lerchphi`, `li`, `ln`, `log`, `loggamma`, `lowergamma`, `lucas`, `marcumq`, `mathieuc`, `mathieucprime`, `mathieus`, `mathieusprime`, `meijerg`, `mobius`, `motzkin`, `multigamma`, `partition`, `periodic_argument`, `piecewise_exclusive`, `piecewise_fold`, `polar_lift`, `polarify`, `polygamma`, `polylog`, `primenu`, `primeomega`, `primepi`, `principal_branch`, `re`, `real_root`, `reduced_totient`, `rf`, `riemann_xi`, `root`, `sec`, `sech`, `sign`, `sin`, `sinc`, `sinh`, `sqrt`, `stieltjes`, `subfactorial`, `tan`, `tanh`, `totient`, `transpose`, `tribonacci`, `trigamma`, `udivisor_sigma`, `unbranched_argument`, `unpolarify`, `uppergamma`, `yn`, `zeta`.

**`symox.functions.special.beta_functions` — 3/4 actuellement ; `available`.**

- Adaptation / développement à examiner : `betainc_mpmath_fix`.

**`symox.functions.special.elliptic_integrals` — 0/4 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `elliptic_e`, `elliptic_f`, `elliptic_k`, `elliptic_pi`.

**`symox.functions.special.error_functions` — 0/21 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Chi`, `Ci`, `E1`, `Ei`, `Li`, `Shi`, `Si`, `erf`, `erf2`, `erf2inv`, `erfc`, `erfcinv`, `erfi`, `erfinv`, `expint`, `fresnelc`, `fresnels`, `li`.
- Adaptation / développement à examiner : `FresnelIntegral`, `TrigonometricIntegral`, `real_to_real_as_real_imag`.

**`symox.functions.special.polynomials` — 0/14 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `assoc_laguerre`, `assoc_legendre`, `chebyshevt`, `chebyshevt_root`, `chebyshevu`, `chebyshevu_root`, `gegenbauer`, `hermite`, `hermite_prob`, `jacobi`, `jacobi_normalized`, `laguerre`, `legendre`.
- Adaptation / développement à examiner : `OrthogonalPolynomial`.

**`symox.functions.special.zeta_functions` — 0/6 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `dirichlet_eta`, `lerchphi`, `polylog`, `riemann_xi`, `stieltjes`, `zeta`.

</details>

<a id="geometry"></a>
<details>
<summary>geometry</summary>

Réexporter deg/rad dans polygon. are_coplanar existe comme méthode de Point et peut servir à un adaptateur. GeometrySet, LinearEntity2D/3D et les fonctions libres de transformation demandent un examen de hiérarchie et de signatures. Ne pas réutiliser find des expressions pour find de geometry.util.

**`symox.geometry.entity` — 1/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `GeometrySet`, `rotate`, `scale`, `translate`.

**`symox.geometry.line` — 10/12 actuellement ; `available`.**

- Adaptation / développement à examiner : `LinearEntity2D`, `LinearEntity3D`.

**`symox.geometry.polygon` — 3/5 actuellement ; `available`.**

- Source retrouvée dans `symox.geometry.util` : `deg`, `rad`.

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

Priorité élevée : créer le paquet et sa façade transforms. Les 30 exports du paquet existent à la racine ; les 32 de transforms existent dans _transforms. intpoly et meijerint demandent de nouveaux algorithmes/interfaces. decompose des polynômes et intersection de géométrie ne remplacent pas leurs homonymes dans intpoly.

**`symox.integrals` — 0/30 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `CosineTransform`, `FourierTransform`, `HankelTransform`, `Integral`, `InverseCosineTransform`, `InverseFourierTransform`, `InverseHankelTransform`, `InverseLaplaceTransform`, `InverseMellinTransform`, `InverseSineTransform`, `LaplaceTransform`, `MellinTransform`, `SineTransform`, `cosine_transform`, `fourier_transform`, `hankel_transform`, `integrate`, `inverse_cosine_transform`, `inverse_fourier_transform`, `inverse_hankel_transform`, `inverse_laplace_transform`, `inverse_mellin_transform`, `inverse_sine_transform`, `laplace_correspondence`, `laplace_initial_conds`, `laplace_transform`, `line_integrate`, `mellin_transform`, `sine_transform`, `singularityintegrate`.

**`symox.integrals.intpoly` — 0/22 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `best_origin`, `cross_product`, `decompose`, `distance_to_side`, `gradient_terms`, `hyperplane_parameters`, `integration_reduction`, `integration_reduction_dynamic`, `intersection`, `is_vertex`, `left_integral2D`, `left_integral3D`, `lineseg_integrate`, `main_integrate`, `main_integrate3d`, `norm`, `plot_polynomial`, `plot_polytope`, `point_sort`, `polygon_integrate`, `polytope_integrate`, `strip`.

**`symox.integrals.meijerint` — 0/3 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `meijerint_definite`, `meijerint_indefinite`, `meijerint_inversion`.

**`symox.integrals.meijerint_doc` — 0/0 actuellement ; `missing_module`.**

- Aucun nom compté ; vérifier le chemin et les objets non mesurés (constantes ou paquet de documentation).

**`symox.integrals.transforms` — 0/32 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `CosineTransform`, `FourierTransform`, `HankelTransform`, `IntegralTransform`, `IntegralTransformError`, `InverseCosineTransform`, `InverseFourierTransform`, `InverseHankelTransform`, `InverseLaplaceTransform`, `InverseMellinTransform`, `InverseSineTransform`, `LaplaceTransform`, `MellinTransform`, `SineTransform`, `cosine_transform`, `fourier_transform`, `hankel_transform`, `inverse_cosine_transform`, `inverse_fourier_transform`, `inverse_hankel_transform`, `inverse_laplace_transform`, `inverse_mellin_transform`, `inverse_sine_transform`, `laplace_correspondence`, `laplace_initial_conds`, `laplace_transform`, `mellin_transform`, `sine_transform`.
- Source retrouvée dans `symox._transforms` : `FourierTypeTransform`, `HankelTypeTransform`, `MellinTransformStripError`, `SineCosineTypeTransform`.

</details>

<a id="interactive"></a>
<details>
<summary>interactive</summary>

Raccorder interactive_traversal, init_printing et init_session aux chemins compatibles. Les transformations automatiques de saisie et initialisations spécifiques IPython/Python demandent des adaptateurs de session.

**`symox.interactive` — 2/3 actuellement ; `available`.**

- Source retrouvée dans `symox` : `interactive_traversal`.

**`symox.interactive.printing` — 0/1 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `init_printing`.

**`symox.interactive.session` — 0/6 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `init_session`.
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

Créer logic et logic.inference : opérateurs booléens, formes normales et satisfiable existent à la racine. gateinputcount et les interfaces de base de connaissances KB/PropKB demandent des ajouts.

**`symox.logic` — 0/20 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `And`, `Equivalent`, `ITE`, `Implies`, `Nand`, `Nor`, `Not`, `Or`, `POSform`, `SOPform`, `Xor`, `bool_map`, `false`, `satisfiable`, `simplify_logic`, `to_cnf`, `to_dnf`, `to_nnf`, `true`.
- Adaptation / développement à examiner : `gateinputcount`.

**`symox.logic.inference` — 0/7 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `satisfiable`.
- Adaptation / développement à examiner : `KB`, `PropKB`, `entails`, `literal_symbol`, `pl_true`, `valid`.

</details>

<a id="matrices"></a>
<details>
<summary>matrices</summary>

Priorité élevée : compléter le paquet et les façades immutable, sparse, matrixbase et sparsetools. 28 des 31 exports manquants du paquet existent à la racine. MatrixKind/MatrixSet, CompanionMatrix, les classes intermédiaires et les réécritures blockmatrix nécessitent des adaptations ou implémentations.

**`symox.matrices` — 42/73 actuellement ; `available`.**

- Source retrouvée dans `symox` : `DeferredVector`, `GramSchmidt`, `MatrixBase`, `MutableMatrix`, `NonSquareMatrixError`, `Permanent`, `ShapeError`, `banded`, `blockcut`, `casoratian`, `det`, `diagonalize_vector`, `jordan_cell`, `list2numpy`, `matrix2numpy`, `matrix_multiply_elementwise`, `matrix_symbols`, `per`, `randMatrix`, `rot_axis1`, `rot_axis2`, `rot_axis3`, `rot_ccw_axis1`, `rot_ccw_axis2`, `rot_ccw_axis3`, `rot_givens`, `symarray`, `trace`.
- Adaptation / développement à examiner : `MatrixKind`, `MatrixSet`, `dotprodsimp`.

**`symox.matrices.expressions` — 30/40 actuellement ; `available`.**

- Source retrouvée dans `symox` : `Permanent`, `blockcut`, `det`, `diagonalize_vector`, `matrix_symbols`, `per`, `trace`.
- Adaptation / développement à examiner : `CompanionMatrix`, `MatrixSet`, `combine_kronecker`.

**`symox.matrices.expressions.blockmatrix` — 3/16 actuellement ; `available`.**

- Source retrouvée dans `symox` : `blockcut`.
- Adaptation / développement à examiner : `bc_block_plus_ident`, `bc_dist`, `bc_inverse`, `bc_matadd`, `bc_matmul`, `bc_transpose`, `bc_unpack`, `blockinverse_1x1`, `blockinverse_2x2`, `bounds`, `deblock`, `reblock_2x2`.

**`symox.matrices.immutable` — 0/6 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `ImmutableDenseMatrix`, `ImmutableMatrix`, `ImmutableSparseMatrix`.
- Adaptation / développement à examiner : `ImmutableRepMatrix`, `sympify_matrix`, `sympify_mpmath_matrix`.

**`symox.matrices.kind` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `MatrixKind`.

**`symox.matrices.matrixbase` — 0/7 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `DeferredVector`, `MatrixBase`, `MatrixError`, `NonSquareMatrixError`, `ShapeError`.
- Adaptation / développement à examiner : `a2idx`, `classof`.

**`symox.matrices.sparse` — 0/3 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `MutableSparseMatrix`, `SparseMatrix`.
- Adaptation / développement à examiner : `SparseRepMatrix`.

**`symox.matrices.sparsetools` — 0/1 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `banded`.

</details>

<a id="ntheory"></a>
<details>
<summary>ntheory</summary>

Raccorder digits de ntheory.digits dans ntheory.factor_. Le reste des chemins mesurés est complet en noms.

**`symox.ntheory.factor_` — 37/38 actuellement ; `available`.**

- Source retrouvée dans `symox.ntheory.digits` : `digits`.

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

Créer la façade pygletplot vers PygletPlot existant, en conservant le backend réellement proposé. plot_factory et check_arguments nécessitent des adaptateurs ; vérifier le contrat avant de leur donner le nom d’une fonction existante.

**`symox.plotting.plot` — 10/12 actuellement ; `available`.**

- Adaptation / développement à examiner : `check_arguments`, `plot_factory`.

**`symox.plotting.pygletplot` — 0/1 actuellement ; `missing_module`.**

- Source retrouvée dans `symox.plotting.plot` : `PygletPlot`.

</details>

<a id="polys"></a>
<details>
<summary>polys</summary>

Priorité élevée : 81 des 86 exports manquants du paquet existent à la racine. Raccorder aussi DM à domainmatrix. Les cinq constructeurs de polynômes spéciaux, DDM/SDM/DFM, routines ddm_*/sdm_*, endomorphismes de modules et configuration des algorithmes sont un travail distinct du moteur DomainMatrix existant.

**`symox.polys` — 90/176 actuellement ; `available`.**

- Source retrouvée dans `symox` : `AlgebraicField`, `CC`, `CRootOf`, `ComplexField`, `ComplexRootOf`, `Domain`, `EX`, `EXRAW`, `ExpressionDomain`, `FF`, `FF_gmpy`, `FF_python`, `FiniteField`, `FractionField`, `GF`, `GMPYFiniteField`, `GMPYIntegerRing`, `GMPYRationalField`, `IntegerRing`, `Options`, `PolynomialRing`, `PurePoly`, `PythonFiniteField`, `PythonIntegerRing`, `PythonRational`, `QQ`, `QQ_I`, `QQ_gmpy`, `QQ_python`, `RR`, `RationalField`, `RealField`, `RootOf`, `RootSum`, `ZZ`, `ZZ_I`, `ZZ_gmpy`, `ZZ_python`, `apart`, `cancel`, `cofactors`, `compose`, `degree`, `div`, `exquo`, `field`, `field_isomorphism`, `galois_group`, `gcd_list`, `gcdex`, `half_gcdex`, `invert`, `is_zero_dimensional`, `isolate`, `lcm_list`, `nroots`, `pdiv`, `pexquo`, `poly`, `pquo`, `prem`, `prime_decomp`, `prime_valuation`, `primitive`, `quo`, `rem`, `ring`, `rootof`, `roots`, `round_two`, `sfield`, `sqf`, `sqf_part`, `sring`, `terms_gcd`, `to_number_field`, `together`, `vfield`, `vring`, `xfield`, `xring`.
- Adaptation / développement à examiner : `andre_poly`, `bernoulli_c_poly`, `bernoulli_poly`, `euler_poly`, `genocchi_poly`.

**`symox.polys.matrices._dfm` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `DFM`.

**`symox.polys.matrices._typing` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `RingElement`.

**`symox.polys.matrices.ddm` — 0/1 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `DDM`.

**`symox.polys.matrices.dense` — 0/17 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `R`, `T`, `ddm_berk`, `ddm_iadd`, `ddm_idet`, `ddm_iinv`, `ddm_ilu`, `ddm_ilu_solve`, `ddm_ilu_split`, `ddm_imatmul`, `ddm_imul`, `ddm_ineg`, `ddm_irmul`, `ddm_irref`, `ddm_irref_den`, `ddm_isub`, `ddm_transpose`.

**`symox.polys.matrices.domainmatrix` — 1/2 actuellement ; `available`.**

- Source retrouvée dans `symox.polys.matrices` : `DM`.

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

Mixte : raccorder les fonctions print_* et les exports pretty existants ; créer la façade codeprinter pour les fonctions déjà disponibles. Les classes C89/C99/C11/FCode/Rust/Python, les exceptions et protocoles communs de CodePrinter, la symbologie et stringPict demandent davantage qu’un réexport.

**`symox.printing.aesaracode` — 0/4 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `AesaraPrinter`, `aesara_code`, `aesara_function`, `dim_handling`.

**`symox.printing.c` — 1/10 actuellement ; `available`.**

- Source retrouvée dans `symox.printing` : `print_ccode`.
- Adaptation / développement à examiner : `C11CodePrinter`, `C89CodePrinter`, `C99CodePrinter`, `CodePrinter`, `get_math_macros`, `known_functions_C89`, `known_functions_C99`, `requires`.

**`symox.printing.codeprinter` — 0/11 actuellement ; `missing_module`.**

- Source retrouvée dans `symox.printing` : `ccode`, `cxxcode`, `fcode`, `print_ccode`, `print_fcode`, `rust_code`.
- Adaptation / développement à examiner : `AssignmentError`, `CodePrinter`, `PrintMethodNotImplementedError`, `print_rust_code`, `requires`.

**`symox.printing.conventions` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `requires_partial`, `split_super_sub`.

**`symox.printing.fortran` — 1/3 actuellement ; `available`.**

- Source retrouvée dans `symox.printing` : `print_fcode`.
- Adaptation / développement à examiner : `FCodePrinter`.

**`symox.printing.jscode` — 2/4 actuellement ; `available`.**

- Source retrouvée dans `symox.printing` : `print_jscode`.
- Adaptation / développement à examiner : `known_functions`.

**`symox.printing.julia` — 2/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `known_fcns_src1`, `known_fcns_src2`, `print_julia_code`.

**`symox.printing.lambdarepr` — 0/6 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `IntervalPrinter`, `LambdaPrinter`, `MpmathPrinter`, `NumPyPrinter`, `PythonCodePrinter`, `lambdarepr`.

**`symox.printing.latex` — 4/7 actuellement ; `available`.**

- Adaptation / développement à examiner : `accepted_latex_functions`, `latex_escape`, `translate`.

**`symox.printing.maple` — 2/3 actuellement ; `available`.**

- Source retrouvée dans `symox.printing` : `print_maple_code`.

**`symox.printing.mathematica` — 2/3 actuellement ; `available`.**

- Adaptation / développement à examiner : `known_functions`.

**`symox.printing.mathml` — 5/6 actuellement ; `available`.**

- Adaptation / développement à examiner : `MathMLPrinter`.

**`symox.printing.octave` — 2/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `known_fcns_src1`, `known_fcns_src2`, `print_octave_code`.

**`symox.printing.precedence` — 2/12 actuellement ; `available`.**

- Adaptation / développement à examiner : `PRECEDENCE_FUNCTIONS`, `PRECEDENCE_VALUES`, `precedence_Float`, `precedence_FracElement`, `precedence_Integer`, `precedence_Mul`, `precedence_PolyElement`, `precedence_Rational`, `precedence_UnevaluatedExpr`, `precedence_traditional`.

**`symox.printing.pretty` — 3/6 actuellement ; `available`.**

- Source retrouvée dans `symox.printing` : `pager_print`, `pprint_try_use_unicode`, `pprint_use_unicode`.

**`symox.printing.pretty.pretty` — 6/19 actuellement ; `available`.**

- Source retrouvée dans `symox.printing` : `pager_print`.
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

**`symox.printing.rcode` — 2/4 actuellement ; `available`.**

- Source retrouvée dans `symox.printing` : `print_rcode`.
- Adaptation / développement à examiner : `known_functions`.

**`symox.printing.rust` — 0/5 actuellement ; `available`.**

- Adaptation / développement à examiner : `RustCodePrinter`, `TypeCast`, `float_ceiling`, `float_floor`, `known_functions`.

**`symox.printing.str` — 2/4 actuellement ; `available`.**

- Source retrouvée dans `symox.printing` : `sstrrepr`.
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

Priorité élevée : raccorder les cinq ensembles numériques du paquet et créer conditionset, fancysets, powerset, sets à partir des objets existants. SetKind, les variantes CartesianComplexRegion/PolarComplexRegion et les auxiliaires d’arithmétique des ensembles nécessitent des adaptateurs ou développements.

**`symox.sets` — 21/26 actuellement ; `available`.**

- Source retrouvée dans `symox` : `Integers`, `Naturals`, `Naturals0`, `Rationals`, `Reals`.

**`symox.sets.conditionset` — 0/3 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `ConditionSet`, `Contains`.
- Adaptation / développement à examiner : `SetKind`.

**`symox.sets.fancysets` — 0/12 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `ComplexRegion`, `Complexes`, `ImageSet`, `Integers`, `Naturals`, `Naturals0`, `Range`, `Rationals`, `Reals`.
- Adaptation / développement à examiner : `CartesianComplexRegion`, `PolarComplexRegion`, `normalize_theta_set`.

**`symox.sets.powerset` — 0/1 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `PowerSet`.

**`symox.sets.sets` — 0/22 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Complement`, `DisjointUnion`, `EmptySet`, `FiniteSet`, `Intersection`, `Interval`, `ProductSet`, `Set`, `SymmetricDifference`, `Union`, `UniversalSet`, `imageset`.
- Adaptation / développement à examiner : `SetKind`, `is_function_invertible_in_set`, `set_add`, `set_div`, `set_function`, `set_mul`, `set_pow`, `set_sub`, `simplify_intersection`, `simplify_union`.

</details>

<a id="simplify"></a>
<details>
<summary>simplify</summary>

Priorité élevée : créer le paquet et les façades vers _simptools, _epath et les fonctions natives ; les 32 exports du paquet existent à la racine. Les 31 exports de simplify.fu, y compris les règles TR*, existent via _simptools et son objet fu_rules. Les internes CSE et les collections/formules hypergéométriques demandent encore des ajouts. Ne pas importer Operator/apply_operators de quantum ou secondquant dans hyperexpand.

**`symox.simplify` — 0/32 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `EPath`, `FU`, `besselsimp`, `collect`, `collect_const`, `combsimp`, `cse`, `denom`, `epath`, `exptrigsimp`, `fraction`, `fu`, `gammasimp`, `hyperexpand`, `hypersimilar`, `hypersimp`, `kroneckersimp`, `logcombine`, `nsimplify`, `numer`, `posify`, `powdenest`, `powsimp`, `radsimp`, `ratsimp`, `ratsimpmodprime`, `rcollect`, `separatevars`, `signsimp`, `simplify`, `sqrtdenest`, `trigsimp`.

**`symox.simplify.combsimp` — 0/1 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `combsimp`.

**`symox.simplify.cse_main` — 0/11 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `cse`.
- Adaptation / développement à examiner : `FuncArgTracker`, `Unevaluated`, `cse_release_variables`, `cse_separate`, `match_common_args`, `opt_cse`, `postprocess_for_cse`, `preprocess_for_cse`, `reps_toposort`, `tree_cse`.

**`symox.simplify.epathtools` — 0/2 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `EPath`, `epath`.

**`symox.simplify.fu` — 0/31 actuellement ; `missing_module`.**

- Source retrouvée dans `symox._simptools.fu_rules` : `L`, `TR0`, `TR1`, `TR10`, `TR10i`, `TR11`, `TR111`, `TR12`, `TR12i`, `TR13`, `TR14`, `TR15`, `TR16`, `TR2`, `TR22`, `TR2i`, `TR3`, `TR4`, `TR5`, `TR6`, `TR7`, `TR8`, `TR9`, `TRmorrie`, `TRpower`, `as_f_sign_1`, `fu`, `hyper_as_trig`, `process_common_addends`, `sincos_to_sum`, `trig_split`.

**`symox.simplify.hyperexpand` — 0/37 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `hyperexpand`.
- Adaptation / développement à examiner : `Formula`, `FormulaCollection`, `G_Function`, `Hyper_Function`, `MeijerFormula`, `MeijerFormulaCollection`, `MeijerShiftA`, `MeijerShiftB`, `MeijerShiftC`, `MeijerShiftD`, `MeijerUnShiftA`, `MeijerUnShiftB`, `MeijerUnShiftC`, `MeijerUnShiftD`, `MultOperator`, `Operator`, `ReduceOrder`, `ShiftA`, `ShiftB`, `UnShiftA`, `UnShiftB`, `add_formulae`, `add_meijerg_formulae`, `apply_operators`, `build_hypergeometric_formula`, `debug`, `devise_plan`, `devise_plan_meijer`, `hyperexpand_special`, `make_derivative_operator`, `make_simp`, `reduce_order`, `reduce_order_meijer`, `try_lerchphi`, `try_polynomial`, `try_shifted_sum`.

**`symox.simplify.hyperexpand_doc` — 0/0 actuellement ; `missing_module`.**

- Aucun nom compté ; vérifier le chemin et les objets non mesurés (constantes ou paquet de documentation).

**`symox.simplify.powsimp` — 0/2 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `powdenest`, `powsimp`.

**`symox.simplify.radsimp` — 0/17 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `collect`, `collect_const`, `denom`, `fraction`, `numer`, `radsimp`, `rcollect`.
- Adaptation / développement à examiner : `collect_abs`, `collect_sqrt`, `denom_expand`, `expand_denom`, `expand_fraction`, `expand_numer`, `fraction_expand`, `numer_expand`, `rad_rationalize`, `split_surds`.

**`symox.simplify.ratsimp` — 0/2 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `ratsimp`, `ratsimpmodprime`.

**`symox.simplify.simplify` — 0/21 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `besselsimp`, `hypersimilar`, `hypersimp`, `kroneckersimp`, `logcombine`, `nsimplify`, `posify`, `separatevars`, `signsimp`, `simplify`.
- Adaptation / développement à examiner : `clear_coefficients`, `dotprodsimp`, `factor_sum`, `inversecombine`, `nc_simplify`, `nthroot`, `product_mul`, `product_simplify`, `sum_add`, `sum_combine`, `sum_simplify`.

**`symox.simplify.sqrtdenest` — 0/6 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `sqrtdenest`.
- Adaptation / développement à examiner : `SqrtdenestStopIteration`, `is_algebraic`, `is_sqrt`, `sqrt_biquadratic_denest`, `sqrt_depth`.

**`symox.simplify.trigsimp` — 0/5 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `exptrigsimp`, `trigsimp`.
- Source retrouvée dans `symox._simptools` : `futrig`.
- Adaptation / développement à examiner : `trigsimp_groebner`, `trigsimp_old`.

</details>

<a id="solvers"></a>
<details>
<summary>solvers</summary>

Raccorder les récurrences, solve_poly_system, solve_triangulated, factor_system, reduce_inequalities, diophantine et les autres fonctions existantes. simplex et ses erreurs nécessitent un véritable solveur LP. Les internes des ODE/systèmes et les inversions de solveset demandent des ajouts ; homogeneous_order est déjà disponible.

**`symox.solvers` — 29/46 actuellement ; `available`.**

- Source retrouvée dans `symox` : `Complexes`, `check_assumptions`, `diophantine`, `factor_system`, `failing_assumptions`, `reduce_inequalities`, `rsolve`, `rsolve_hyper`, `rsolve_poly`, `rsolve_ratio`, `solve_poly_system`, `solve_triangulated`, `solve_univariate_inequality`, `substitution`.
- Adaptation / développement à examiner : `linprog`, `lpmax`, `lpmin`.

**`symox.solvers.ode.ode` — 4/23 actuellement ; `available`.**

- Source retrouvée dans `symox` : `homogeneous_order`.
- Adaptation / développement à examiner : `check_linear_2eq_order1`, `check_nonlinear_2eq_order1`, `check_nonlinear_2eq_order2`, `check_nonlinear_3eq_order1`, `check_nonlinear_3eq_order2`, `classify_sysode`, `constant_renumber`, `get_numbered_constants`, `iter_numbered_constants`, `ode_1st_power_series`, `ode_2nd_power_series_ordinary`, `ode_2nd_power_series_regular`, `ode_sol_simplicity`, `odesimp`, `solve_ics`, `sysode_linear_2eq_order1`, `sysode_nonlinear_2eq_order1`, `sysode_nonlinear_3eq_order1`.

**`symox.solvers.recurr` — 0/4 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `rsolve`, `rsolve_hyper`, `rsolve_poly`, `rsolve_ratio`.

**`symox.solvers.simplex` — 0/6 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `InfeasibleLPError`, `UnboundedLPError`, `linprog`, `lpmax`, `lpmin`, `show_linprog`.

**`symox.solvers.solveset` — 5/14 actuellement ; `available`.**

- Source retrouvée dans `symox` : `NonlinearError`, `linear_eq_to_matrix`, `substitution`.
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

Raccorder Idx, Indexed, IndexedBase, shape et les méthodes d’indices. Les expressions abstraites ArraySymbol/ArrayContraction/ArrayTensorProduct et leurs conversions sont distinctes des tableaux concrets existants. PartialDerivative tensoriel requiert aussi une interface dédiée.

**`symox.tensor` — 13/19 actuellement ; `available`.**

- Source retrouvée dans `symox` : `Idx`, `Indexed`, `IndexedBase`, `get_contraction_structure`, `get_indices`, `shape`.

**`symox.tensor.array.expressions` — 0/16 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `ArrayAdd`, `ArrayContraction`, `ArrayDiagonal`, `ArrayElement`, `ArrayElementwiseApplyFunc`, `ArraySymbol`, `ArrayTensorProduct`, `OneArray`, `PermuteDims`, `Reshape`, `ZeroArray`, `array_derive`, `convert_array_to_indexed`, `convert_array_to_matrix`, `convert_indexed_to_array`, `convert_matrix_to_array`.

**`symox.tensor.index_methods` — 0/3 actuellement ; `missing_module`.**

- Source retrouvée dans `symox._public_matrix` : `IndexConformanceException`, `get_contraction_structure`, `get_indices`.

**`symox.tensor.indexed` — 0/4 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `Idx`, `Indexed`, `IndexedBase`.
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

Raccorder les six exports manquants du paquet et créer enumerative depuis _enumerative (six objets existants). lambdify existe mais lambdastr/implemented_function demandent un contrat. autowrap/ufuncify et les objets CodeGen/Routine exigent un backend/adaptateur supplémentaire. symox.codegen est un module, pas la fonction utilities.codegen.codegen.

**`symox.utilities` — 19/25 actuellement ; `available`.**

- Source retrouvée dans `symox` : `filldedent`, `lambdify`, `public`, `threaded`, `timed`, `xthreaded`.

**`symox.utilities.autowrap` — 0/9 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `CodeWrapError`, `CodeWrapper`, `CythonCodeWrapper`, `DummyWrapper`, `F2PyCodeWrapper`, `UfuncifyCodeWrapper`, `autowrap`, `binary_function`, `ufuncify`.

**`symox.utilities.codegen` — 0/16 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `Argument`, `CCodeGen`, `CodeGen`, `DataType`, `FCodeGen`, `InputArgument`, `JuliaCodeGen`, `OctaveCodeGen`, `OutputArgument`, `Result`, `Routine`, `RustCodeGen`, `codegen`, `default_datatypes`, `get_default_datatype`, `make_routine`.

**`symox.utilities.decorator` — 3/10 actuellement ; `available`.**

- Source retrouvée dans `symox` : `memoize_property`.
- Adaptation / développement à examiner : `T`, `conserve_mpmath_dps`, `deprecated`, `doctest_depends_on`, `no_attrs_in_subclass`, `threaded_factory`.

**`symox.utilities.enumerative` — 0/6 actuellement ; `missing_module`.**

- Source retrouvée dans `symox.utilities._enumerative` : `MultisetPartitionTraverser`, `PartComponent`, `factoring_visitor`, `list_visitor`, `multiset_partitions_taocp`, `part_key`.

**`symox.utilities.exceptions` — 0/3 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `SymPyDeprecationWarning`, `ignore_warnings`, `sympy_deprecation_warning`.

**`symox.utilities.lambdify` — 0/3 actuellement ; `missing_module`.**

- Source retrouvée dans `symox` : `lambdify`.
- Adaptation / développement à examiner : `implemented_function`, `lambdastr`.

**`symox.utilities.memoization` — 0/2 actuellement ; `missing_module`.**

- Adaptation / développement à examiner : `assoc_recurrence_memo`, `recurrence_memo`.

**`symox.utilities.misc` — 1/13 actuellement ; `available`.**

- Adaptation / développement à examiner : `Undecidable`, `as_int`, `debug`, `debug_decorator`, `debugf`, `find_executable`, `func_name`, `ordinal`, `rawlines`, `replace`, `strlines`, `translate`.

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
