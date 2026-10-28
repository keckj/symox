# Derivative class query repair

Private repair stage: `branches/symbolic-derivative-atoms-fix`; private Demeter stage: `/home/keckj/symox-ci/wt/next-derivative-atoms-fix`.

The new native `SymbolicDerivative` registers the public name `Derivative` in `Func::from_name`. The atoms adapter previously interpreted the public Derivative class as an ordinary callable function head. It therefore omitted the established machine-count `Kind::Derivative`, breaking ODE/PDE auto-detection and downstream algorithms that call `atoms(Derivative)`.

The repair changes only `crates/symox-py/src/expression_inspection.rs` and adds `tests/test_derivative_class_queries.py`:

- Keep Derivative class selectors on the existing Python isinstance protocol, covering both representations and custom metaclass protocols. A Python callable merely named Derivative retains its previous atom-filter protocol.
- Exclude SymbolicDerivative from `function_node`, consistently with `is_Function=False`. This also restores the instance-selector convention: `atoms(an_instance_of_Derivative)` selects all Derivative instances, across both native representations.

The tests compare reference1.14 classes, instances, multiple filters, tuple atom filters, machine/symbolic/BigInt/Float/rational counts, complete selected sets, find/count/group/has, custom instancecheck callbacks and errors, Function exclusion, and ODE/PDE automatic function detection. No mathematical derivative evaluator changes. No main source was edited by this agent.

Two has gaps remain explicitly unchanged: `expr.has((Derivative,Symbol))` and `expr.has(a_python_callable)` raise TypeError, where SymPy returns False. Both refusals are reproduced in the older native-algebra/after package (native SHA9405dccb462e37ed095dc19b597595b74662ce45d23556e96bf8d52e836636a2), before the new representation. The tests preserve these errors; no has implementation is changed.

Private final validation:

- DEBUG binding build succeeded; core sources were not edited.
- 940 passed in `pytest-final-inspection.log`: new13 + existing expression_inspection + symbolic_derivative_orders + solver_checks + the complete differential_equation_interfaces module.
- 15 passed in `pytest-final-query.log`, including new13, solver_checks held-derivative regression and holo_rest_parity::test_multiple_derivative_final_step. Combined with the preceding run: 941 unique tests passed.
- The earlier exploratory `pytest-focused.log` loaded an initial test file with six assertions for the two preexisting has gaps, and eventually hit a diagnostic timeout on unrelated late solve/nonlinsolve tests in the holo module. It is not a gate or a claim of complete success. Final logs above use the final test file.

Sources are frozen. `source-manifest.json` records exact base/final hashes. Patch `derivative-class-queries.patch` has SHA16e89141eabfc481ca6da6a4b7a6efa0bb27c9acf3eb71ec3800654c5d6212ce. It is based on main release02fa2def and does not include the parent's independent derivative_vars lookahead fix. Parent is responsible for the combined release and full-suite gate.
