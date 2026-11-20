# Current gap entries to reconcile after limits2

Reconciled after the final replay: the active TODO and series README now
record 1,240/1,293, exactly three original nseries corrections, and the
retained branch/lazy/multiscale gaps. The proposals below are retained as
review history. The 13 remaining selected cases are individually recorded
in `remaining-selected-gaps.json`: seven explicit errors and six output/order
differences, rather than 13 declines.

Status: the current-entry corrections were applied to TODO.md after the final
1240/1293 replay. The remaining 13-case class was checked individually:
seven explicit errors and six expression/order differences; see
`remaining-selected-gaps.json`. The proposal below is retained as review history.

Read-only proposal for the integrating agent. No shared TODO/parity file
was edited. Historical wave summaries and their measured scores must stay
unchanged. Apply only after the final limits2 receipts have been reviewed.

## TODO.md: current pending-deviation list

At the active list containing `piecewise_fold w23rust`, remove the obsolete
claim that `limit(subfactorial(x)/factorial(x), x, oo)` is held. Its present
result is the proved `1/E`, delivered in wave 46. Keep the intentional
reference disagreement in the existing completed wave-46 certificate
entry; this is not an unresolved numerical limit or a new wave-47 gain.

Proposed replacement for that fragment:

```text
piecewise_fold w23rust ; règle Abs/sign(c) de
limitsrest r2 identique à SymPy même pour y complexe ; trois asserts de
```

The completed wave-46 entry in section 4 already links the correct
certificate and retains the unresolved subfactorial domains. No change to
historical wave-40 results in `docs/algorithmic-parity.md` is proposed.

## TODO.md: current limits/series gap inventory

The active gap list still treats all of `nseries` bug2, issue3463 and
log_power1 as unsupported. The original 32-test file proves 25 actual
finite-order baseline failures are now fixed, with its seven existing
guards preserved. Reconcile the named original corpus assertions against
`raw/limits2-deltas.json` before claiming literal corpus matches.

Replace the blanket symbolic/log exponent bullet with this bounded text:

```text
  - séries à exposants symboliques ou logarithmiques : un noyau Rust
    supplémentaire traite des développements sélectionnés d'ordre fini
    en 0+, dont les régressions bug2, issue3463 et log_power1 (25 gains
    vérifiés sur la base figée). Restent les cas multiscales aseries
    simple 2, 6, 7, series 14384:3 et 6350 ; les exposants non ordonnables,
    branches non certifiées et erreurs multipliées par un logarithme non
    borné gardent les anciens replis ;
  - les itérateurs series(n=None), nseries(n=None) et lseries conservent
    tous le moteur baseline. Le tri natif certifié des exposants
    généralisés reste ouvert : log(2**53 - 1)/log(2) est strictement
    inférieur à log(2**53 + 1)/log(2), mais les deux valent 53 en binary64 ;
    l'ancien tri inversait les termes. Les routes cdir non certifiées
    conservent également le moteur antérieur, sans nouveau calcul sur
    une branche positive supposée ;
```

Do not infer new values for `56/1293`, 16 general declines, 17 Gruntz
declines or 30 form-only differences from the 25 focused gains. Refresh
only figures/categories supported by the final replay. Value deltas and
source-line deltas are separate; only elapsed time (row index 1) is ignored.
Trace/stat changes likewise remain distinct from changed mathematical
outputs or newly equal assertion sides.

## TODO.md section 4: explicit completed subset and remaining work

Insert the following near the existing completed subfactorial entry:

```text
- [x] Ajouter un noyau Rust de séries généralisées d'ordre fini en 0+
  pour des exposants réels exacts non rationnels certifiables (vague 47).
  Les 25 gains finis du fichier initial et les gardes d'admission font
  partie des 66 tests séries, tous réussis dans les 596 tests ciblés.
- [ ] Porter le tri natif certifié des nouvelles séries lazy et étendre
  leurs branches/directions ; les itérateurs conservent pour l'instant
  le moteur baseline. Compléter les aseries multiscales simple 2/6/7 et
  les cas series 14384:3/6350. Les réponses indécidables restent explicites.
```

Keep the existing open Ei/erf, Stirling, harmonic, Piecewise/infinity and
subfactorial-outside-certificate entries. This lot does not close them.

## docs/algorithmic-parity.md: current wave-47 section

The current wave-47 section already states finite-order admission, the
complex branch guard, baseline lazy routes, the exact ordering gap,
multiscale gaps and the 596/66 successful checks. No rewrite is needed.
After Rust/corpus receipts, replace only its pending-validation sentences
with the actual completed checks and limits2 results. Do not alter the
historical limits2 wave-43 snapshot, where the older gaps/scores were true
at that time; the wave-47 section supplies the current superseding scope.

The proved boundary distinction stays explicit: mathematical success 2/2,
literal SymPy Expr.series agreement 1/2; SymPy Expr.nseries supplies the
stronger O(x**alpha), while the native result certifies requested O(x**53).
