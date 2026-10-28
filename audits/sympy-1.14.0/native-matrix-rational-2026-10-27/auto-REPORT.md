# DomainMatrix AUTO : statistiques et choix ZZ/QQ natifs

Livraison des six fichiers source figée dans `native-domainmatrix-auto-20261027/work`. Aucun principal ni ancien paquet figé modifié, aucun commit, maturin ou installation partagée. Le binaire de cette copie demeure **DEBUG privé**, réservé aux preuves de correction. Une release des mêmes sources est en compilation séparée ; aucune mesure ni accélération revendiquée à ce stade.

## Chemin public et portée

`DomainMatrix.rref()` et `.rref_den()` appellent `_dm_rref_choose_method`. Le nouveau raccord AUTO stock ZZ/QQ exécute en Rust le comptage des lignes/entrées, les tailles binaires, le LCM des dénominateurs et les seuils GJ/FF/CD. Il réutilise les noyaux FF/GJ et les préparations QQ/clear_denoms déjà natifs ; pas de duplication des éliminations. Python ajoute seulement un raccord et l'authentification des protocoles. Les anciens calculs Python restent complets pour les domaines, formats, données et protocoles non reconnus. Toutes les méthodes explicites, dont `_dense`, conservent leur dispatch antérieur.

- Sparse : les lignes littéralement vides et les zéros stockés participent exactement aux statistiques originales. Les branches courtes n'inspectent pas les coefficients ou clés que l'ancien sélecteur ne lisait pas.
- Dense : les zéros et lignes nulles sont exclus comme dans `to_sdm`. Seuls les coefficients natifs déjà typés pour ce domaine sont admis ; des entiers Python injectés ou autres valeurs conservent la conversion complète.
- Les seuils reproduisent l'ordre des opérations binary64 de Python. `int/int` passe par la rationnelle exacte et son arrondi indépendant certifié existant ; aucun double arrondi des opérandes. Les produits de tailles binaires utilisent BigInt, les comptes d'entrées u128.
- Le LCM QQ conserve l'ordre des valeurs et le même arrêt anticipé, les seuils `5*numer_bits` et 50 bits. Les données ZZ mal formées de bits maximum 0 gardent précisément le `ZeroDivisionError` original via le repli, sans méthode inventée.
- Les dimensions supportées par usize sont lues sans cast i64. Au-delà, le calcul Python original reste possible ou lève son erreur propre ; aucun quota nouveau, refus mathématique ou densification.

## Protocoles et capacités

La garde vise les méthodes et données utilisées : `to_sdm`, `to_flat_nz`, `to_dok`, le constructeur SDM et `_native_elements`, flags/domaines, bool/int et numérateur/dénominateur. Elle ne capture aucun dictionnaire entier d'Expr : Stats/vector/Ore et les attributs de matrices sans rapport ne changent pas l'éligibilité. Une recharge de rref de code canonique est authentifiée avec code marshal v0, namespace, défauts et cellules ; les recharges ne bénissent pas des callbacks modifiés.

La construction dense lit `convert`, mais `of_type` natif réussit sur tous les scalaires admis. La lecture d'une fonction ordinaire est donc inerte ; un descripteur utilisateur qui exécute `__get__` garde son parcours et ses callbacks. Le binding refuse les valeurs non typées avant toute conversion. Les callbacks des classes, instances, fonctions, builtins, clés, conteneurs et coefficients sont conservés sans préflight observable ; `KeyboardInterrupt` est propagé. Aucun domaine composite, RR, CC, EX ou GF n'est remplacé par ZZ/QQ. Les classes DDM/SDM/DomainMatrix recréées par recharge non authentifiée peuvent conserver le repli complet : aucun gain n'est annoncé sur ce cas.

## Identités et intégration

Base immuable : `native-literal-derivatives-r1-20261027/after`, source `2eb95e7788af6673c6bf1f8dfd2faed142d976880d62f577948ffc2a9b4702d7`, release de base `3243129bfe3e4ec0611ebe7761b5988f5bba3312df8a9e9911fea5c7a785eaee`.

Six fichiers, et seulement ces six dans le patch :

1. `crates/symox-core/src/poly/fraction_free.rs` : statistiques, seuils, LCM, quatre tests Rust.
2. `crates/symox-py/src/poly_fraction_free.rs` : lecture CAPI conservatrice et `_dm_rref_auto`.
3. `python/symox/polys/matrices/rref.py` : adaptation AUTO, ancien sélecteur conservé.
4. `python/symox/polys/matrices/domainmatrix.py` : enregistrement des classes complètes.
5. Nouveau `python/symox/polys/matrices/_rref_native.py` : protocoles seulement.
6. Nouveau `tests/test_domainmatrix_native_auto.py` : 505 cas.

Les quatre bases existantes ont été vérifiées identiques au principal au gel. Les hashes base/final complets sont dans `identity.json`. Patch : `ef6aa7352ed38558a1536063cb378cbad76cd4999121514444c178dd2dc39339`. Manifeste DEBUG des 2771 entrées : `6c1f5123b94f4dc58b9753a344d568b69f25313a38a1990f2f45d01875a04240`. Binaire DEBUG privé : `323aad58aadf1ce7dca6e2db390f5eac401075053d442eea3cb21d05c9036532`. La différence du binaire est documentée séparément et n'appartient pas au patch source.

## Preuves

- **505 nouveaux tests Python** : 384 comparaisons de décisions avec le parcours Python original et SymPy1.14 ; 30 véritables appels publics RREF/RREF_DEN à GJ/FF/CD, domaines/types/formats/dénominateurs/pivots et certificat exact de rang/espace de lignes ; autres cas protocoles, données invalides, zéros, dimensions, interruptions et domaines hors portée.
- Six permutations Stats/vector/Ore en processus frais, trois recharges de rref à chaque fois ; maintien de la route native stock. Les callbacks avant/après recharge restent observables.
- **1248 ciblés locaux verts**, dont les 505 nouveaux et 743 existants DD(M)/DomainMatrix/clear_denoms, Stats/holonomic/vector.
- **1248 ciblés Demeter verts**, liste collectée exactement égale à la liste terminée, sans doublon. 2771 hashes vérifiés avant et après, puis recontrôlés localement après rapatriement.
- **4 tests Rust verts**, dont 35 fixtures de statistiques/seuils calculées indépendamment par CPython, incluant des entiers au-delà de 2^53, et les interruptions du LCM.
- Cargo check et clippy exit 0 ; aucun diagnostic clippy dans les ajouts. Les diagnostics existants du dépôt sont visibles dans le journal intégral.
- `git apply --check` contre le principal : succès. Aucun test existant modifié.

Logs : `coupling-local.log`, `cargo-auto-tests-final.log`, `cargo-check.log`, `clippy.log`, `demeter/targeted.log`, `demeter/execution.json`, `demeter/command.json`, `demeter-run.log`. Rejeu distant : `run_demeter.py`, Python 3.10 privé `/home/keckj/symox-ci/venv/bin/python`, copie `/home/keckj/symox-ci/wt/domainmatrix-auto-20261027`.

## Contrôle de performance séparé

La release privée en cours doit produire une nouvelle copie `after-release`, sans remplacer le DEBUG figé. Les mêmes entrées et certificats exacts serviront aux appels publics `.rref()`/`.rref_den()` avant/après, dans une fenêtre calme autorisée. Aucun ratio de composant, aucun résultat DEBUG, aucune réduction de capacité ne sera présenté comme accélération. La suite complète du paquet intégré reste à la charge du parent sur Demeter.
