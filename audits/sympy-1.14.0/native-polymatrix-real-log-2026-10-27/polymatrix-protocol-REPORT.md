# PolyMatrix : protocoles vector/Ore et recharge de rationaltools

Livraison finale exclusive : `followup-213/work`. Deux fichiers Python ; aucun Rust, noyau mathématique, domaine ou précision changé. Le principal est resté intact.

## Cause et correction

Le garde PolyMatrix capturait tout le dictionnaire de `_core.Expr`. Les imports officiels `vector` et `holonomic.recurrence` ajoutent respectivement `_coordinate_vector_dispatch_installed` et `_ore_dispatch_installed`, et installent quatre wrappers arithmétiques. Ces différences faisaient refuser la route native même pour les seuls scalaires stock. Reproduction baseline : **57 échecs / 66 passes**, tous les échecs de route, valeurs correctes.

Le garde conserve la vérification des classes, bases, fonctions, dictionnaires, constructeurs et convertisseurs. Seuls ces deux marqueurs vrais et ces quatre méthodes peuvent différer. La composition des wrappers est authentifiée par l'API stable existante `rationaltools._atan_stock_operators(None, None, operators, states)` : code des wrappers, cellules, délégués stock, registries, classes scalaires et absence de callbacks. Les fonctions/dépendances de ce vérificateur sont capturées au bootstrap. Aucun import vector/Ore pendant la sélection et aucune copie d'un ancien vérificateur.

Les tests LRT rechargent légitimement rationaltools. Les helpers ainsi recréés peuvent être reconnus avec le même namespace, les mêmes défauts/cellules et le même code sérialisé canonique. Le format **marshal version 0** évite les références partagées du format 4, dont les marqueurs dépendent de l'état de références au moment de la sérialisation. Aucune égalité ni sérialisation Python d'objet hôte n'est appelée. Un helper modifié ou un code contenant des constantes hôtes refuse la sélection. Les interruptions ne sont pas avalées.

## Identités

- `python/symox/polys/_polymatrix_native.py` : base `a8e9d3d2a2b9e2ce7dfe19a8546abc9325e4fdd31629cf91dc190257e22abd30`, final `68ec973b02d152e2cde41b995325c1b7fda587aa4b11f2c272d68e87af1b3304`.
- Nouveau `tests/test_polymatrix_vector_protocol.py` : `ff27a7f95223347bd3f5edd951a8acfc81f4d95beb08d3e17638f428c02b6c30`.
- Patch : `a0ad1e8d83a81c014e5dcc5e59c7aff23cd1a63bffcb9c9fe9a0eebef53d30aa`.
- Manifeste des 792 entrées Python/tests/extension : `b26f45b6596c460ab7ef01507149354568cd5f10c5726d94b26c92ae4cd8b4f2`.
- RT optimisé existant conservé : `c966bb72e3f6b83b6ff2fbff151f4096a7ba04eb4955517b0ac0bcd990b6ec27`.
- Extension release R5 existante : `532386bda46046d3d21c99f49c1b863fac62bdee5c357d03cf5123bde6913d4b`, aucune reconstruction privée ni installation partagée.

## Preuves

- **213 passent localement**, ordre LRT avant PolyMatrix ; `local-final.log`.
- **213 passent localement**, ordre PolyMatrix avant LRT ; `local-opposite-order.log`.
- **213 passent sur Demeter**, ordre LRT avant PolyMatrix, manifeste des 792 entrées vérifié avant exécution ; `demeter-final.log`.
- `git apply --check` contre le principal courant : succès.
- 86 tests PolyMatrix existants, 37 tests coordinate/vector existants, 64 tests LRT/Atan existants, **26 tests nouveaux**. Assertions existantes inchangées.
- Processus frais pour vector seul, Ore seul, vector→Ore, Ore→vector ; même résultat natif exact avant/après. Repli et ordre de callbacks utilisateur, codes/cellules/bases/predicats/proofs altérés, attribut sans rapport, recharge RT et constantes hôtes contrôlés.

Copie Demeter : `/home/keckj/symox-ci/wt/polymatrix-vector-protocol-final-20261027`. Le script `replay-demeter.sh` rejoue seulement ce ciblé. Aucun benchmark ni grande suite répétés ; validation complète combinée à la charge du parent.

## Chronologie et limites

Le premier paquet annoncé figé à 210 tests avait été réutilisé pour les diagnostics avant l'échec d'ordre sur Demeter. Ses octets ont été récupérés depuis le script de génération et le préfixe des tests ; les deux hashes de fichiers, le patch `1455dd6f…` et le manifeste `28316bec…` sont certifiés identiques à l'annonce initiale dans `../preliminary-210/proof`. Ses logs local210pass et Demeter59fail/151pass sont conservés. Les diagnostics intermédiaires restent visibles dans le dossier parent, notamment l'échec marshal v4. Cette livraison distincte remplace exclusivement ce suivi préliminaire.

Cette correction porte sur l'éligibilité des protocoles stock. Les autres changements de classes/protocoles et les domaines hors noyau existant gardent le calcul Python complet. Aucune extension de parité mathématique ni accélération n'est revendiquée. AUTO DomainMatrix demeure différé et n'a reçu aucune édition.
