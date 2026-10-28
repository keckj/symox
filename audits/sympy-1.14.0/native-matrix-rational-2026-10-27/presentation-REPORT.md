# Suivi R6 : présentation des expressions et éligibilité PolyMatrix

Livraison finale 44 : deux fichiers Python figés dans `native-polymatrix-protocol-r6-20261027/work`. Aucun Rust, paquet antérieur ou principal modifié par cet agent. Base R6 initiale encore en attente d'acceptation : source `04d221ece837895da0cdafd2ca80a3a737faf848f73ea5275dbfd5613b846634`, natif `cc6d021087f8c3d8404011ecf59947d071dfa7ef9bb865c2ae7bd44431cc09ab` ; le snapshot initial et ses échecs restent immuables.

## Diagnostic vérifié

Le shard19 complet reproduit **57 échecs / 931 passes**, assertions de route uniquement. Le probe est exécuté à la fin de la collecte puis avant le premier PolyMatrix. Dans les deux cas :

- Un seul enregistrement faux : `_core.Expr`, dont seuls `__str__` et `__repr__` ont changé.
- Aucun lookup ou état de fonction capturé modifié ; vérificateur vector/Ore True.
- Registre converter rigoureusement identique, objet/ordre/contenu compris : 0 entrée avant/après.

Minimal frais : l'import officiel de `symox.physics.units` installe les wrappers `_install_abbreviated_printing` et suffit à faire refuser l'ancienne garde. La collecte de `test_physics_vector_mechanics` importe ces méthodes avant l'exécution du shard. Le balayage des exports lazy du fichier Gaussian n'est pas requis pour reproduire le refus.

## Correction et preuve de non-consommation

L'élimination ground admissible, exacte sur QQ et fonctions rationnelles formelles, ne présente jamais ses expressions. Les **valeurs de deux seuls attributs**, `CoreExpr.__str__` et `CoreExpr.__repr__`, sont donc exclues de la comparaison des enregistrements concernés. Leurs éventuels records FunctionType sont aussi exclus dès la capture, évitant une vérification accidentelle de leur code/closure si les méthodes étaient déjà Python au bootstrap.

Les deux clés restent obligatoirement présentes. La suppression d'une méthode d'impression ne peut ainsi pas masquer l'ajout d'un autre hook à longueur de dictionnaire égale. Toutes les autres valeurs de classe, bases, fonctions et closures mathématiques, les attributs, constructeurs, convertisseurs et wrappers arithmétiques conservent leurs gardes précédentes. Aucun whitelist de module ajouté, aucun calcul mathématique modifié.

Les nouveaux tests substituent des impressions à compteur et des impressions qui lèvent, pour QQ/RF/QQ(y,z)/ZZ(y), RREF et nullspace : **zéro appel**, aussi bien sur le candidat natif que sur l'ancien calcul complet forcé. La vérification indépendante du parent observe le même résultat et est archivée dans `independent-root-*`.

Les chemins où l'impression est réellement observable restent intacts : erreur de générateurs dupliqués, annulation utilisateur qui imprime puis échoue, conteneur d'items mal formé dont le getter imprime puis échoue. Types/messages d'erreur et ordre/nombre des callbacks correspondent au parcours ancien ; une impression hostile conserve précisément son exception.

## Identités finales

| Fichier | Base SHA256 | Final SHA256 |
|---|---|---|
| `python/symox/polys/_polymatrix_native.py` | `272f7d4522cdb3c4c614ed227538cf05c1073a3624a5c470978dfd9e6f447047` | `35d51e2a5d2d667066dd18ffbfb62707a6b52188f9161e07e768b5c0dd481146` |
| nouveau `tests/test_polymatrix_presentation_protocol.py` | absent | `6b94078b4e8968563dc7bbbde54e372f52a3ecaa5097b5f0953e358d72f04ab4` |

Patch `ffe15698a39f732a1341a53d323a3ab65a3f3aebde56d61b86de001329fc1c2e` ; manifeste des 2994 entrées `ee296b15bbc273b7b736868eab6c5bf07116f1a8c67173021efdb78e9066ba51`. Natif cc6d inchangé, aucune reconstruction requise.

Au contrôle d'intégration, le parent avait déjà copié les deux fichiers : hashes finaux strictement égaux et `git apply --reverse --check` réussi. Le patch avant intégration ne s'applique donc plus en sens direct ; aucune action d'écriture n'a été effectuée sur le principal par ce contrôle.

## Gates

- **184 ciblés locaux passent** : 44 nouveaux +86 ground +26 vector/Ore +28 Stats, toutes assertions existantes inchangées.
- Shard19 complet corrigé localement : **988 passes**.
- Demeter, shard19 exact avec le plugin d'enregistrement du parent : **988 passes**.
- Demeter, shard19 +44 nouveaux +26 vector/Ore +28 Stats +37 coordinate/vector +64 LRT/Atan : **1187 passes**.
- Les listes collectées et terminées sont égales, sans doublon ; 2994 hashes vérifiés avant/après les deux gates distants puis à nouveau localement.
- Six familles de contrôles : absence d'impression sur les quatre domaines, imports unités/Stats/vector/Ore et recharges RT, bootstrap avec impressions FunctionType, erreurs avec impression réelle, refus arithmétique/attribut inchangé, refus suppression-de-clé masquant un hook.

Logs et commandes : `baseline-shard19.log`, `shard19-fixed.log`, `focused-final.log`, `demeter-final.log`, `demeter/{shard19.log,coupling.log,*-command.json,*-execution.json}`. Rejeu : `run_demeter.py`. Copie privée : `/home/keckj/symox-ci/wt/polymatrix-protocol-r6-20261027`, venv Python 3.10 du parent.

Le premier candidat 42 (`8cee65be…`) passait déjà 988/1185 mais n'était pas annoncé figé. Ses deux fichiers, patch/manifeste et tous ses logs bruts ont été conservés séparément dans `preliminary-42`. La version finale 44 ajoute uniquement la présence obligatoire des clés et ses deux tests négatifs.

Ce suivi corrige une sélection de protocole, sans promesse nouvelle de parité globale ou de vitesse. La suite complète R6-r1 sera rejouée par le parent. Nouveau port sparse GJ toujours suspendu après lectures seulement ; replay AUTO R6 préparé mais aucun benchmark autorisé/exécuté pendant cette correction.
