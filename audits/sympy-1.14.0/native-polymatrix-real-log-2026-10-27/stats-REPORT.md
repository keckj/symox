# Suivi R5 : installation tardive des méthodes Stats et PolyMatrix

Livraison distincte et figée : `native-polymatrix-protocol-r5-20261027/work`. Le principal, le suivi 213 et les autres paquets figés n'ont pas été modifiés. Base complète vérifiée sur les 2867 entrées du manifeste R5 : source `ceeaa3a15a3dab2b1961f41bc766cadfa4e11c1b3fd862c33623d02fa4133123`, extension release `532386bda46046d3d21c99f49c1b863fac62bdee5c357d03cf5123bde6913d4b`.

## Diagnostic exact

Le shard42 complet reproduit **4 échecs / 963 passes**. Un probe séparé observe un seul enregistrement de classe refusé : `node_class('Expr')`. L'import officiel de Stats y ajoute `doit`, `rewrite`, `_stats_doit_installed`. Aucune dépendance globale ni fonction capturée ne diffère, et le vérificateur vector/Ore retourne True. Cette classe nominale est distincte de la façade publique `S.Expr`, dont les écritures ciblent justement la classe nominale.

Le parcours ground admis est commutatif, exact et composé de coefficients rationnels ou de fonctions rationnelles formelles. Son annulation passe par `cancel_entry` vers le noyau natif ; il n'appelle ni `doit` ni `rewrite`. Les nœuds statistiques et coefficients hors ce parcours restent refusés par le binding existant et suivent le calcul complet.

## Correction ciblée

- Le bootstrap capture la classe nominale et les deux méthodes héritées d'origine.
- La différence admise est exactement l'ajout des deux méthodes officielles et du marqueur vrai. Toutes les anciennes entrées/bases/fonctions restent vérifiées.
- Les fonctions doivent avoir les codes lexicaux canoniques de `stats._symbolic.install`, le namespace de ce module, aucun défaut ajouté, la cellule attendue et le délégué d'origine exact. La sélection n'importe pas Stats et n'exécute pas ses méthodes.
- Stats ajoute seulement trois constantes privées de code et leur extraction. **L'AST de `install` est strictement identique avant/après**, branches mathématiques et registre HELD compris.
- Les overrides publics, changements de code, cellules, défauts, marqueur, module, templates, installateur et attributs sans rapport conservent le repli complet.

## Trois fichiers à intégrer

| Fichier | SHA256 base | SHA256 final |
|---|---|---|
| `python/symox/polys/_polymatrix_native.py` | `68ec973b02d152e2cde41b995325c1b7fda587aa4b11f2c272d68e87af1b3304` | `272f7d4522cdb3c4c614ed227538cf05c1073a3624a5c470978dfd9e6f447047` |
| `python/symox/stats/_symbolic.py` | `dc9a2cae04abf321036cf67705ecf370cf38356361c1382a28890fa6cc132395` | `0619cc85813800018b0854ce6b820345850403c47faa607538673d5898a1f468` |
| nouveau `tests/test_polymatrix_stats_protocol.py` | absent | `63156a3a470aab371f2032e8ea3f56b698eb05abd28d98d24e6de90238dec7a3` |

Patch `d2131e4498336794edbfe7fcbb43935a1e8c28bee93803601093cce973592739` ; manifeste final des 2868 entrées `30391129d7bf0a746b988668d6f6318aa376a7943437b53fbfd052096c7632a3`. RT optimisé `c966bb72e3f6b83b6ff2fbff151f4096a7ba04eb4955517b0ac0bcd990b6ec27` et binaire 532 inchangés. Aucun Rust ni reconstruction nécessaire. `git apply --check` sur le principal : succès.

## Preuves

- Local, shard42 complet : **967 passes**. Base avec les quatre refus conservée dans `baseline-shard42.log` ; correction dans `shard42-fixed.log`.
- Local, nouveaux 28 + anciens 112 PolyMatrix : **140 passes** (`focused-final-v2.log`).
- Demeter, même shard42 et plugin du parent : **967 passes** (`demeter/shard42.log`).
- Demeter, shard42 + nouveaux 28 + PolyMatrix 86 + coordinate 37 + LRT/Atan 64 : **1182 passes** (`demeter/coupling.log`).
- Pour les deux runs distants : liste collectée exactement égale à la liste terminée, aucun doublon ;2868 hashes vérifiés avant et après, puis de nouveau localement après rapatriement.
- Six ordres d'import Stats/vector/Ore dans des processus frais, sans import de Stats pendant la sélection ; mêmes résultats natifs exacts avant/après.
- QQ et QQ(y,z), RREF et base complète de noyau comparés à SymPy1.14 et au parcours Python conservé ; certificats exacts `M*v=0`.
- Registre HELD hostile non consulté sur ce parcours, ni natif ni Python complet ; un vrai coefficient statistique conserve précisément son erreur `CoercionFailed`.
- Toutes les assertions des tests existants, notamment les 26 du suivi 213 et les fichiers du shard42, restent inchangées.

Copie privée Demeter : `/home/keckj/symox-ci/wt/polymatrix-protocol-r5-20261027`. `run_demeter.py`, les commandes JSON, journaux et listes exact-once sont archivés dans `proof`. Le rejeu utilise le venv3.10 et le binaire R5 existant. Aucun benchmark réalisé. Ce correctif concerne l'éligibilité, sans nouvelle prétention de parité globale ; la suite complète combinée après intégration reste gérée par le parent.

AUTO est suspendu dans sa copie distincte `native-domainmatrix-auto-20261027` : aucun benchmark ; son premier cargo check local privé était vert. Aucun de ses changements n'est inclus ici.
