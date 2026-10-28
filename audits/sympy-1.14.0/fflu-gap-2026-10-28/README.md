# FFLU : erreurs héritées et correction mathématique démontrée

Les implémentations Python SymPy 1.14.0 et Symox produisent des facteurs faux ou une division par zéro sur certains pivots retardés par des colonnes nulles. Le [rapport gelé](PUBLIC-REPORT.md) donne les deux contre-exemples, les identités exactes, les domaines, les messages et la dérivation par mineurs.

**203 identifiants de fixtures**, dont 201 entiers en ZZ/QQ et 2 fractionnels en QQ, donnent **404 cas de domaine** ; les 30 formes m=0..4/n=0..5 et tous leurs rangs sont présents. Chaque moteur exécute 1 616 appels : **1 408 certifiés, 84 factorisations fausses et 124 exceptions**. Les doublons de permutation sont conservés ; ce corpus représente 164 matrices numériques distinctes.

Le prototype mathématique indépendant corrige les 404 cas, vérifie 1 582 divisions exactes et leurs déterminants de mineurs. L’injection du packed LU corrigé dans la seconde phase originale inchangée donne 404 certificats complets dans chaque moteur. Ces contrôles sont distincts des appels réels qui échouent encore. Les 1 056 résultats publics déjà corrects par moteur conservent exactement leurs facteurs.

**Aucune correction de production, aucun port et aucune mesure de vitesse ne sont revendiqués.** Le futur noyau Rust doit couvrir les deux phases, et le repli générique doit recevoir la correction mathématique. Les cas actuellement faux ou exceptionnels ne seront pas des baselines de vitesse.

[L’archive de 34 fichiers](public-proof.tar.gz), sa [carte d’empreintes](proof-archive-map.json), les scripts et les sorties brutes sont vérifiés. Le snapshot Symox de reproduction est le R11 initial immuable, native 532c4b08/source 1801af86 ; son gate complet échoue séparément sur quatre assertions d’admission PolyMatrix, sans lien avec ces erreurs FFLU héritées. Aucun port FFLU n’a modifié ce snapshot.
