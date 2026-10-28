# Payload Float0 des polynômes EX — livraison isolée v18

## Base, périmètre et hashes

Stage : `expression_zero_payload`, nouvelle copie du main incluant calcul sparse v21, FF et radsimp0/0. Base immuable : `expression_zero_payload_base`. Aucun changement de main ou d'un ancien stage ; aucun commit, maturin ou installation partagée.

La baseline a été réellement compilée en DEBUG privé sur Demeter depuis cette base ; aucune extension transitoire du main n'est utilisée comme preuve.

- Manifeste base : 2776 fichiers, SHA256 `b5ddbcb082de35cc26ee640f954de93136b19eb91abc5597905c39f3dffc7655`.
- Binaire baseline : `7658355be59d944da48273c22eb90ea9db7bf6d987d46317d1dd748531cc16ca`.
- Manifeste final `source-v18.sha256` : 2777 fichiers, SHA256 `27690544f92b547bc5292395dd58ff2e16fe837e626340fcc1390809b7a9d62b`.
- Binaire final : `7d39429b63a408e6bdb71a30126621daf0308764c2e7cd63c0678621ff3b71e8`.
- Patch `changes.patch` : SHA256 `c5f387732b48c31e048737b2b67dd6da15dbf2aae5a60c56a04a336d9198d488`.

`git apply --check` passe sur la base immuable. `file-hashes.json` donne les hashes base/final des neuf fichiers : expression_domain.rs, fraction.rs, bindings poly.rs ; _poly_constructor.py, _purepoly.py, polys/_evaluation.py, polys/_poly_methods.py, polys/polyclasses.py ; nouveau test_expression_zero_payload_parity.py. Aucun TODO/contrat partagé n'est modifié.

## Comportement corrigé

`Poly.from_dict({(5,): 1, (0,): Float(0, precision=139)}, x, domain=EX)` garde désormais le Float0 brut. Les listes de coefficients, TC/nth et eval(0) rendent le Float avec précision 139, alors que as_dict/as_expr l'omettent correctement du support mathématique. Un zéro Float en tête et une feuille multivariée entièrement nulle sont éliminés comme dans le DMP de référence ; seul un zéro intérieur ou trailing d'une feuille non nulle survit.

La construction et cette normalisation sont natives. Un constructeur Rust reçoit les coefficients structurés avant qu'un Add/Mul d'expression puisse absorber Float0 ou fusionner des générateurs formels. Le noyau retire les têtes/lignes zéro légitimes avant d'estimer Dense/Sparse : un monôme nul d'exposant 2**50 ne force plus une allocation dense. Les doublons gardent la dernière valeur, comme la construction ordinaire antérieure.

Les payloads survivants restent sur le parcours Dense complet établi. Un export natif distinct du support non nul alimente les protocoles de représentation et pickle. TC/nth lisent les coefficients réellement stockés ; all_coeffs garde le contrôle dense antérieur puis rétablit leur payload. La préparation des listes/dictionnaires DMP, la validation des coefficients natifs et les bornes sont dans le binding Rust. Python transporte les objets et restaure les classes ; aucun nouvel algorithme mathématique Python.

Les transitions gardent les conventions de la référence : ajout de zéro et négation peuvent conserver Float0 ; le produit numérique fini Float0*1 en EX produit le zéro entier. Ce dernier cas était masqué par la perte au constructeur : EX.mul normalise désormais ce seul produit numérique nul, sans remplacer `0*inf` par `0`. Une addition avec un opérande sparse conserve le repli Dense lorsque le payload l'exige ; après normalisation arithmétique du payload, le chemin sparse redevient admissible. Les cinq précisions sont vérifiées exactement.

Les conversions identitaires EX, to_field/to_ring et copies stock conservent le contenu. Une conversion EX→QQ convertit le zéro en zéro entier ; QQ→EX ne réinvente pas de Float0. Le pickle stock Poly/PurePoly garde classe, domaine et payload. Il s'agit volontairement d'un round trip plus fidèle que celui de SymPy 1.14 : le pickle de référence peut perdre le domaine EX et Float0 ; son résultat réel reste enregistré dans les preuves, sans le présenter comme une égalité de pickle.

## Raccord d'évaluation et protocoles

L'adaptateur univarié reconstruisait sa valeur depuis as_expr(), qui omet légitimement Float0. L'évaluation en zéro numérique stock passe maintenant au noyau natif avec l'index déjà résolu. Cela évite aussi de prendre un générateur formel numérique comme 2 pour l'index 2. Les générateurs formels `0`, `2`, `1/2`, `Float('1.25')` et les trois variables sont comparés à la référence.

Le raccourci reste réservé au Poly stock, au domaine EX stock, aux méthodes/domaines non remplacés et aux coefficients/générateurs natifs déterministes. La sélection des méthodes utilise leurs dictionnaires de classe sans invoquer les descripteurs de remplacement. Les sous-classes, overrides de as_expr/as_dict, conversions de domaine remplacées et descripteurs personnalisés gardent l'ancien parcours et ses callbacks. Les évaluations non nulles gardent aussi leurs chemins.

Les coefficients host, fonctions numériques susceptibles d'appeler le host, contextes de budget/interruption et non-finis prennent les anciennes routes complètes. Ils n'effectuent pas de nouvelle passe is_zero et ne sont pas déclarés couverts par le correctif de payload. Le détenteur DMP historique garde également les dictionnaires à coordonnées plus larges que u32, les métadonnées de niveau très grandes et les listes/protocoles étrangers. Refuser un raccourci n'est jamais présenté comme un gain de capacité ni comme un résultat nul.

## Preuves et validation privée Demeter

- 88 nouveaux cas permanents : précision 1/53/139/140/211 ; Float0 et coefficients non nuls comparés par leurs valeurs binaires exactes ; une à trois variables ; têtes/lignes zéro, doublons, grands exposants, transitions Dense/Sparse ; LC/TC/nth/all_coeffs, as_expr, eval de zéros numériques, conversion, vues DMP, pickle et PurePoly ; générateurs formels numériques ; callbacks et descripteurs ; absence d'import runtime de SymPy.
- Comparaison avant/après de 15 jeux d'entrées : tous les inputs binaires sont identiques. Les 10 divergences de payload de la baseline sont corrigées ; les cinq cas où une tête zéro doit disparaître restent corrects. Voir baseline.json, final-v18.json et baseline-comparison.json.
- Les traces d'un véritable coefficient host sont identiques à la baseline pour construction, représentation, eval(0), conversion EX et copie. Voir host_probe.py, host-baseline.json et host-final-v18.json.
- Une sélection de régressions existantes : 1133 passages avant les dernières gardes conservatrices (`regressions-v9.log`). Elle n'est pas relancée en boucle.
- Tests natifs de polynômes : 221 passages (`native-poly-final.log`), avant la dernière garde portant sur les générateurs du raccourci. Le noyau final repasse les deux tests ciblés de payload (`native-v18.log`).
- Sur le dernier binaire/source : 437 passages (`focused-v18.log`), dont les 88 nouveaux et les contrôles existants des constructeurs, domaines et hooks/lifetime/reuse host. Aucun échec, skip ou xfail ajouté.

Les sources sont vérifiées avant compilation et après les derniers tests ; le SHA du binaire est vérifié après exécution. `validation-link.json` lie base, source, binaire et résultats. Les extensions existent seulement dans les deux arbres privés `/home/keckj/symox-ci/wt/expression-zero-payload[-base]-20261026/python/symox/`. Aucun benchmark privé. La suite complète de l'instantané intégré reste coordonnée par le parent.

## Limites précises conservées

Ce lot ferme le cas stock EX décrit ci-dessus, pas tous les comportements EX :

- Les mélanges avec coefficients/générateurs host, fonctions numériques et autres contextes exclus gardent les anciennes limites de payload et leurs protocoles complets.
- Les sous-classes, notamment l'évaluation PurePoly, conservent les replis d'évaluation antérieurs. La construction, les vues et le pickle PurePoly sont testés ; son alignement/égalité spécifique n'est pas porté ici.
- RR→EX peut créer des Float0 implicites dans la référence. Ce port de conversion n'est pas fait ; les chemins existants restent en place. La sélection implicite d'un domaine influencée par Float0 reste également distincte du stockage EX explicite/déjà inféré.
- Le constructeur DMP brut __new__ conserve ses types et son chemin antérieur ; son filtre historique de zéros n'est pas réécrit. Les vues Poly.rep et les constructeurs DMP de conversion sont le périmètre réparé.
- Un Float0 intérieur survivant peut nécessiter la représentation Dense existante, donc une mémoire proportionnelle au degré. Aucune capacité sparse pour ces payloads bruts n'est revendiquée, et aucun quota arbitraire ne masque ce coût.
- Les algèbres tenues sous evaluate(False) ne sont pas globalement portées : la référence peut conserver Mul(Float0,1) en EX, alors que l'arithmétique du noyau reste canonique. Le constructeur public précédent perdait déjà Float0 dans ce contexte ; aucune parité globale de ces opérations n'est annoncée.

Aucune revendication de parité globale EX/EXRAW ou de tous les domaines de paramètres ne découle de ces résultats.
