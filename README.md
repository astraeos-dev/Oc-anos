# Océanos

Calcul harmonique des marées pour les côtes françaises, **100 % hors ligne** : aucun appel réseau, aucune API, aucune prédiction embarquée. Les hauteurs et les coefficients sont calés sur l'annuaire officiel du SHOM. Le moteur `tildes.py` donne pour Cherbourg-en-Cotentin (et Brest) les pleines mers (PM), les basses mers (BM), les hauteurs d'eau au-dessus du zéro hydrographique et les coefficients de marée. Il n'utilise que la bibliothèque standard de Python (3.9 ou plus récent).

## Utilisation

```bash
python3 tildes.py                    # aujourd'hui et demain, Cherbourg
python3 tildes.py 2027-03-21 -n 7    # 7 jours à partir du 21/03/2027
python3 tildes.py -p BREST           # port de Brest
python3 -m unittest                  # tests
```

Comme bibliothèque :

```python
from datetime import date, datetime, timedelta
from zoneinfo import ZoneInfo
import tildes

tildes.extrema_journee(date(2027, 3, 21))        # [{"type", "heure", "hauteur", "coeff"}, ...]
t = datetime(2027, 3, 21, 12, 0, tzinfo=ZoneInfo("Europe/Paris"))
tildes.hauteur_maree(t)                          # hauteur en m au-dessus du zéro hydrographique
tildes.courbe(t, t + timedelta(days=1), 10)      # marégramme au pas de 10 min
```

Les heures sont rendues en heure légale (`Europe/Paris`), mais tout le calcul se fait en UTC. Un `datetime` naïf est lu comme de l'UTC.

Sous Windows, Python n'embarque pas la base des fuseaux horaires. `tildes.py` applique alors sa propre règle de l'heure légale française (`HeureLegaleFrance`, règle de l'UE en vigueur depuis 1996), qui donne les mêmes heures que `zoneinfo` de 2020 à 2030. Rien n'est donc à installer pour `Europe/Paris`. Pour un autre fuseau, installez `tzdata` (`py -m pip install tzdata`).

## Méthode

h(t) = Z0 + Σ f · H · cos(V(t) + u − G)

- **V(t)** est l'argument astronomique de chaque onde (nombres de Doodson, convention de Schureman). Il est calculé à l'instant t à partir des longitudes moyennes de la Lune, du Soleil, du périgée lunaire, du nœud et du périhélie (polynômes de Meeus).
- **f et u** sont les corrections nodales de l'OHI (liste des constituantes du TWCWG, annexe A). Les ondes composées héritent de celles de leurs ondes mères (annexe B).
- **Satellites nodaux.** À Cherbourg, la modulation de M2 sur le cycle nodal de 18,6 ans est plus faible que ne le prévoit la théorie d'équilibre, comme ailleurs en Manche et en mer du Nord. Avec les seules formules de l'OHI, l'heure des PM oscillait de ± 2 à 3 min au rythme de ce cycle. Les raies situées à ± N' de M2, N2 et K2 sont donc ajustées librement, ce qui supprime cette oscillation.
- **H et G** viennent de notre propre analyse harmonique par moindres carrés des hauteurs d'eau observées par les marégraphes REFMAR. Les données utilisées sont les mesures validées en temps différé, au pas de 10 min : 1,61 million de mesures sur 1995-2025 à Cherbourg, 1,36 million sur 2000-2025 à Brest. Le modèle ajusté est exactement celui de `tildes.py` (mêmes V, f et u). Il compte 145 ondes à Cherbourg et 127 à Brest :
  - les 50 ondes de TICON-4 ;
  - 40 ondes de petits fonds et astronomiques (ter-, quart-, sixième- à douzième-diurnes, MO3, MK3, ML4, 3MS4, M10, etc.) ;
  - les 6 satellites nodaux ;
  - les raies de marée trouvées dans le spectre du résidu (amplitude ≥ 2 mm) : 49 à Cherbourg et 31 à Brest. Ce sont des ondes de petits fonds d'ordre élevé (SKM2, 3M2S2, 2MK2, MSL6…), des modulations saisonnières et des satellites. Chacune reçoit la correction nodale de la combinaison d'ondes mères la plus plausible.
- **Z0** est le niveau moyen officiel du SHOM : 3,87 m à Cherbourg et 4,14 m à Brest. C'est celui des Références altimétriques maritimes, publié dans les métadonnées des marégraphes REFMAR, et celui de l'annuaire officiel, dont les hauteurs sont ainsi reproduites. La mer réelle est aujourd'hui plus haute : le niveau moyen observé attendu en 2027 est de 3,938 m à Cherbourg et 4,278 m à Brest (`NIVEAU_MOYEN_OBSERVE_2027`, tendance des 21 dernières années : +6,0 et +6,9 mm/an). Pour prédire la mer réelle plutôt que l'annuaire, il suffit de donner ces valeurs à `z0` et, à Cherbourg, de reprendre les constantes REFMAR sans calage (`CONSTANTES_REFMAR_CHERBOURG`).
- **Calage du marnage (Cherbourg).** Le marnage de l'annuaire est un peu plus grand que celui de nos constantes REFMAR : de 0,9 % en mortes-eaux et de 1,3 % en vives-eaux. Les amplitudes de M2 et de S2 sont donc multipliées par 1,0135 et 1,0116 (`CALAGE_ANNUAIRE_CHERBOURG`). Ces facteurs sont ajustés sur 54 PM et BM officielles du 15 au 28 novembre 2026. Des réglages plus riches (N2, K2, K1, O1, phases, ondes de petits fonds) collaient mieux à la semaine ajustée, mais se dégradaient sur l'autre : ce sont les seuls qui tiennent en validation croisée.
- **Extrema** : ce sont les zéros de la dérivée analytique dh/dt. Ils sont encadrés toutes les 10 min puis affinés par bissection. Deux extrema voisins séparés de moins de 2 cm (simple ondulation) sont écartés.
- **Coefficient** : il suit la définition française. Il est calculé à Brest avec les seules ondes semi-diurnes : C = 100 · A / U, où A est l'amplitude de la pleine mer semi-diurne et U = 3,05 m l'unité de hauteur. Chaque pleine mer locale reçoit le coefficient de la pleine mer de Brest de la même onde de marée, environ 4 h 07 plus tôt pour Cherbourg. Le SHOM calcule ce coefficient avec ses propres constantes de Brest. Les écarts entre nos coefficients et les siens oscillaient au rythme du mois anomalistique (27,5 jours), signe que N2 surtout y diffère de l'analyse TICON-4. `BREST_COEFFICIENT` part donc de TICON-4 et remplace M2, S2, N2, K2, L2 et NU2 par des valeurs effectives. Elles sont ajustées sur 232 coefficients officiels (novembre 2026 à février 2027) et sur la liste ouverte des grandes marées de Saint-Malo (2014-2028).

## Précision

### Accord avec l'annuaire officiel (SHOM)

C'est la référence. Prédictions officielles de Cherbourg du 15 au 28 novembre 2026 : une semaine de mortes-eaux puis une de vives-eaux (coefficients de 34 à 97), soit 54 PM et BM :

| Écart à l'annuaire officiel | Version précédente | Version actuelle |
| --- | --- | --- |
| Heure des PM (biais / RMS / max) | +0,0 / 2,4 / 5,6 min | +0,7 / 2,4 / 4,7 min |
| Heure des BM (biais / RMS / max) | +0,6 / 1,5 / 2,4 min | +0,6 / 1,5 / 2,3 min |
| Hauteur des PM (biais / RMS / max) | +4,4 / 4,8 / 8,1 cm | +0,0 / 1,7 / 3,5 cm |
| Hauteur des BM (biais / RMS / max) | +9,0 / 9,2 / 11,4 cm | −0,2 / 1,2 / 3,1 cm |
| Hauteurs à moins de 1 % de la valeur officielle | 20 sur 54 | 52 sur 54 |
| Coefficients exacts | 9 sur 27 | 20 sur 27 |

Les deux hauteurs au-delà de 1 % sont des BM, à 1,15 % (3,1 cm) et 1,22 % (1,7 cm sur une BM de 1,36 m). Une BM de vives-eaux d'environ 1 m ne tolère que 1 cm d'écart, à peine plus que l'arrondi au centimètre de l'annuaire.

Coefficients officiels de Cherbourg, de novembre 2026 à février 2027 (232 PM, de 32 à 110) :

| Écart aux coefficients officiels | Version précédente | Version actuelle |
| --- | --- | --- |
| Coefficients exacts | 41 % | 74 % (71 % en validation croisée par mois) |
| Coefficients à ± 1 | 89 % | 100 % |
| Biais | +0,50 | −0,02 |
| Grandes marées 2014-2028 (coefficient ≥ 100, Saint-Malo) | 443 trouvées sur 443, 93 jours en trop | 442 sur 443, 47 jours en trop |

La validation croisée cale les constantes du coefficient sur trois mois et les teste sur le quatrième. Les jours de grandes marées en trop sont des jours à 99 officiels que le calcul donne à 100.

### Accord avec la mer observée

La mesure se fait sur la mer réelle, sur des années que l'ajustement n'a jamais vues : constantes calculées sur 1995-2023 (Brest : 2000-2023), puis comparées aux observations du 1er janvier 2024 au 12 novembre 2026. Sur ces 34 mois, Cherbourg compte environ 2 000 PM et 2 000 BM. Les PM et BM observées sont obtenues sans aucun modèle : un polynôme de degré 6 est ajusté sur ± 2 h de mesures autour de chaque extremum. Sur une courbe synthétique, ce repérage a un biais inférieur à 0,2 min.

**Cherbourg**, écart aux observations 2024-2026 :

| Écart aux observations | Version précédente (TICON-4) | Version actuelle (REFMAR) |
| --- | --- | --- |
| Heure des PM (biais / RMS) | +4,3 min / 8,1 min | +0,5 min / 4,6 min |
| Heure des BM (biais / RMS) | +0,7 min / 4,8 min | +0,7 min / 4,2 min |
| PM à ± 5 min / ± 10 min | 51 % / 82 % | 84 % / 96 % |
| BM à ± 5 min / ± 10 min | 77 % / 96 % | 84 % / 97 % |
| PM de mortes-eaux, coefficient < 55 (biais / RMS) | +9,1 min / 12,2 min | +0,7 min / 6,9 min |
| PM de vives-eaux, coefficient ≥ 75 (RMS) | 4,7 min | 2,5 min |
| Hauteur des PM, hors surcote (RMS) | 6,2 cm | 5,1 cm |
| Hauteur des BM, hors surcote (RMS) | 5,7 cm | 4,5 cm |
| Niveau moyen observé − Z0 | +6,4 cm | +7,1 cm (Z0 officiel) |
| Résidu de hauteur, hors surcote (écart-type) | 6,4 cm | 5,1 cm |

**Brest**, écart aux observations 2024-2026 :

| Écart aux observations | Version précédente (TICON-4) | Version actuelle (REFMAR) |
| --- | --- | --- |
| Heure des PM (biais / RMS) | +1,5 min / 4,5 min | −0,3 min / 3,8 min |
| Heure des BM (biais / RMS) | −0,2 min / 4,7 min | −0,2 min / 3,7 min |
| PM et BM à ± 10 min | 96 % | 98 % |
| Hauteur des PM / BM, hors surcote (RMS) | 4,3 cm / 5,1 cm | 3,4 cm / 4,0 cm |
| Niveau moyen observé − Z0 | +11,7 cm | +14,0 cm (Z0 officiel) |

Les constantes publiées dans `tildes.py` ont ensuite été recalculées sur toutes les années validées (jusqu'à 2025). Sur janvier-novembre 2026, toujours hors apprentissage, elles donnent :
- **Cherbourg** : PM à 4,6 min RMS, BM à 4,3 min ;
- **Brest** : PM et BM à 3,7 min RMS.

Avec le niveau moyen observé (`NIVEAU_MOYEN_OBSERVE_2027`) au lieu du niveau officiel, le niveau moyen de 2026 est reproduit à −0,7 cm près à Cherbourg et à −1,2 cm près à Brest.

Ces tableaux portent sur les constantes REFMAR sans calage. Calées sur l'annuaire, les hauteurs de Cherbourg s'écartent un peu plus de la mer réelle de 2026 : 6,1 cm RMS au lieu de 5,1 pour les PM, 4,7 cm au lieu de 4,1 pour les BM, hors surcote. Le marnage de l'annuaire est en effet plus grand que celui observé. Les horaires, eux, ne changent pas.

Dans ces écarts, la « hauteur hors surcote » retire la surcote météorologique, estimée par la moyenne glissante sur 25 h du résidu.

L'erreur restante vient surtout de la météo, qu'aucun modèle astronomique ne prévoit. À Cherbourg, une surcote de plus de 30 cm retarde la BM de 6 min en moyenne, et les plus gros écarts (30 à 60 min) surviennent pendant les tempêtes, par exemple en janvier 2025 et en janvier 2026. Les mortes-eaux restent les plus difficiles (6,9 min RMS), parce que la courbe y est très plate autour de la pleine mer.

Pour mémoire, la version à constantes TICON-4 avait été validée sur l'annuaire officiel d'octobre 2026 à octobre 2027 (prédictions SHOM publiées par Météo Consult, 1 531 PM et BM à Cherbourg). Elle donnait 4,7 min RMS sur l'heure des PM, 4,3 min sur celle des BM, 3,4 cm sur la hauteur des PM et 4,6 cm sur celle des BM, avec des PM de mortes-eaux 4 à 5 min trop tard.

## Limites

- Ce sont des prédictions astronomiques. Les effets météorologiques n'y figurent pas : surcote ou décote, environ 1 cm par hPa d'écart à 1013 hPa, plus le vent. Ils décalent aussi les heures (voir plus haut).
- Z0 est le niveau moyen officiel, comme dans l'annuaire. La mer réelle est en moyenne 7 cm plus haute à Cherbourg et 14 cm plus haute à Brest, et l'écart grandit d'environ 6 mm par an. Le Z0 de Brest n'a pas été contrôlé contre des hauteurs officielles de Brest.
- L'échantillon officiel d'heures et de hauteurs couvre deux semaines de novembre 2026 (un cycle mortes-eaux puis vives-eaux). Le calage du marnage est à revérifier sur d'autres saisons et dans quelques années : on ne sait pas si l'écart avec le SHOM est constant ou s'il suit le cycle nodal.
- Les horaires restent à 2,4 min RMS de l'annuaire pour les PM (4,7 min au pire) et à 1,5 min pour les BM. Nos PM sont un peu en avance en mortes-eaux et en retard en vives-eaux : cela tient à la forme de la courbe près de la pleine mer (ondes de petits fonds), qu'il faudrait caler sur beaucoup plus de semaines officielles.
- Les satellites nodaux sont estimés sur 31 ans à Cherbourg et 26 ans à Brest, soit 1,4 à 1,7 cycle nodal. Il faut au moins 18,6 ans d'observations pour les séparer de leur onde mère.
- Environ un coefficient sur quatre diffère encore d'un point de la valeur officielle. Les constantes du coefficient sont calées sur quatre mois d'hiver et sur les seuils de grandes marées ; des coefficients officiels de printemps et d'été affineraient l'ajustement.

## Outils de développement

Le répertoire `outils/` sert à recalculer les constantes et à mesurer la précision. Il demande `numpy`, contrairement à `tildes.py`, qui n'en dépend pas.

```bash
python3 outils/refmar.py 13 3 1995 2025          # observations validées de Cherbourg (3 pour Brest)
python3 outils/refmar.py 13 2 2026 2026          # observations brutes de l'année en cours
python3 outils/ajuster.py 13 1995 2025           # analyse harmonique → code à reporter dans tildes.py
python3 outils/evaluer.py 2024-01-01 2026-11-12  # écart de tildes.py aux observations
```

- `refmar.py` télécharge les observations du service data.shom.fr dans `donnees/refmar/`. Ce répertoire est ignoré par git : les observations ne sont pas redistribuées.
- `harmonique.py` réalise l'ajustement par moindres carrés avec les mêmes arguments astronomiques et corrections nodales que `tildes.py`.
- `ajuster.py` enchaîne nettoyage, ajustement, recherche des raies dans le spectre du résidu, ajustement final et choix de Z0.
- `evaluer.py` compare les PM et BM de `tildes.py`, ou d'un jeu de constantes candidat (`--candidat`), à celles observées.

## Historique des corrections

Première série, corrections de la version initiale :

1. **Jours de changement d'heure.** L'arithmétique `datetime + timedelta` en heure locale faisait reculer le temps d'une heure. Le 29/03/2026, le programme créait une fausse PM à 02 h 52 et une fausse BM à 03 h 02. Le calcul se fait désormais en UTC.
2. **Coefficient.** L'ancienne version utilisait la marée de Cherbourg et une unité de hauteur locale de 2,62 m, et évaluait la composante semi-diurne à l'heure de la PM complète. Le coefficient est maintenant calculé à Brest selon la définition officielle.
3. **Constantes harmoniques.** Les valeurs arrondies, présentées comme « officielles » sans source, ont été remplacées par 48 ondes TICON-4 sourcées et 12 ondes de petits fonds inférées, au lieu de 16 ondes. Par exemple, la situation de S2 était fausse de 4°.
4. **Corrections nodales.** L2 utilisait celles de M2, alors que la sienne dépend du périgée. Les formules complètes de l'OHI remplacent les approximations au premier ordre.
5. **Recherche des extrema.** L'interpolation parabolique sur un pas de 5 min est remplacée par la recherche des zéros de la dérivée analytique.
6. **Affichage.** L'en-tête « Heure (CET) » était faux l'été ; le décalage UTC réel est maintenant affiché à chaque ligne.

Seconde série, recalcul sur observations :

7. **Constantes recalculées sur les observations REFMAR** de Cherbourg et de Brest, avec 145 et 127 ondes, dont les ondes de petits fonds mesurées au lieu d'être inférées, les satellites nodaux et les raies trouvées dans le résidu. Le retard des PM de mortes-eaux à Cherbourg (+9 min face à la mer observée) disparaît.
8. **Niveau moyen.** Z0 était calculé sur 2016-2021 à Cherbourg et sur 2006-2009 à Brest. C'est désormais le niveau attendu en 2027 : +6,2 cm à Cherbourg et +11,5 cm à Brest.
9. **Décalage Brest-port.** L'écart de phase de M2 était pris modulo 360°. Pour Brest elle-même, un écart légèrement négatif aurait donné 12,4 h au lieu de 0, et donc le coefficient de la marée précédente. Il est maintenant ramené entre −60° et 300°.
10. **Windows.** Sans base de fuseaux horaires, `tildes.py` applique sa propre règle de l'heure légale française.

Troisième série, calage sur l'annuaire officiel :

11. **Niveau moyen officiel.** Z0 passe au niveau moyen du SHOM : 3,87 m à Cherbourg et 4,14 m à Brest. À Cherbourg, l'écart de hauteur à l'annuaire tombe de +6 à +8 cm à environ ± 1 cm.
12. **Coefficients.** Les constantes de Brest servant au coefficient (`BREST_COEFFICIENT`) sont ajustées sur les coefficients officiels : 74 % d'exacts au lieu de 41 %, 100 % à ± 1 au lieu de 89 %, et aucun biais.
13. **Marnage de Cherbourg.** M2 et S2 sont relevés de 1,35 % et 1,16 %. Sur deux semaines officielles, 52 hauteurs sur 54 sont à moins de 1 %, au lieu de 20.

## Sources et licences

- Observations marégraphiques : réseau REFMAR (SHOM et partenaires), marégraphes de Cherbourg (n° 13) et de Brest (n° 3), données validées en temps différé, téléchargées sur data.shom.fr. Le site refmar.shom.fr indique que ses contenus sont sous licence Etalab 2.0 « sauf mention contraire ». Seules les constantes harmoniques que nous en avons tirées figurent dans ce dépôt.
- Niveaux moyens officiels : SHOM, Références altimétriques maritimes, via les métadonnées des marégraphes REFMAR (services.data.shom.fr).
- Calage du marnage de Cherbourg : 54 PM et BM officielles du 15 au 28 novembre 2026 (pages de maree.info, © SHOM), utilisées pour l'ajustement seulement et non reproduites ici.
- Calage des coefficients :
  - 232 coefficients officiels de Cherbourg, de novembre 2026 à février 2027 (calendrier de maree.info, « d'après les prédictions officielles © SHOM »), utilisés pour l'ajustement seulement et non reproduits ici ;
  - liste des grandes marées de la Ville de Saint-Malo, 2014-2028 (data.gouv.fr, Licence Ouverte).
- Liste des ondes de base et point de départ des constantes du coefficient (Brest) : Hart-Davis, M., Dettmering, D., Seitz, F. (2025). *TICON-4: TIdal CONstants based on GESLA-4 sea-level records*. SEANOE. https://doi.org/10.17882/109129 — licence CC BY 4.0. Elles sont extraites via [slackwater-database](https://github.com/openwatersio/slackwater-database), qui ramène aussi les phases de SGM, T3, R3, 3L2, 2MK5 et 2MO5 à la convention de Schureman.
- Corrections nodales et numéros de Doodson : IHO Tidal and Water Level Working Group, *Constituent List* (annexes A et B) ; Schureman, *Manual of Harmonic Analysis and Prediction of Tides* (1958).
- Astronomie : J. Meeus, *Astronomical Algorithms*.
- Code : licence MIT (voir `LICENSE`).
