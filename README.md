# Océanos

Calcul harmonique des marées pour les côtes françaises, **100 % hors ligne** : aucun appel réseau, aucune API, aucune dépendance au SHOM. Le moteur `tildes.py` donne pour Cherbourg-en-Cotentin les pleines mers (PM), les basses mers (BM), les hauteurs d'eau au-dessus du zéro hydrographique et les coefficients de marée. Il n'utilise que la bibliothèque standard de Python (3.9 ou plus récent).

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

## Méthode

h(t) = Z0 + Σ f · H · cos(V(t) + u − G)

- **V(t)** est l'argument astronomique de chaque onde (nombres de Doodson, convention de Schureman). Il est calculé à l'instant t à partir des longitudes moyennes de la Lune, du Soleil, du périgée lunaire, du nœud et du périhélie (polynômes de Meeus).
- **f et u** sont les corrections nodales de l'OHI (liste des constituantes du TWCWG, annexe A). Les ondes composées héritent de celles de leurs ondes mères (annexe B).
- **H et G** proviennent de TICON-4. Ce jeu de constantes est issu de l'analyse harmonique des enregistrements marégraphiques GESLA-4 : 1943-2025 pour Cherbourg, 1846-2025 pour Brest. 48 de ses 50 ondes sont retenues : MA2 et MB2, modulations saisonnières de M2, sont écartées parce que ce ne sont pas des ondes astronomiques et que les retirer améliore l'accord avec l'annuaire officiel sur chacun des deux semestres testés.
- **12 ondes de petits fonds** absentes de TICON-4 sont inférées (MK4, SN4, SK4, 2MN6, MSN6, 2MK6, 2SM6, MSK6, 3MN8, 2MSN8, 3MS8, 2MS8). L'amplitude d'une onde composée est proportionnelle au produit des amplitudes de ses ondes mères, et sa phase est la somme de leurs phases plus un déphasage propre à son espèce. Chacune est donc déduite de l'onde observée de même espèce et de fréquence la plus proche. Le détail est dans `INFERENCES`.
- **Z0** est le niveau moyen observé au-dessus du zéro hydrographique : 3,876 m à Cherbourg, mesuré sur 2016-2021.
- **Extrema** : ce sont les zéros de la dérivée analytique dh/dt. Ils sont encadrés toutes les 10 min puis affinés par bissection. Deux extrema voisins séparés de moins de 2 cm (simple ondulation) sont écartés.
- **Coefficient** : il suit la définition française. Il est calculé à Brest avec les seules ondes semi-diurnes : C = 100 · A / U, où A est l'amplitude de la pleine mer semi-diurne et U = 3,05 m l'unité de hauteur. Chaque pleine mer locale reçoit le coefficient de la pleine mer de Brest de la même onde de marée, environ 4 h 07 plus tôt pour Cherbourg.

## Précision

La validation a porté sur l'annuaire officiel (prédictions SHOM publiées par Météo Consult) d'octobre 2026 à octobre 2027 : 1 531 PM et BM à Cherbourg, jours de changement d'heure compris. Ces tables ne sont pas embarquées et aucune constante n'a été ajustée dessus : elles ont servi à mesurer l'écart et à trancher deux choix de modèle (ondes inférées, retrait de MA2 et MB2).

| Écart avec l'annuaire officiel | Version initiale | Version actuelle |
| --- | --- | --- |
| Heure des PM (RMS / max) | 12,1 min / 54 min | 4,7 min / 19 min |
| Heure des BM (RMS / max) | 7,4 min / 22 min | 4,3 min / 19 min |
| Heures à ± 5 min | 39 % | 77 % |
| Heures à ± 10 min | 70 % | 97 % |
| Hauteur des PM (RMS / max) | 10,5 cm / 27 cm | 3,4 cm / 10 cm |
| Hauteur des BM (RMS / max) | 10,8 cm / 27 cm | 4,6 cm / 12 cm |
| Coefficient (biais / RMS / max) | +3,7 / 4,7 / 11 | +0,8 / 1,2 / 3 |
| Coefficient exact / à ± 1 | 9 % / 24 % | 30 % / 78 % |
| Extrema parasites | 2 (changement d'heure) | 0 |

Les ondes de petits fonds inférées réduisent à elles seules l'erreur RMS de 25 % sur les heures et de 22 % sur les hauteurs. En vives-eaux (coefficient ≥ 75), l'erreur RMS sur les heures tombe à environ 3 min. Les plus gros écarts (14 à 19 min) surviennent en mortes-eaux, quand la courbe est très plate autour de l'extremum ; les pleines mers y sortent en moyenne 4 à 5 min trop tard.

## Limites

- Ce sont des prédictions astronomiques. Les effets météorologiques n'y figurent pas : surcote ou décote, environ 1 cm par hPa d'écart à 1013 hPa, plus le vent.
- Z0 est un niveau moyen observé ; la hausse du niveau marin (quelques mm par an) n'est pas extrapolée.
- Le coefficient reste en moyenne 0,8 point au-dessus de la valeur officielle. Les amplitudes semi-diurnes TICON de Brest sont en effet environ 1 % plus fortes que celles qu'impliquent les coefficients publiés. Une unité de hauteur « effective » de 3,088 m annulerait ce biais (46 % de coefficients exacts, 93 % à ± 1), mais il faudrait l'étalonner sur les coefficients publiés. C'est pourquoi la valeur officielle de 3,05 m est conservée.

## Corrections apportées à la version initiale

1. **Jours de changement d'heure.** L'arithmétique `datetime + timedelta` en heure locale faisait reculer le temps d'une heure. Le 29/03/2026, le programme créait une fausse PM à 02 h 52 et une fausse BM à 03 h 02. Le calcul se fait désormais en UTC.
2. **Coefficient.** L'ancienne version utilisait la marée de Cherbourg et une unité de hauteur locale de 2,62 m, et évaluait la composante semi-diurne à l'heure de la PM complète. Le coefficient est maintenant calculé à Brest selon la définition officielle.
3. **Constantes harmoniques.** Les valeurs arrondies, présentées comme « officielles » sans source, sont remplacées par 48 ondes TICON-4 sourcées et 12 ondes de petits fonds inférées, au lieu de 16 ondes. Par exemple, la situation de S2 était fausse de 4°.
4. **Corrections nodales.** L2 utilisait celles de M2, alors que la sienne dépend du périgée. Les formules complètes de l'OHI remplacent les approximations au premier ordre.
5. **Recherche des extrema.** L'interpolation parabolique sur un pas de 5 min est remplacée par la recherche des zéros de la dérivée analytique.
6. **Affichage.** L'en-tête « Heure (CET) » était faux l'été ; le décalage UTC réel est maintenant affiché à chaque ligne.

## Sources et licences

- Constantes harmoniques : Hart-Davis, M., Dettmering, D., Seitz, F. (2025). *TICON-4: TIdal CONstants based on GESLA-4 sea-level records*. SEANOE. https://doi.org/10.17882/109129 — licence CC BY 4.0. Elles sont extraites via [slackwater-database](https://github.com/openwatersio/slackwater-database), qui ramène aussi les phases de SGM, T3, R3, 3L2, 2MK5 et 2MO5 à la convention de Schureman.
- Observations sous-jacentes : GESLA-4 (Global Extreme Sea Level Analysis), marégraphes REFMAR de Cherbourg et de Brest.
- Corrections nodales et numéros de Doodson : IHO Tidal and Water Level Working Group, *Constituent List* (annexes A et B) ; Schureman, *Manual of Harmonic Analysis and Prediction of Tides* (1958).
- Astronomie : J. Meeus, *Astronomical Algorithms*.
- Code : licence MIT (voir `LICENSE`).
