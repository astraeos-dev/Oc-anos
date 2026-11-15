# Océanos

Calcul harmonique des marées pour les côtes françaises, **100 % hors ligne** : aucun appel réseau, aucune API, aucune dépendance au SHOM. Le moteur `tildes.py` donne pour Cherbourg-en-Cotentin (et Brest) les pleines mers (PM), les basses mers (BM), les hauteurs d'eau au-dessus du zéro hydrographique et les coefficients de marée. Il n'utilise que la bibliothèque standard de Python (3.9 ou plus récent).

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
- **Satellites nodaux.** À Cherbourg, la modulation de M2 sur le cycle nodal de 18,6 ans est plus faible que ne le prévoit la théorie d'équilibre, comme ailleurs en Manche et en mer du Nord. Avec les seules formules de l'OHI, l'heure des PM oscillait de ± 2 à 3 min au rythme de ce cycle. Les raies situées à ± N' de M2, N2 et K2 sont donc ajustées librement, ce qui supprime cette oscillation.
- **H et G** viennent de notre propre analyse harmonique par moindres carrés des hauteurs d'eau observées par les marégraphes REFMAR. Les données utilisées sont les mesures validées en temps différé, au pas de 10 min : 1,61 million de mesures sur 1995-2025 à Cherbourg, 1,36 million sur 2000-2025 à Brest. Le modèle ajusté est exactement celui de `tildes.py` (mêmes V, f et u). Il compte 145 ondes à Cherbourg et 127 à Brest :
  - les 50 ondes de TICON-4 ;
  - 40 ondes de petits fonds et astronomiques (ter-, quart-, sixième- à douzième-diurnes, MO3, MK3, ML4, 3MS4, M10, etc.) ;
  - les 6 satellites nodaux ;
  - les raies de marée trouvées dans le spectre du résidu (amplitude ≥ 2 mm) : 49 à Cherbourg et 31 à Brest. Ce sont des ondes de petits fonds d'ordre élevé (SKM2, 3M2S2, 2MK2, MSL6…), des modulations saisonnières et des satellites. Chacune reçoit la correction nodale de la combinaison d'ondes mères la plus plausible.
- **Z0** est le niveau moyen attendu en 2027 : 3,938 m à Cherbourg et 4,278 m à Brest. Il est déduit de la tendance des moyennes annuelles observées sur les 21 dernières années (+6,0 et +6,9 mm/an). La hausse s'est accélérée : une tendance sur 30 ans sous-estimerait le niveau actuel.
- **Extrema** : ce sont les zéros de la dérivée analytique dh/dt. Ils sont encadrés toutes les 10 min puis affinés par bissection. Deux extrema voisins séparés de moins de 2 cm (simple ondulation) sont écartés.
- **Coefficient** : il suit la définition française. Il est calculé à Brest avec les seules ondes semi-diurnes : C = 100 · A / U, où A est l'amplitude de la pleine mer semi-diurne et U = 3,05 m l'unité de hauteur. Chaque pleine mer locale reçoit le coefficient de la pleine mer de Brest de la même onde de marée, environ 4 h 07 plus tôt pour Cherbourg. Ce calcul garde les constantes TICON-4 de Brest (`BREST_TICON`), sur lesquelles il a été validé contre les coefficients officiels (voir plus bas). Les coefficients sont donc inchangés.

## Précision

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
| Niveau moyen observé − Z0 | +6,4 cm | +1,9 cm |
| Résidu de hauteur, hors surcote (écart-type) | 6,4 cm | 5,1 cm |

**Brest**, écart aux observations 2024-2026 :

| Écart aux observations | Version précédente (TICON-4) | Version actuelle (REFMAR) |
| --- | --- | --- |
| Heure des PM (biais / RMS) | +1,5 min / 4,5 min | −0,3 min / 3,8 min |
| Heure des BM (biais / RMS) | −0,2 min / 4,7 min | −0,2 min / 3,7 min |
| PM et BM à ± 10 min | 96 % | 98 % |
| Hauteur des PM / BM, hors surcote (RMS) | 4,3 cm / 5,1 cm | 3,4 cm / 4,0 cm |
| Niveau moyen observé − Z0 | +11,7 cm | +2,3 cm |

Les constantes publiées dans `tildes.py` ont ensuite été recalculées sur toutes les années validées (jusqu'à 2025). Sur janvier-novembre 2026, toujours hors apprentissage, elles donnent :
- **Cherbourg** : PM à 4,6 min RMS, BM à 4,3 min, niveau moyen à −0,7 cm ;
- **Brest** : PM et BM à 3,7 min RMS, niveau moyen à −1,2 cm.

Dans ces écarts, la « hauteur hors surcote » retire la surcote météorologique, estimée par la moyenne glissante sur 25 h du résidu.

L'erreur restante vient surtout de la météo, qu'aucun modèle astronomique ne prévoit. À Cherbourg, une surcote de plus de 30 cm retarde la BM de 6 min en moyenne, et les plus gros écarts (30 à 60 min) surviennent pendant les tempêtes, par exemple en janvier 2025 et en janvier 2026. Les mortes-eaux restent les plus difficiles (6,9 min RMS), parce que la courbe y est très plate autour de la pleine mer.

**Comparaison avec l'annuaire officiel.** La version précédente avait été validée sur l'annuaire officiel (prédictions SHOM publiées par Météo Consult) d'octobre 2026 à octobre 2027. Sur 1 531 PM et BM à Cherbourg, elle donnait 4,7 min RMS sur l'heure des PM et 4,3 min sur celle des BM, 3,4 cm sur la hauteur des PM et 4,6 cm sur celle des BM. Les PM de mortes-eaux y sortaient 4 à 5 min trop tard. Cette comparaison n'a pas été refaite, car les tables ne sont pas dans le dépôt. Par rapport à la version précédente, la version actuelle avance les PM de 3,5 min en moyenne et relève toutes les hauteurs d'environ 6 cm (nouveau Z0). Elle se rapproche ainsi de la mer observée, mais l'écart avec l'annuaire dépend du niveau moyen retenu par le SHOM.

## Limites

- Ce sont des prédictions astronomiques. Les effets météorologiques n'y figurent pas : surcote ou décote, environ 1 cm par hPa d'écart à 1013 hPa, plus le vent. Ils décalent aussi les heures (voir plus haut).
- Z0 est le niveau moyen attendu en 2027 et reste constant dans le calcul. Au rythme actuel (environ +6 mm/an), il faudra le relever d'environ 3 cm tous les cinq ans.
- Les satellites nodaux sont estimés sur 31 ans à Cherbourg et 26 ans à Brest, soit 1,4 à 1,7 cycle nodal. Il faut au moins 18,6 ans d'observations pour les séparer de leur onde mère.
- Le coefficient reste en moyenne 0,8 point au-dessus de la valeur officielle. Les amplitudes semi-diurnes TICON de Brest sont en effet environ 1 % plus fortes que celles qu'impliquent les coefficients publiés. Une unité de hauteur « effective » de 3,088 m annulerait ce biais (46 % de coefficients exacts, 93 % à ± 1), mais il faudrait l'étalonner sur les coefficients publiés. C'est pourquoi la valeur officielle de 3,05 m est conservée.

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

## Sources et licences

- Observations marégraphiques : réseau REFMAR (SHOM et partenaires), marégraphes de Cherbourg (n° 13) et de Brest (n° 3), données validées en temps différé, téléchargées sur data.shom.fr. Le site refmar.shom.fr indique que ses contenus sont sous licence Etalab 2.0 « sauf mention contraire ». Seules les constantes harmoniques que nous en avons tirées figurent dans ce dépôt.
- Liste des ondes de base et constantes du coefficient (Brest) : Hart-Davis, M., Dettmering, D., Seitz, F. (2025). *TICON-4: TIdal CONstants based on GESLA-4 sea-level records*. SEANOE. https://doi.org/10.17882/109129 — licence CC BY 4.0. Elles sont extraites via [slackwater-database](https://github.com/openwatersio/slackwater-database), qui ramène aussi les phases de SGM, T3, R3, 3L2, 2MK5 et 2MO5 à la convention de Schureman.
- Corrections nodales et numéros de Doodson : IHO Tidal and Water Level Working Group, *Constituent List* (annexes A et B) ; Schureman, *Manual of Harmonic Analysis and Prediction of Tides* (1958).
- Astronomie : J. Meeus, *Astronomical Algorithms*.
- Code : licence MIT (voir `LICENSE`).
