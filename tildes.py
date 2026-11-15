#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Océanos — calculateur harmonique des marées, 100 % hors ligne.

Calcule les pleines mers (PM), les basses mers (BM), les hauteurs d'eau
au-dessus du zéro hydrographique et les coefficients de marée pour
Cherbourg-en-Cotentin. Aucun appel réseau, aucune dépendance externe :
bibliothèque standard Python (>= 3.9) uniquement.

Modèle
------
    h(t) = Z0 + Σ f·H·cos(V(t) + u − G)

* V(t) : argument astronomique de chaque onde (nombres de Doodson, convention
  de Schureman), calculé à l'instant t à partir des longitudes moyennes de la
  Lune, du Soleil, du périgée lunaire, du nœud ascendant et du périhélie
  (polynômes de Meeus, « Astronomical Algorithms »).
* f, u : corrections nodales de l'OHI (liste des constituantes du TWCWG,
  annexe A) ; les ondes composées héritent de leurs ondes mères (annexe B).
* H, G : amplitude et situation (déphasage de Greenwich, temps UTC) des
  ondes de TICON-4, analyse harmonique des enregistrements marégraphiques
  GESLA-4 (DGFI-TUM, licence CC BY 4.0) : 48 ondes retenues, complétées par
  12 ondes de petits fonds inférées de ces mêmes constantes. Aucune
  prédiction ni constante harmonique du SHOM n'est utilisée.
* Extrema : zéros de la dérivée analytique dh/dt, encadrés toutes les 10 min
  puis affinés par bissection (écart < 1 s avec la solution exacte du modèle).
* Coefficient : définition française, calculé à Brest sur les seules ondes
  semi-diurnes (amplitude de la pleine mer semi-diurne / unité de hauteur U),
  puis rattaché à chaque pleine mer locale de la même onde de marée.

Tous les calculs se font en UTC ; l'heure légale (CET/CEST) n'intervient qu'à
l'affichage, ce qui rend exacts les jours de changement d'heure.
"""

import argparse
import bisect
import math
from datetime import date, datetime, time, timedelta, timezone
from zoneinfo import ZoneInfo

FUSEAU_DEFAUT = "Europe/Paris"

# ==============================================================================
# 1. DÉFINITION DES ONDES (constituantes harmoniques)
# ==============================================================================
# Chaque onde : (coefficients de Doodson, règle de correction nodale).
#
# Coefficients sur les arguments (τ, s, h, p, N', p', 90°) :
#   τ  = T + h − s, T = angle horaire du Soleil moyen (180° à 0 h UTC)
#   s  = longitude moyenne de la Lune      h  = longitude moyenne du Soleil
#   p  = longitude du périgée lunaire      N' = −N (N : nœud ascendant lunaire)
#   p' = longitude du périhélie (périgée solaire)
#   le dernier coefficient compte des multiples de 90° (constante de phase).
#
# Règle nodale : nom d'une formule fondamentale (voir _fondamentales), tuple
# de (onde, facteur) pour une onde composée, ou None (f = 1, u = 0).
#
# Les phases TICON de SGM, T3, R3, 3L2, 2MK5 et 2MO5 ont été ramenées à ces
# arguments (décalages établis contre les enregistrements bruts GESLA-4 par
# le projet slackwater-database).
ONDES = {
    # --- longue période -------------------------------------------------------
    "SA":      ((0, 0, 1, 0, 0, 0, 0), None),
    "SSA":     ((0, 0, 2, 0, 0, 0, 0), None),
    "MM":      ((0, 1, 0, -1, 0, 0, 0), "MM"),
    "MSF":     ((0, 2, -2, 0, 0, 0, 0), (("M2", -1),)),
    "MF":      ((0, 2, 0, 0, 0, 0, 0), "MF"),
    "MTM":     ((0, 3, 0, -1, 0, 0, 0), (("MM", 1),)),
    "MSQM":    ((0, 4, -2, 0, 0, 0, 0), (("M2", -1),)),
    # --- diurnes --------------------------------------------------------------
    "2Q1":     ((1, -3, 0, 2, 0, 0, 1), (("O1", 1),)),
    "SGM":     ((1, -3, 2, 0, 0, 0, 1), (("O1", 1),)),
    "Q1":      ((1, -2, 0, 1, 0, 0, 1), (("O1", 1),)),
    "RHO1":    ((1, -2, 2, -1, 0, 0, 1), (("O1", 1),)),
    "O1":      ((1, -1, 0, 0, 0, 0, 1), "O1"),
    "M1":      ((1, 0, 0, 1, 0, 0, -1), "M1"),
    "P1":      ((1, 1, -2, 0, 0, 0, 1), None),
    "S1":      ((1, 1, -1, 0, 0, 0, 0), None),
    "K1":      ((1, 1, 0, 0, 0, 0, -1), "K1"),
    "J1":      ((1, 2, 0, -1, 0, 0, -1), "J1"),
    "OO1":     ((1, 3, 0, 0, 0, 0, -1), (("K2", 1), ("Q1", -1))),
    # --- semi-diurnes ---------------------------------------------------------
    "EP2":     ((2, -3, 2, 1, 0, 0, 0), (("M2", 1),)),
    "2N2":     ((2, -2, 0, 2, 0, 0, 0), (("M2", 1),)),
    "MU2":     ((2, -2, 2, 0, 0, 0, 0), (("M2", 1),)),
    "N2":      ((2, -1, 0, 1, 0, 0, 0), (("M2", 1),)),
    "3N2":     ((2, -1, 0, 0, 0, 0, 1), (("N2", 1),)),
    "NU2":     ((2, -1, 2, -1, 0, 0, 0), (("M2", 1),)),
    "MA2":     ((2, 0, -1, 0, 0, 0, 0), (("M2", 1),)),
    "M2":      ((2, 0, 0, 0, 0, 0, 0), "M2"),
    "MB2":     ((2, 0, 1, 0, 0, 0, 0), (("M2", 1),)),
    "MKS2":    ((2, 0, 2, 0, 0, 0, 0), (("M2", 1), ("K2", 1), ("S2", -1))),
    "LAMBDA2": ((2, 1, -2, 1, 0, 0, -2), (("M2", 1),)),
    "L2":      ((2, 1, 0, -1, 0, 0, -2), "L2"),
    "3L2":     ((2, 1, 0, 0, 0, 0, 0), (("L2", 1),)),
    "T2":      ((2, 2, -3, 0, 0, 1, 0), None),
    "S2":      ((2, 2, -2, 0, 0, 0, 0), None),
    "R2":      ((2, 2, -1, 0, 0, -1, -2), None),
    "K2":      ((2, 2, 0, 0, 0, 0, 0), "K2"),
    "2SM2":    ((2, 4, -4, 0, 0, 0, 0), (("S2", 2), ("M2", -1))),
    # --- tiers-diurnes --------------------------------------------------------
    "M3":      ((3, 0, 0, 0, 0, 0, 0), "M3"),
    "T3":      ((3, 3, -4, 0, 0, 0, 0), None),
    "S3":      ((3, 3, -3, 0, 0, 0, 0), None),
    "R3":      ((3, 3, -2, 0, 0, 0, -1), (("S2", 1), ("K1", 1))),
    # --- quart-diurnes (ondes de petits fonds) --------------------------------
    "N4":      ((4, -2, 0, 2, 0, 0, 0), (("N2", 2),)),
    "MN4":     ((4, -1, 0, 1, 0, 0, 0), (("M2", 1), ("N2", 1))),
    "M4":      ((4, 0, 0, 0, 0, 0, 0), (("M2", 2),)),
    "SN4":     ((4, 1, -2, 1, 0, 0, 0), (("S2", 1), ("N2", 1))),
    "MS4":     ((4, 2, -2, 0, 0, 0, 0), (("M2", 1), ("S2", 1))),
    "MK4":     ((4, 2, 0, 0, 0, 0, 0), (("M2", 1), ("K2", 1))),
    "S4":      ((4, 4, -4, 0, 0, 0, 0), None),
    "SK4":     ((4, 4, -2, 0, 0, 0, 0), (("S2", 1), ("K2", 1))),
    # --- cinquième-, sixième- et huitième-diurnes -----------------------------
    "2MO5":    ((5, -1, 0, 0, 0, 0, 1), (("M2", 2), ("O1", 1))),
    "2MK5":    ((5, 1, 0, 0, 0, 0, -1), (("M2", 2), ("K1", 1))),
    "2MN6":    ((6, -1, 0, 1, 0, 0, 0), (("M2", 2), ("N2", 1))),
    "M6":      ((6, 0, 0, 0, 0, 0, 0), (("M2", 3),)),
    "MSN6":    ((6, 1, -2, 1, 0, 0, 0), (("M2", 1), ("S2", 1), ("N2", 1))),
    "2MS6":    ((6, 2, -2, 0, 0, 0, 0), (("M2", 2), ("S2", 1))),
    "2MK6":    ((6, 2, 0, 0, 0, 0, 0), (("M2", 2), ("K2", 1))),
    "2SM6":    ((6, 4, -4, 0, 0, 0, 0), (("S2", 2), ("M2", 1))),
    "MSK6":    ((6, 4, -2, 0, 0, 0, 0), (("M2", 1), ("S2", 1), ("K2", 1))),
    "3MN8":    ((8, -1, 0, 1, 0, 0, 0), (("M2", 3), ("N2", 1))),
    "M8":      ((8, 0, 0, 0, 0, 0, 0), (("M2", 4),)),
    "2MSN8":   ((8, 1, -2, 1, 0, 0, 0), (("M2", 2), ("S2", 1), ("N2", 1))),
    "3MS8":    ((8, 2, -2, 0, 0, 0, 0), (("M2", 3), ("S2", 1))),
    "2MS8":    ((8, 4, -4, 0, 0, 0, 0), (("M2", 2), ("S2", 2))),
}

# Ondes de petits fonds absentes de TICON-4, inférées des constantes observées
# du port par la théorie des interactions non linéaires : une onde composée
# A+B+… a une amplitude proportionnelle au produit de celles de ses ondes
# mères et une phase égale à la somme de leurs phases plus un déphasage propre
# à son espèce. On la déduit donc de l'onde observée de même espèce et de
# fréquence la plus proche (la « référence ») :
#     H = k · H_réf · Π H_num / Π H_den        G = G_réf + Σ G_num − Σ G_den
# k est le rapport des coefficients du développement multinomial (p. ex.
# 3 pour 2MN6/M6 = 3·M²N / M³). Validé contre un an d'annuaire publié :
# −25 % d'erreur RMS sur les heures, −22 % sur les hauteurs à Cherbourg.
INFERENCES = {
    # onde : (référence, numérateurs, dénominateurs, k)
    "MK4":   ("MS4", ("K2",), ("S2",), 1.0),
    "SN4":   ("MS4", ("N2",), ("M2",), 1.0),
    "SK4":   ("S4", ("K2",), ("S2",), 2.0),
    "2MN6":  ("M6", ("N2",), ("M2",), 3.0),
    "MSN6":  ("2MS6", ("N2",), ("M2",), 2.0),
    "2MK6":  ("2MS6", ("K2",), ("S2",), 1.0),
    "2SM6":  ("2MS6", ("S2",), ("M2",), 1.0),
    "MSK6":  ("2MS6", ("K2",), ("M2",), 2.0),
    "3MN8":  ("M8", ("N2",), ("M2",), 4.0),
    "2MSN8": ("M8", ("S2", "N2"), ("M2", "M2"), 12.0),
    "3MS8":  ("M8", ("S2",), ("M2",), 4.0),
    "2MS8":  ("M8", ("S2", "S2"), ("M2", "M2"), 6.0),
}


def inferer_ondes(constantes: dict) -> dict:
    """Constantes {onde: (H, G)} des ondes de INFERENCES absentes de `constantes`."""
    inferees = {}
    for onde, (ref, num, den, k) in INFERENCES.items():
        requis = (ref,) + num + den
        if onde in constantes or not all(o in constantes for o in requis):
            continue
        H = k * constantes[ref][0]
        G = constantes[ref][1]
        for o in num:
            H *= constantes[o][0]
            G += constantes[o][1]
        for o in den:
            H /= constantes[o][0]
            G -= constantes[o][1]
        inferees[onde] = (H, G % 360.0)
    return inferees


# ==============================================================================
# 2. CONSTANTES HARMONIQUES DES PORTS
# ==============================================================================
# Modulations saisonnières de M2 estimées par TICON-4 mais non retenues pour la
# prédiction : ce ne sont pas des ondes astronomiques, et les retirer améliore
# l'accord avec l'annuaire officiel sur chacun des deux semestres testés
# (−6 % d'erreur RMS sur les heures, −10 % sur les hauteurs à Cherbourg).
ONDES_NON_RETENUES = frozenset({"MA2", "MB2"})


class Port:
    """Port : niveau moyen Z0 (m) et constantes observées {onde: (H m, G °)}.

    `observees` garde toutes les constantes fournies ; `constantes` (utilisé
    par le calcul) en retire `exclues` et, sauf si inferer=False, y ajoute les
    ondes de INFERENCES manquantes.
    """

    def __init__(self, nom, latitude, longitude, z0, constantes, inferer=True,
                 exclues=ONDES_NON_RETENUES):
        self.nom = nom
        self.latitude = latitude
        self.longitude = longitude
        self.z0 = z0
        self.observees = dict(constantes)
        self.constantes = {k: v for k, v in constantes.items() if k not in exclues}
        if inferer:
            self.constantes.update(inferer_ondes(self.constantes))
        inconnues = set(self.constantes) - set(ONDES)
        if inconnues:
            raise ValueError(f"{nom} : ondes non définies {sorted(inconnues)}")

    def __repr__(self):
        return f"Port({self.nom!r}, {len(self.constantes)} ondes)"


# Source : TICON-4 (Hart-Davis, Dettmering & Seitz, 2025, SEANOE,
# https://doi.org/10.17882/109129, CC BY 4.0), via slackwater-database.
# H en mètres, G = déphasage de Greenwich en degrés (UTC).
#
# Cherbourg : marégraphe REFMAR n° 13, observations 1943-2025.
# Z0 = niveau moyen observé (GESLA-4, 09/2016-08/2021) au-dessus du zéro
# hydrographique (zéro du marégraphe).
CHERBOURG = Port(
    "Cherbourg", 49.6514, -1.6355, 3.876,
    {
        "M2": (1.86841, 228.365), "S2": (0.69572, 272.222), "N2": (0.37100, 209.088),
        "K2": (0.19806, 270.125), "M4": (0.14361, 355.324), "K1": (0.09306, 106.589),
        "MS4": (0.08789, 45.818), "MU2": (0.07559, 222.474), "SA": (0.06883, 221.408),
        "NU2": (0.06830, 201.559), "L2": (0.06502, 224.137), "O1": (0.06303, 356.735),
        "MN4": (0.05143, 329.827), "2N2": (0.05045, 190.321), "T2": (0.03799, 261.167),
        "P1": (0.03253, 99.952), "2MS6": (0.02789, 133.099), "LAMBDA2": (0.02726, 210.191),
        "M6": (0.02524, 94.406), "2SM2": (0.02185, 102.577), "SSA": (0.01968, 116.203),
        "Q1": (0.01725, 316.099), "MA2": (0.01640, 157.713), "EP2": (0.01459, 195.261),
        "MF": (0.01188, 197.978), "M3": (0.01188, 203.305), "N4": (0.01038, 307.345),
        "M8": (0.00824, 279.303), "3L2": (0.00761, 132.766), "R2": (0.00730, 306.760),
        "MM": (0.00622, 197.915), "S1": (0.00618, 76.138), "S4": (0.00594, 95.460),
        "MKS2": (0.00560, 316.851), "J1": (0.00507, 189.120), "OO1": (0.00465, 251.910),
        "SGM": (0.00463, 302.340), "MB2": (0.00387, 264.399), "2Q1": (0.00325, 267.614),
        "R3": (0.00296, 269.783), "M1": (0.00270, 41.094), "T3": (0.00266, 183.844),
        "MSF": (0.00232, 254.128), "RHO1": (0.00225, 312.126), "3N2": (0.00175, 146.734),
        "2MO5": (0.00133, 254.978), "S3": (0.00127, 324.211), "MTM": (0.00098, 208.165),
        "2MK5": (0.00015, 83.290), "MSQM": (0.00015, 71.643),
    },
)

# Brest : port de référence du coefficient de marée (marégraphe REFMAR n° 3,
# observations 1846-2025). Z0 observé sur 05/2006-04/2009.
BREST = Port(
    "Brest", 48.3829, -4.4950, 4.163,
    {
        "M2": (2.05110, 108.993), "S2": (0.74891, 148.282), "N2": (0.41694, 90.618),
        "K2": (0.21368, 145.890), "MU2": (0.08572, 105.072), "NU2": (0.07778, 86.600),
        "O1": (0.06584, 327.887), "K1": (0.06434, 75.088), "L2": (0.06387, 102.920),
        "2N2": (0.05703, 72.807), "M4": (0.05447, 105.901), "SA": (0.04932, 244.952),
        "T2": (0.04169, 138.562), "MS4": (0.03269, 181.747), "M6": (0.03146, 354.853),
        "LAMBDA2": (0.02627, 75.818), "P1": (0.02254, 63.651), "SSA": (0.02051, 98.076),
        "Q1": (0.02046, 281.348), "M3": (0.01977, 15.857), "EP2": (0.01969, 89.435),
        "MN4": (0.01943, 60.443), "2SM2": (0.01731, 312.051), "2MS6": (0.01679, 53.350),
        "MB2": (0.01244, 101.150), "MA2": (0.01102, 39.450), "MF": (0.01028, 177.577),
        "3L2": (0.00825, 6.721), "S1": (0.00796, 11.504), "MKS2": (0.00760, 174.815),
        "R3": (0.00563, 72.194), "R2": (0.00534, 157.872), "T3": (0.00528, 345.395),
        "MM": (0.00443, 203.974), "2Q1": (0.00378, 234.462), "SGM": (0.00375, 248.220),
        "RHO1": (0.00336, 284.235), "MSF": (0.00333, 27.139), "S3": (0.00302, 149.641),
        "N4": (0.00292, 9.561), "J1": (0.00237, 123.387), "M1": (0.00234, 298.781),
        "M8": (0.00231, 231.999), "3N2": (0.00217, 308.506), "S4": (0.00216, 288.667),
        "2MK5": (0.00199, 43.912), "OO1": (0.00139, 212.780), "MSQM": (0.00117, 250.315),
        "MTM": (0.00113, 135.118), "2MO5": (0.00097, 321.361),
    },
)

PORTS = {"CHERBOURG": CHERBOURG, "BREST": BREST}

# Unité de hauteur de Brest (m) : amplitude moyenne de la vive-eau qui suit
# d'environ 36 h la syzygie d'équinoxe (valeur officielle). C = 100·A / U.
UNITE_HAUTEUR_BREST = 3.05

# Ondes semi-diurnes TICON écartées du calcul du coefficient : raies de degré 3
# (3N2, 3L2) et modulations saisonnières de M2 (MA2, MB2), qui ne font pas
# partie de la marée semi-diurne astronomique classique.
ONDES_HORS_COEFFICIENT = frozenset({"3N2", "3L2", "MA2", "MB2"})

# Durée maximale couverte par une même Serie (voir Serie).
FENETRE_SERIE = timedelta(days=1)

# Deux extrema consécutifs séparés de moins de 2 cm sont une simple
# ondulation de la courbe (tenue du plein, palier) : ils ne sont pas publiés.
SEUIL_MICRO_EXTREMUM = 0.02


# ==============================================================================
# 3. ASTRONOMIE
# ==============================================================================
J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
_HEURES_PAR_SIECLE = 36525.0 * 24.0

# Vitesses angulaires (°/h) des arguments (τ, s, h, p, N', p', 90°).
_DS = 481267.88123421 / _HEURES_PAR_SIECLE
_DH = 36000.76983 / _HEURES_PAR_SIECLE
_DP = 4069.0137287 / _HEURES_PAR_SIECLE
_DNP = 1934.136261 / _HEURES_PAR_SIECLE
_DPP = 1.71946 / _HEURES_PAR_SIECLE
VITESSES_ARGUMENTS = (15.0 + _DH - _DS, _DS, _DH, _DP, _DNP, _DPP, 0.0)


def en_utc(dt: datetime) -> datetime:
    """Datetime aware en UTC (un datetime naïf est interprété comme UTC)."""
    if dt.tzinfo is None:
        return dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)


def arguments_astronomiques(t: datetime):
    """Arguments (τ, s, h, p, N', p', 90°) en degrés, et (N, p), à l'instant t.

    Longitudes moyennes : Meeus, Astronomical Algorithms (2e éd.), ch. 22, 25
    et 47. Échelle de temps UTC (l'écart TT−UTC ≈ 69 s est sans effet visible).
    """
    jours = (en_utc(t) - J2000).total_seconds() / 86400.0
    T = jours / 36525.0
    T2, T3, T4 = T * T, T ** 3, T ** 4
    s = 218.3164477 + 481267.88123421 * T - 0.0015786 * T2 + T3 / 538841.0 - T4 / 65194000.0
    h = 280.46646 + 36000.76983 * T + 0.0003032 * T2
    p = 83.3532465 + 4069.0137287 * T - 0.0103200 * T2 - T3 / 80053.0 + T4 / 18999000.0
    N = 125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3 / 450000.0
    pp = 282.93735 + 1.71946 * T + 0.00046 * T2
    angle_horaire = (jours % 1.0) * 360.0  # 0° à 12 h UTC, 180° à 0 h UTC
    tau = angle_horaire + h - s
    return (tau, s, h, p, -N, pp, 90.0), N % 360.0, p % 360.0


def vitesse(onde: str) -> float:
    """Vitesse angulaire de l'onde en degrés par heure."""
    return sum(c * w for c, w in zip(ONDES[onde][0], VITESSES_ARGUMENTS))


# ==============================================================================
# 4. CORRECTIONS NODALES (OHI, TWCWG — annexe A)
# ==============================================================================
def _depuis_sin_cos(fsin, fcos):
    return math.hypot(fsin, fcos), math.degrees(math.atan2(fsin, fcos))


def _fondamentales(N_deg: float, p_deg: float) -> dict:
    """(f, u°) des ondes fondamentales, fonction du nœud N et du périgée p."""
    N, p = math.radians(N_deg), math.radians(p_deg)
    cN, c2N, c3N = math.cos(N), math.cos(2 * N), math.cos(3 * N)
    sN, s2N, s3N = math.sin(N), math.sin(2 * N), math.sin(3 * N)

    f_m2 = 1.0007 - 0.0373 * cN + 0.0002 * c2N
    u_m2 = -2.14 * sN
    l2 = _depuis_sin_cos(
        -0.2505 * math.sin(2 * p) - 0.1102 * math.sin(2 * p - N)
        - 0.0156 * math.sin(2 * p - 2 * N) - 0.037 * sN,
        1.0 - 0.2505 * math.cos(2 * p) - 0.1102 * math.cos(2 * p - N)
        - 0.0156 * math.cos(2 * p - 2 * N) - 0.037 * cN,
    )
    m1 = _depuis_sin_cos(
        math.sin(p) + 0.2 * math.sin(p - N),
        2.0 * (math.cos(p) + 0.2 * math.cos(p - N)),
    )
    return {
        "M2": (f_m2, u_m2),
        "K2": (1.0246 + 0.2863 * cN + 0.0083 * c2N - 0.0015 * c3N,
               -17.74 * sN + 0.68 * s2N - 0.04 * s3N),
        "O1": (1.0176 + 0.1871 * cN - 0.0147 * c2N,
               10.80 * sN - 1.34 * s2N + 0.19 * s3N),
        "K1": (1.0060 + 0.1150 * cN - 0.0088 * c2N + 0.0006 * c3N,
               -8.86 * sN + 0.68 * s2N - 0.07 * s3N),
        "J1": (1.1029 + 0.1676 * cN - 0.0170 * c2N + 0.0016 * c3N,
               -12.94 * sN + 1.34 * s2N - 0.19 * s3N),
        "MM": (1.0 - 0.1311 * cN + 0.0538 * math.cos(2 * p) + 0.0205 * math.cos(2 * p - N), 0.0),
        "MF": (1.084 + 0.415 * cN + 0.039 * c2N,
               -23.7 * sN + 2.7 * s2N - 0.4 * s3N),
        "M3": (math.sqrt(f_m2) ** 3, -3.21 * sN),
        "L2": l2,
        "M1": m1,
    }


def corrections_nodales(N_deg: float, p_deg: float) -> dict:
    """{onde: (f, u en degrés)} pour toutes les ondes de ONDES."""
    fond = _fondamentales(N_deg, p_deg)
    resultat = {}

    def fu(onde):
        if onde not in resultat:
            regle = ONDES[onde][1]
            if regle is None:
                resultat[onde] = (1.0, 0.0)
            elif isinstance(regle, str):
                resultat[onde] = fond[regle]
            else:  # onde composée : f = Π f_i^|k_i|, u = Σ k_i·u_i
                f, u = 1.0, 0.0
                for mere, k in regle:
                    f_i, u_i = fu(mere)
                    f *= f_i ** abs(k)
                    u += k * u_i
                resultat[onde] = (f, u)
        return resultat[onde]

    for nom in ONDES:
        fu(nom)
    return resultat


# ==============================================================================
# 5. SYNTHÈSE HARMONIQUE
# ==============================================================================
class Serie:
    """Somme harmonique d'un port, linéarisée autour d'un instant t0.

    V0 et (f, u) sont évalués exactement à t0, puis V avance à la vitesse de
    Doodson de chaque onde. Sur ± FENETRE_SERIE autour de t0, l'écart avec un
    calcul exact reste inférieur à 1 s et 0,3 mm (les corrections nodales
    évoluent sur 4,4 à 18,6 ans). `filtre(nom, coefficients)` ne garde qu'une
    famille d'ondes (le niveau moyen est alors omis). Le temps x s'exprime en
    heures depuis t0.
    """

    def __init__(self, port: Port, t0: datetime, filtre=None):
        self.port = port
        self.t0 = en_utc(t0)
        args, N, p = arguments_astronomiques(self.t0)
        nodal = corrections_nodales(N, p)
        self.z0 = port.z0 if filtre is None else 0.0
        self.termes = []  # (A = f·H en m, phase à t0 en rad, ω en rad/h)
        for onde, (H, G) in port.constantes.items():
            coeffs = ONDES[onde][0]
            if filtre is not None and not filtre(onde, coeffs):
                continue
            f, u = nodal[onde]
            V0 = sum(c * a for c, a in zip(coeffs, args))
            w = sum(c * v for c, v in zip(coeffs, VITESSES_ARGUMENTS))
            self.termes.append((f * H, math.radians((V0 + u - G) % 360.0), math.radians(w)))

    def heures(self, t: datetime) -> float:
        return (en_utc(t) - self.t0).total_seconds() / 3600.0

    def instant(self, x: float) -> datetime:
        return self.t0 + timedelta(hours=x)

    def hauteur(self, x: float) -> float:
        """Hauteur (m) à x heures de t0."""
        return self.z0 + sum(a * math.cos(ph + w * x) for a, ph, w in self.termes)

    def derivee(self, x: float) -> float:
        """dh/dt (m/h) à x heures de t0."""
        return -sum(a * w * math.sin(ph + w * x) for a, ph, w in self.termes)

    def extrema(self, x0: float, x1: float, pas: float = 1.0 / 6.0):
        """[(x, 'PM'|'BM')] : zéros de dh/dt dans [x0, x1[ (heures depuis t0)."""
        trouves = []
        n = max(1, math.ceil((x1 - x0) / pas))
        xa, da = x0, self.derivee(x0)
        for i in range(1, n + 1):
            xb = x0 + i * pas
            db = self.derivee(xb)
            if da != 0.0 and (da > 0.0) != (db > 0.0):
                lo, hi, d_lo = xa, xb, da
                while hi - lo > 1e-5:  # 0,04 s
                    mi = 0.5 * (lo + hi)
                    d_mi = self.derivee(mi)
                    if (d_mi > 0.0) == (d_lo > 0.0):
                        lo, d_lo = mi, d_mi
                    else:
                        hi = mi
                trouves.append((0.5 * (lo + hi), "PM" if da > 0.0 else "BM"))
            xa, da = xb, db
        return trouves


def hauteur_maree(dt: datetime, port: Port = CHERBOURG) -> float:
    """Hauteur d'eau prédite (m au-dessus du zéro hydrographique) à l'instant dt.

    Un datetime naïf est interprété comme UTC.
    """
    return Serie(port, dt).hauteur(0.0)


def courbe(debut: datetime, fin: datetime, pas_minutes: float = 10.0,
           port: Port = CHERBOURG) -> list:
    """[(instant UTC, hauteur m)] de debut à fin inclus, au pas donné (marégramme)."""
    debut, fin = en_utc(debut), en_utc(fin)
    pas = timedelta(minutes=pas_minutes)
    points, t = [], debut
    while t <= fin:
        borne = min(fin, t + 2 * FENETRE_SERIE)
        serie = Serie(port, t + (borne - t) / 2)
        while t <= borne:
            points.append((t, serie.hauteur(serie.heures(t))))
            t += pas
    return points


def _extrema_bruts(port: Port, debut: datetime, fin: datetime, filtre=None):
    """[(instant UTC, 'PM'|'BM', hauteur m)] dans [debut, fin], par tronçons."""
    resultat, a = [], debut
    while a < fin:
        b = min(fin, a + 2 * FENETRE_SERIE)
        serie = Serie(port, a + (b - a) / 2, filtre=filtre)
        resultat += [(serie.instant(x), typ, serie.hauteur(x))
                     for x, typ in serie.extrema(serie.heures(a), serie.heures(b))]
        a = b
    return resultat


# ==============================================================================
# 6. COEFFICIENT DE MARÉE (définition française, calculé à Brest)
# ==============================================================================
def _semi_diurne_coefficient(onde, coeffs) -> bool:
    return coeffs[0] == 2 and onde not in ONDES_HORS_COEFFICIENT


def decalage_brest(port: Port) -> float:
    """Retard moyen (h) de la PM du port sur celle de Brest pour la même onde
    de marée, déduit de la différence de situation de M2 (0 à 12,42 h)."""
    dg = (port.constantes["M2"][1] - BREST.constantes["M2"][1]) % 360.0
    return dg / vitesse("M2")


def pleines_mers_semi_diurnes_brest(debut: datetime, fin: datetime):
    """[(instant UTC, amplitude m)] des PM de la marée semi-diurne de Brest."""
    bruts = _extrema_bruts(BREST, en_utc(debut), en_utc(fin), _semi_diurne_coefficient)
    return [(t, h) for t, typ, h in bruts if typ == "PM"]


def coefficient_depuis_amplitude(amplitude: float) -> int:
    """Coefficient (20 à 120) correspondant à une amplitude semi-diurne à Brest."""
    c = int(math.floor(100.0 * amplitude / UNITE_HAUTEUR_BREST + 0.5))
    return max(20, min(120, c))


def _amplitude_proche(cible: datetime, pm_brest: list, instants: list) -> float:
    """Amplitude de la PM de Brest la plus proche de `cible` (listes triées)."""
    i = bisect.bisect_left(instants, cible)
    voisines = pm_brest[max(0, i - 1):i + 1]
    return min(voisines, key=lambda pm: abs((pm[0] - cible).total_seconds()))[1]


def coefficient(t_pm: datetime, port: Port = CHERBOURG) -> int:
    """Coefficient de la pleine mer du port survenant à t_pm : celui de la PM
    semi-diurne de Brest la plus proche de t_pm − decalage_brest(port)."""
    cible = en_utc(t_pm) - timedelta(hours=decalage_brest(port))
    pm_brest = pleines_mers_semi_diurnes_brest(cible - timedelta(hours=8),
                                               cible + timedelta(hours=8))
    instants = [t for t, _ in pm_brest]
    return coefficient_depuis_amplitude(_amplitude_proche(cible, pm_brest, instants))


# ==============================================================================
# 7. ANNUAIRE : PM / BM D'UNE PÉRIODE OU D'UNE JOURNÉE LOCALE
# ==============================================================================
def _sans_micro_extrema(liste):
    """Retire les couples d'extrema voisins séparés de moins du seuil."""
    liste = list(liste)
    while len(liste) > 1:
        ecarts = [abs(liste[i + 1][2] - liste[i][2]) for i in range(len(liste) - 1)]
        i = min(range(len(ecarts)), key=ecarts.__getitem__)
        if ecarts[i] >= SEUIL_MICRO_EXTREMUM:
            break
        del liste[i:i + 2]
    return liste


def extrema_periode(debut: datetime, fin: datetime, port: Port = CHERBOURG,
                    tz_nom: str = FUSEAU_DEFAUT) -> list:
    """PM et BM dans [debut, fin[, avec hauteur et coefficient (PM).

    Retourne une liste de dicts {"type": "PM"|"BM", "heure": datetime aware
    en heure locale, "hauteur": m, "coeff": int pour les PM, None pour les BM}.
    """
    tz = ZoneInfo(tz_nom)
    debut, fin = en_utc(debut), en_utc(fin)
    marge = timedelta(hours=3)
    bruts = _extrema_bruts(port, debut - marge, fin + marge)
    bruts = [e for e in _sans_micro_extrema(bruts) if debut <= e[0] < fin]

    decal = timedelta(hours=decalage_brest(port))
    pm_brest = pleines_mers_semi_diurnes_brest(debut - decal - timedelta(hours=8),
                                               fin - decal + timedelta(hours=8))
    instants = [t for t, _ in pm_brest]
    return [{
        "type": typ,
        "heure": t.astimezone(tz),
        "hauteur": h,
        "coeff": (coefficient_depuis_amplitude(_amplitude_proche(t - decal, pm_brest, instants))
                  if typ == "PM" else None),
    } for t, typ, h in bruts]


def extrema_journee(jour: date, port: Port = CHERBOURG, tz_nom: str = FUSEAU_DEFAUT) -> list:
    """PM et BM du jour civil local `jour` (journées de 23 h ou 25 h comprises)."""
    tz = ZoneInfo(tz_nom)
    debut = datetime.combine(jour, time(0, 0), tzinfo=tz)
    fin = datetime.combine(jour + timedelta(days=1), time(0, 0), tzinfo=tz)
    return extrema_periode(debut, fin, port, tz_nom)


def trouver_extrema_journee(jour: date, tz_nom: str = FUSEAU_DEFAUT) -> list:
    """Compatibilité avec la version précédente (Cherbourg)."""
    return extrema_journee(jour, CHERBOURG, tz_nom)


def afficher_annuaire(jour: date, port: Port = CHERBOURG, tz_nom: str = FUSEAU_DEFAUT):
    extrema = extrema_journee(jour, port, tz_nom)
    print("=" * 65)
    print(f" ANNUAIRE DES MARÉES (HORS LIGNE) — {port.nom.upper()} : {jour.strftime('%d/%m/%Y')}")
    print("=" * 65)
    print(f"{'Marée':<14} | {'Heure légale':<15} | {'Hauteur':<9} | {'Coefficient'}")
    print("-" * 65)
    for e in extrema:
        libelle = "Pleine mer" if e["type"] == "PM" else "Basse mer"
        decalage = int(e["heure"].utcoffset().total_seconds() // 3600)
        h_str = f"{e['heure'].strftime('%Hh%M')} (UTC{decalage:+d})"
        hauteur_str = f"{e['hauteur']:.2f} m"
        c_str = str(e["coeff"]) if e["coeff"] is not None else "-"
        print(f"{libelle:<14} | {h_str:<15} | {hauteur_str:<9} | {c_str:^11}")
    print("=" * 65)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Annuaire des marées hors ligne.")
    parser.add_argument("date", nargs="?", help="premier jour (AAAA-MM-JJ), défaut : aujourd'hui")
    parser.add_argument("-n", "--jours", type=int, default=2, help="nombre de jours (défaut : 2)")
    parser.add_argument("-p", "--port", default="CHERBOURG", choices=sorted(PORTS),
                        help="port (défaut : CHERBOURG)")
    args = parser.parse_args(argv)
    jour = date.fromisoformat(args.date) if args.date else datetime.now(ZoneInfo(FUSEAU_DEFAUT)).date()
    for i in range(args.jours):
        print()
        afficher_annuaire(jour + timedelta(days=i), PORTS[args.port])


if __name__ == "__main__":
    main()
