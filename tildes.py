#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Océanos — calculateur harmonique des marées, 100 % hors ligne.

Calcule les pleines mers (PM), les basses mers (BM), les hauteurs d'eau
au-dessus du zéro hydrographique et les coefficients de marée pour
Cherbourg-en-Cotentin (et Brest). Aucun appel réseau, aucune dépendance
externe : bibliothèque standard Python (>= 3.9) uniquement.

Modèle
------
    h(t) = Z0 + Σ f·H·cos(V(t) + u − G)

* V(t) : argument astronomique de chaque onde (nombres de Doodson, convention
  de Schureman), calculé à l'instant t à partir des longitudes moyennes de la
  Lune, du Soleil, du périgée lunaire, du nœud ascendant et du périhélie
  (polynômes de Meeus, « Astronomical Algorithms »).
* f, u : corrections nodales de l'OHI (liste des constituantes du TWCWG,
  annexe A) ; les ondes composées héritent de leurs ondes mères (annexe B).
  Les écarts réels à ces formules sont pris en compte par des satellites
  nodaux (± N') ajustés librement pour M2, N2 et K2.
* H, G : amplitude et situation (déphasage de Greenwich, temps UTC) de
  chaque onde, issues de notre analyse harmonique par moindres carrés des
  hauteurs d'eau observées par les marégraphes REFMAR de Cherbourg et de
  Brest (données validées au pas de 10 min, voir outils/ajuster.py) : ondes
  astronomiques et de petits fonds, plus les raies de marée trouvées dans le
  spectre du résidu. Aucune prédiction ni constante harmonique du SHOM n'est
  utilisée. Z0 est le niveau moyen officiel du SHOM, celui de l'annuaire.
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
from datetime import date, datetime, time, timedelta, timezone, tzinfo
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

FUSEAU_DEFAUT = "Europe/Paris"


class HeureLegaleFrance(tzinfo):
    """Heure légale de la France métropolitaine, sans base de fuseaux horaires.

    Règle de l'Union européenne, en vigueur depuis 1996 : heure d'été (CEST,
    UTC+2) du dernier dimanche de mars à 1 h UTC au dernier dimanche d'octobre
    à 1 h UTC, heure d'hiver (CET, UTC+1) le reste de l'année. Sert de repli
    quand Python n'a pas la base IANA des fuseaux (Windows sans le paquet
    tzdata). Les heures locales inexistantes ou répétées suivent le PEP 495
    (attribut fold), comme zoneinfo.
    """

    HIVER, ETE = timedelta(hours=1), timedelta(hours=2)

    @staticmethod
    def bascules(annee: int):
        """Instants UTC (naïfs) du passage à l'heure d'été et à l'heure d'hiver."""
        def dernier_dimanche(mois):
            fin_mois = date(annee, mois, 31)
            return datetime.combine(fin_mois - timedelta(days=(fin_mois.weekday() + 1) % 7), time(1))
        return dernier_dimanche(3), dernier_dimanche(10)

    def utcoffset(self, dt):
        if dt is None:
            return self.HIVER
        ete, hiver = self.bascules(dt.year)
        local = dt.replace(tzinfo=None, fold=0)
        h = timedelta(hours=1)
        if ete + 2 * h <= local < hiver + h:
            return self.ETE
        if ete + h <= local < ete + 2 * h:      # 2 h-3 h le jour du passage à l'été : inexistante
            return self.ETE if dt.fold else self.HIVER
        if hiver + h <= local < hiver + 2 * h:  # 2 h-3 h le jour du retour : répétée
            return self.HIVER if dt.fold else self.ETE
        return self.HIVER

    def dst(self, dt):
        return self.utcoffset(dt) - self.HIVER

    def tzname(self, dt):
        return "CEST" if self.dst(dt) else "CET"

    def fromutc(self, dt):
        if dt.tzinfo is not self:
            raise ValueError("fromutc : dt.tzinfo n'est pas self")
        utc = dt.replace(tzinfo=None)
        ete, hiver = self.bascules(utc.year)
        if ete <= utc < hiver:
            return (utc + self.ETE).replace(tzinfo=self)
        repetee = hiver <= utc < hiver + timedelta(hours=1)
        return (utc + self.HIVER).replace(tzinfo=self, fold=int(repetee))

    def __repr__(self):
        return "HeureLegaleFrance()"


HEURE_LEGALE_FRANCE = HeureLegaleFrance()


def fuseau(nom: str = FUSEAU_DEFAUT) -> tzinfo:
    """Fuseau horaire `nom` (base IANA de zoneinfo).

    Sans base IANA (cas de Python sous Windows sans le paquet tzdata),
    « Europe/Paris » est remplacé par HEURE_LEGALE_FRANCE ; tout autre fuseau
    demande d'installer tzdata.
    """
    try:
        return ZoneInfo(nom)
    except ZoneInfoNotFoundError:
        if nom == "Europe/Paris":
            return HEURE_LEGALE_FRANCE
        raise ZoneInfoNotFoundError(
            f"fuseau {nom!r} introuvable : installez le paquet tzdata "
            "(py -m pip install tzdata)") from None

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
    "MSM":     ((0, 1, -2, 1, 0, 0, 0), (("MM", 1),)),
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
    "OQ2":     ((2, -3, 0, 1, 0, 0, 2), (("O1", 1), ("Q1", 1))),
    "2NS2":    ((2, -4, 2, 2, 0, 0, 0), (("N2", 2), ("S2", -1))),
    "MSK2":    ((2, 0, -2, 0, 0, 0, 0), (("M2", 1), ("S2", 1), ("K2", -1))),
    "2SK2":    ((2, 2, -4, 0, 0, 0, 0), (("S2", 2), ("K2", -1))),
    "MSN2":    ((2, 3, -2, -1, 0, 0, 0), (("M2", 1), ("S2", 1), ("N2", -1))),
    # --- tiers-diurnes --------------------------------------------------------
    "MO3":     ((3, -1, 0, 0, 0, 0, 1), (("M2", 1), ("O1", 1))),
    "NO3":     ((3, -2, 0, 1, 0, 0, 1), (("N2", 1), ("O1", 1))),
    "SO3":     ((3, 1, -2, 0, 0, 0, 1), (("S2", 1), ("O1", 1))),
    "MK3":     ((3, 1, 0, 0, 0, 0, -1), (("M2", 1), ("K1", 1))),
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
    "2MNS4":   ((4, -3, 2, 1, 0, 0, 0), (("M2", 2), ("N2", 1), ("S2", -1))),
    "3MS4":    ((4, -2, 2, 0, 0, 0, 0), (("M2", 3), ("S2", -1))),
    "MNU4":    ((4, -1, 2, -1, 0, 0, 0), (("M2", 1), ("NU2", 1))),
    "ML4":     ((4, 1, 0, -1, 0, 0, -2), (("M2", 1), ("L2", 1))),
    "2MSN4":   ((4, 3, -2, -1, 0, 0, 0), (("M2", 2), ("S2", 1), ("N2", -1))),
    # --- cinquième-, sixième- et huitième-diurnes -----------------------------
    "2MO5":    ((5, -1, 0, 0, 0, 0, 1), (("M2", 2), ("O1", 1))),
    "2MK5":    ((5, 1, 0, 0, 0, 0, -1), (("M2", 2), ("K1", 1))),
    "2NM6":    ((6, -2, 0, 2, 0, 0, 0), (("N2", 2), ("M2", 1))),
    "4MS6":    ((6, -2, 2, 0, 0, 0, 0), (("M2", 4), ("S2", -1))),
    "2MN6":    ((6, -1, 0, 1, 0, 0, 0), (("M2", 2), ("N2", 1))),
    "M6":      ((6, 0, 0, 0, 0, 0, 0), (("M2", 3),)),
    "MSN6":    ((6, 1, -2, 1, 0, 0, 0), (("M2", 1), ("S2", 1), ("N2", 1))),
    "2MS6":    ((6, 2, -2, 0, 0, 0, 0), (("M2", 2), ("S2", 1))),
    "2MK6":    ((6, 2, 0, 0, 0, 0, 0), (("M2", 2), ("K2", 1))),
    "2SM6":    ((6, 4, -4, 0, 0, 0, 0), (("S2", 2), ("M2", 1))),
    "MSK6":    ((6, 4, -2, 0, 0, 0, 0), (("M2", 1), ("S2", 1), ("K2", 1))),
    "3MK7":    ((7, 1, 0, 0, 0, 0, -1), (("M2", 3), ("K1", 1))),
    "2(MN)8":  ((8, -2, 0, 2, 0, 0, 0), (("M2", 2), ("N2", 2))),
    "3MN8":    ((8, -1, 0, 1, 0, 0, 0), (("M2", 3), ("N2", 1))),
    "M8":      ((8, 0, 0, 0, 0, 0, 0), (("M2", 4),)),
    "2MSN8":   ((8, 1, -2, 1, 0, 0, 0), (("M2", 2), ("S2", 1), ("N2", 1))),
    "3MS8":    ((8, 2, -2, 0, 0, 0, 0), (("M2", 3), ("S2", 1))),
    "2MS8":    ((8, 4, -4, 0, 0, 0, 0), (("M2", 2), ("S2", 2))),
    "3MK8":    ((8, 2, 0, 0, 0, 0, 0), (("M2", 3), ("K2", 1))),
    "2MSK8":   ((8, 4, -2, 0, 0, 0, 0), (("M2", 2), ("S2", 1), ("K2", 1))),
    # --- dixième- et douzième-diurnes -----------------------------------------
    "4MN10":   ((10, -1, 0, 1, 0, 0, 0), (("M2", 4), ("N2", 1))),
    "M10":     ((10, 0, 0, 0, 0, 0, 0), (("M2", 5),)),
    "4MS10":   ((10, 2, -2, 0, 0, 0, 0), (("M2", 4), ("S2", 1))),
    "3M2S10":  ((10, 4, -4, 0, 0, 0, 0), (("M2", 3), ("S2", 2))),
    "M12":     ((12, 0, 0, 0, 0, 0, 0), (("M2", 6),)),
    "5MS12":   ((12, 2, -2, 0, 0, 0, 0), (("M2", 5), ("S2", 1))),
    "4M2S12":  ((12, 4, -4, 0, 0, 0, 0), (("M2", 4), ("S2", 2))),
}

# Satellites nodaux : raies à ± N' (± 1/18,61 an) de M2, N2 et K2, ajustées
# librement et sans correction nodale. Elles corrigent l'écart entre la
# modulation nodale réelle et celle des formules de l'OHI : à Cherbourg, la
# modulation de M2 est plus faible que la théorie d'équilibre, comme ailleurs
# en Manche et en mer du Nord. Sans elles, l'heure des PM y oscille de ± 2 à
# 3 min au rythme du cycle nodal.
SATELLITES_NODAUX = {
    "N2-N":         ((2, -1, 0, 1, -1, 0, 0), None),
    "N2+N":         ((2, -1, 0, 1, 1, 0, 0), None),
    "M2-N":         ((2, 0, 0, 0, -1, 0, 0), None),
    "M2+N":         ((2, 0, 0, 0, 1, 0, 0), None),
    "K2-N":         ((2, 2, 0, 0, -1, 0, 0), None),
    "K2+N":         ((2, 2, 0, 0, 1, 0, 0), None),
}

# Raies trouvées dans le spectre du résidu des analyses de Cherbourg et de
# Brest (amplitude ≥ 2 mm, voir outils/ajuster.py) : ondes de petits fonds
# d'ordre élevé, modulations saisonnières et satellites d'ondes existantes.
# Nom usuel quand il existe, sinon X(τ,s,h,p). Correction nodale : celle de la
# combinaison d'ondes mères (M2, S2, N2, K2, K1, O1) d'ordre ≤ 5 la plus
# plausible, ou aucune (None) faute de combinaison.
RAIES_EMPIRIQUES = {
    "X(1,-2,0,0)":  ((1, -2, 0, 0, 0, 0, 0), None),
    "X(1,-1,1,0)":  ((1, -1, 1, 0, 0, 0, 0), None),
    "TAU1":         ((1, -1, 2, 0, 0, 0, 3), (("M2", 1), ("S2", -1), ("K1", 1))),
    "X(1,0,-2,1)":  ((1, 0, -2, 1, 0, 0, 1), (("M2", -1), ("S2", 1), ("N2", 1), ("K1", -1))),
    "X(1,0,0,0)":   ((1, 0, 0, 0, 0, 0, 0), None),
    "X(1,1,-3,0)":  ((1, 1, -3, 0, 0, 0, 0), None),
    "X(1,2,0,0)":   ((1, 2, 0, 0, 0, 0, 0), None),
    "X(1,3,-4,0)":  ((1, 3, -4, 0, 0, 0, 1), (("M2", -1), ("S2", 2), ("K1", -1))),
    "SO1":          ((1, 3, -2, 0, 0, 0, 3), (("S2", 1), ("O1", -1))),
    "X(1,5,-4,0)":  ((1, 5, -4, 0, 0, 0, 3), (("M2", -1), ("S2", 2), ("O1", -1))),
    "X(2,-6,6,0)":  ((2, -6, 6, 0, 0, 0, 0), None),
    "X(2,-5,4,1)":  ((2, -5, 4, 1, 0, 0, 0), (("M2", 2), ("S2", -2), ("N2", 1))),
    "X(2,-4,2,0)":  ((2, -4, 2, 0, 0, 0, 0), (("M2", 3), ("S2", -1), ("K2", -1))),
    "3M2S2":        ((2, -4, 4, 0, 0, 0, 0), (("M2", 3), ("S2", -2))),
    "X(2,-3,0,3)":  ((2, -3, 0, 3, 0, 0, 0), (("M2", -2), ("N2", 3))),
    "X(2,-3,2,-1)": ((2, -3, 2, -1, 0, 0, 0), None),
    "X(2,-3,4,-1)": ((2, -3, 4, -1, 0, 0, 0), None),
    "2MK2":         ((2, -2, 0, 0, 0, 0, 0), (("M2", 2), ("K2", -1))),
    "X(2,-2,3,0)":  ((2, -2, 3, 0, 0, 0, 0), None),
    "X(2,-2,4,-2)": ((2, -2, 4, -2, 0, 0, 0), None),
    "X(2,-1,-1,1)": ((2, -1, -1, 1, 0, 0, 0), None),
    "X(2,-1,0,-2)": ((2, -1, 0, -2, 0, 0, 0), None),
    "X(2,-1,0,-1)": ((2, -1, 0, -1, 0, 0, 0), (("M2", 3), ("N2", -1), ("K2", -1))),
    "X(2,-1,1,0)":  ((2, -1, 1, 0, 0, 0, 0), None),
    "X(2,-1,1,1)":  ((2, -1, 1, 1, 0, 0, 0), None),
    "X(2,-1,3,-1)": ((2, -1, 3, -1, 0, 0, 0), None),
    "X(2,0,-2,2)":  ((2, 0, -2, 2, 0, 0, 0), (("M2", -2), ("S2", 1), ("N2", 2))),
    "X(2,1,-2,-1)": ((2, 1, -2, -1, 0, 0, 0), (("M2", 2), ("S2", 1), ("N2", -1), ("K2", -1))),
    "X(2,1,0,1)":   ((2, 1, 0, 1, 0, 0, 0), (("M2", -1), ("N2", 1), ("K2", 1))),
    "X(2,1,2,-1)":  ((2, 1, 2, -1, 0, 0, 0), (("M2", 2), ("S2", -1), ("N2", -1), ("K2", 1))),
    "X(2,2,1,0)":   ((2, 2, 1, 0, 0, 0, 0), None),
    "X(2,3,-2,1)":  ((2, 3, -2, 1, 0, 0, 0), (("M2", -2), ("S2", 1), ("N2", 1), ("K2", 1))),
    "ETA2":         ((2, 3, 0, -1, 0, 0, 0), (("M2", 1), ("N2", -1), ("K2", 1))),
    "X(2,4,-5,0)":  ((2, 4, -5, 0, 0, 0, 0), None),
    "SKM2":         ((2, 4, -2, 0, 0, 0, 0), (("M2", -1), ("S2", 1), ("K2", 1))),
    "X(2,6,-6,0)":  ((2, 6, -6, 0, 0, 0, 0), (("M2", -2), ("S2", 3))),
    "X(2,6,-4,0)":  ((2, 6, -4, 0, 0, 0, 0), (("M2", -2), ("S2", 2), ("K2", 1))),
    "X(3,-1,0,1)":  ((3, -1, 0, 1, 0, 0, 0), None),
    "X(3,2,0,0)":   ((3, 2, 0, 0, 0, 0, 0), None),
    "X(4,-4,4,0)":  ((4, -4, 4, 0, 0, 0, 0), None),
    "X(4,-3,0,1)":  ((4, -3, 0, 1, 0, 0, 0), (("M2", 2), ("N2", 1), ("K2", -1))),
    "X(4,-3,4,-1)": ((4, -3, 4, -1, 0, 0, 0), None),
    "X(4,-2,0,0)":  ((4, -2, 0, 0, 0, 0, 0), (("M2", 3), ("K2", -1))),
    "X(4,0,-2,0)":  ((4, 0, -2, 0, 0, 0, 0), (("M2", 2), ("S2", 1), ("K2", -1))),
    "X(4,0,-1,0)":  ((4, 0, -1, 0, 0, 0, 0), None),
    "X(4,0,2,0)":   ((4, 0, 2, 0, 0, 0, 0), (("M2", 2), ("S2", -1), ("K2", 1))),
    "X(4,2,-3,0)":  ((4, 2, -3, 0, 0, 0, 0), None),
    "X(4,3,-4,1)":  ((4, 3, -4, 1, 0, 0, 0), (("M2", -1), ("S2", 2), ("N2", 1))),
    "X(6,-3,2,1)":  ((6, -3, 2, 1, 0, 0, 0), (("M2", 3), ("S2", -1), ("N2", 1))),
    "2MNU6":        ((6, -1, 2, -1, 0, 0, 0), None),
    "X(6,0,-2,0)":  ((6, 0, -2, 0, 0, 0, 0), (("M2", 3), ("S2", 1), ("K2", -1))),
    "X(6,1,0,-1)":  ((6, 1, 0, -1, 0, 0, 0), (("M2", 4), ("N2", -1))),
    "X(6,1,0,1)":   ((6, 1, 0, 1, 0, 0, 0), (("M2", 1), ("N2", 1), ("K2", 1))),
    "MSL6":         ((6, 3, -2, -1, 0, 0, 0), (("M2", 3), ("S2", 1), ("N2", -1))),
    "X(8,-1,2,-1)": ((8, -1, 2, -1, 0, 0, 0), None),
    "X(8,3,-2,-1)": ((8, 3, -2, -1, 0, 0, 0), None),
    "X(10,1,-2,1)": ((10, 1, -2, 1, 0, 0, 0), (("M2", 3), ("S2", 1), ("N2", 1))),
}

ONDES.update(SATELLITES_NODAUX)
ONDES.update(RAIES_EMPIRIQUES)

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
# prédiction avec ces constantes : les retirer améliorait l'accord avec
# l'annuaire officiel. Ajustées sur observations (ports REFMAR ci-dessous),
# elles sont en revanche conservées, car elles y améliorent la prédiction.
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


# Constantes REFMAR : analyse harmonique par moindres carrés (outils/ajuster.py)
# des hauteurs d'eau observées, données validées en temps différé au pas de
# 10 min du réseau REFMAR (SHOM, data.shom.fr). H en mètres au-dessus du zéro
# hydrographique, G = déphasage de Greenwich en degrés (UTC).
# Z0 : niveau moyen officiel du SHOM (Références altimétriques maritimes, publié
# dans les métadonnées des marégraphes REFMAR) : c'est celui de l'annuaire
# officiel, dont les hauteurs sont ainsi reproduites. Le niveau moyen observé
# récemment est plus haut (NIVEAU_MOYEN_OBSERVE_2027) ; pour prédire la mer
# réelle plutôt que l'annuaire, prendre cette valeur et, à Cherbourg, les
# constantes REFMAR sans calage (voir CHERBOURG ci-dessous).
#
# Cherbourg : marégraphe REFMAR n° 13, 1995-2025 (1,61 million de mesures),
# 145 ondes ; tendance du niveau moyen +6,03 mm/an.
CONSTANTES_REFMAR_CHERBOURG = {
    "M2": (1.86931, 228.491), "S2": (0.69654, 272.628), "N2": (0.37200, 209.302),
    "K2": (0.19847, 270.539), "M4": (0.14493, 354.623), "K1": (0.09346, 107.120),
    "MS4": (0.08859, 46.897), "MU2": (0.07689, 222.664), "SA": (0.07041, 223.776),
    "L2": (0.06804, 221.924), "NU2": (0.06801, 201.975), "O1": (0.06279, 356.974),
    "MN4": (0.05231, 329.304), "2N2": (0.05042, 189.963), "T2": (0.03670, 261.750),
    "P1": (0.03213, 99.825), "2MS6": (0.02774, 133.809), "LAMBDA2": (0.02697, 209.088),
    "MK4": (0.02597, 46.479), "M6": (0.02507, 95.878), "2SM2": (0.02233, 104.301),
    "SSA": (0.02152, 117.304), "Q1": (0.01738, 316.455), "MA2": (0.01607, 146.813),
    "MSN2": (0.01571, 79.208), "ML4": (0.01531, 7.997), "MNU4": (0.01503, 319.005),
    "EP2": (0.01483, 195.519), "2MN6": (0.01398, 74.155), "X(2,1,0,1)": (0.01274, 4.536),
    "MF": (0.01258, 199.493), "M3": (0.01187, 203.554), "3MS8": (0.01183, 321.114),
    "3MS4": (0.01147, 90.409), "2MSN4": (0.01099, 248.711), "SKM2": (0.01067, 95.418),
    "N4": (0.01066, 306.656), "2MK2": (0.00980, 277.900), "2SM6": (0.00918, 180.827),
    "M8": (0.00902, 277.625), "MSN6": (0.00868, 112.300), "MSK2": (0.00867, 326.882),
    "3M2S2": (0.00838, 14.190), "2MK6": (0.00822, 132.983), "ETA2": (0.00811, 313.648),
    "MM": (0.00803, 221.125), "M2+N": (0.00742, 137.871), "SN4": (0.00707, 39.974),
    "3MN8": (0.00705, 251.935), "M2-N": (0.00697, 143.501), "X(1,0,0,0)": (0.00695, 309.949),
    "3L2": (0.00685, 135.724), "R2": (0.00679, 302.702), "X(2,-3,4,-1)": (0.00674, 215.858),
    "TAU1": (0.00670, 169.935), "X(2,0,-2,2)": (0.00670, 37.993), "SO1": (0.00636, 262.709),
    "2MNS4": (0.00628, 76.344), "X(2,-1,-1,1)": (0.00603, 129.517), "S4": (0.00593, 97.284),
    "MB2": (0.00591, 254.023), "MSM": (0.00577, 19.070), "MSK6": (0.00561, 178.017),
    "J1": (0.00554, 192.387), "MKS2": (0.00549, 310.490), "X(4,2,-3,0)": (0.00530, 103.790),
    "X(2,-3,0,3)": (0.00525, 166.385), "S1": (0.00524, 66.727),
    "X(4,-2,0,0)": (0.00505, 114.331), "X(2,-4,2,0)": (0.00495, 10.918),
    "2MSN8": (0.00487, 296.251), "SGM": (0.00486, 305.754), "OO1": (0.00462, 250.854),
    "X(2,-2,3,0)": (0.00456, 132.784), "2MS8": (0.00440, 4.213), "2NM6": (0.00434, 52.603),
    "X(4,0,-2,0)": (0.00402, 157.463), "SK4": (0.00395, 98.597), "MO3": (0.00390, 224.135),
    "MTM": (0.00361, 216.129), "X(2,-1,1,1)": (0.00361, 136.210),
    "X(4,-4,4,0)": (0.00352, 217.621), "4MS10": (0.00352, 343.510),
    "X(8,3,-2,-1)": (0.00347, 169.484), "X(4,-3,4,-1)": (0.00333, 67.980),
    "MSL6": (0.00331, 322.503), "X(6,1,0,1)": (0.00325, 110.482), "2MNU6": (0.00323, 85.745),
    "X(4,0,2,0)": (0.00318, 159.358), "2Q1": (0.00315, 269.177),
    "X(4,-3,0,1)": (0.00313, 85.015), "X(1,1,-3,0)": (0.00312, 79.125),
    "X(4,3,-4,1)": (0.00308, 252.939), "3MK8": (0.00306, 320.498),
    "X(1,3,-4,0)": (0.00305, 280.736), "X(2,-1,3,-1)": (0.00303, 135.402),
    "2(MN)8": (0.00295, 225.892), "MSF": (0.00295, 195.243), "K2+N": (0.00284, 194.454),
    "X(2,-5,4,1)": (0.00283, 334.143), "X(3,-1,0,1)": (0.00282, 182.035),
    "R3": (0.00279, 268.830), "X(2,6,-6,0)": (0.00279, 324.409),
    "X(2,4,-5,0)": (0.00273, 168.568), "RHO1": (0.00266, 319.021), "T3": (0.00265, 182.226),
    "X(2,-1,0,-1)": (0.00261, 111.306), "X(2,6,-4,0)": (0.00260, 323.113),
    "2MSK8": (0.00254, 1.373), "X(4,0,-1,0)": (0.00250, 307.110),
    "X(2,-3,2,-1)": (0.00247, 226.716), "X(2,1,2,-1)": (0.00246, 252.971),
    "X(2,-6,6,0)": (0.00230, 123.529), "X(8,-1,2,-1)": (0.00226, 254.859),
    "OQ2": (0.00222, 79.417), "X(2,1,-2,-1)": (0.00222, 177.846),
    "X(1,5,-4,0)": (0.00221, 117.261), "X(1,2,0,0)": (0.00212, 300.262),
    "X(1,-1,1,0)": (0.00201, 217.955), "X(2,3,-2,1)": (0.00198, 273.570),
    "X(2,-1,1,0)": (0.00193, 32.617), "X(1,0,-2,1)": (0.00190, 181.243),
    "3M2S10": (0.00189, 24.113), "X(2,-1,0,-2)": (0.00183, 23.840), "M10": (0.00173, 293.435),
    "SO3": (0.00172, 273.364), "4MN10": (0.00163, 265.032), "NO3": (0.00158, 172.611),
    "K2-N": (0.00157, 197.654), "4MS6": (0.00152, 356.851), "N2-N": (0.00147, 133.265),
    "3N2": (0.00146, 91.644), "MK3": (0.00142, 154.481), "2NS2": (0.00135, 170.652),
    "N2+N": (0.00130, 112.325), "2MO5": (0.00128, 260.692), "S3": (0.00123, 4.424),
    "5MS12": (0.00122, 144.409), "4M2S12": (0.00079, 209.405), "2SK2": (0.00079, 30.867),
    "M12": (0.00059, 109.784), "M1": (0.00030, 74.480), "MSQM": (0.00024, 119.053),
    "3MK7": (0.00019, 20.640), "2MK5": (0.00018, 72.692),
}

# Facteurs d'amplitude alignant le marnage de Cherbourg sur l'annuaire officiel :
# avec les seules constantes REFMAR, les PM officielles sortaient 1 à 3,5 cm plus
# haut et les BM 1 à 3,4 cm plus bas (marnage officiel plus grand de 0,9 % en
# mortes-eaux, 1,3 % en vives-eaux). Ajustés sur 54 PM et BM officielles du 15
# au 28/11/2026 ; seul réglage qui tienne en validation croisée entre semaines.
CALAGE_ANNUAIRE_CHERBOURG = {"M2": 1.0135, "S2": 1.0116}


def appliquer_calage(constantes: dict, facteurs: dict) -> dict:
    """Constantes {onde: (H, G)} dont les amplitudes sont multipliées par `facteurs`."""
    return {o: (H * facteurs.get(o, 1.0), G) for o, (H, G) in constantes.items()}


# Port de l'annuaire (niveau moyen officiel, marnage calé). Pour la mer réelle :
# Port("Cherbourg", 49.6514, -1.6355, NIVEAU_MOYEN_OBSERVE_2027["CHERBOURG"],
#      CONSTANTES_REFMAR_CHERBOURG, inferer=False, exclues=frozenset()).
CHERBOURG = Port(
    "Cherbourg", 49.6514, -1.6355, 3.87,
    appliquer_calage(CONSTANTES_REFMAR_CHERBOURG, CALAGE_ANNUAIRE_CHERBOURG),
    inferer=False, exclues=frozenset(),
)

# Brest : marégraphe REFMAR n° 3, 2000-2025 (1,36 million de mesures), 127 ondes ;
# tendance du niveau moyen +6,86 mm/an.
BREST = Port(
    "Brest", 48.3829, -4.4950, 4.14,
    {
        "M2": (2.04967, 108.626), "S2": (0.75249, 148.547), "N2": (0.41728, 90.369),
        "K2": (0.21297, 146.468), "MU2": (0.08785, 104.406), "NU2": (0.07748, 86.131),
        "L2": (0.06545, 99.456), "K1": (0.06490, 75.815), "O1": (0.06427, 327.718),
        "M4": (0.05787, 106.041), "2N2": (0.05759, 72.399), "SA": (0.05107, 243.439),
        "T2": (0.04145, 139.995), "MS4": (0.03557, 179.882), "M6": (0.03120, 355.203),
        "SSA": (0.02975, 99.835), "LAMBDA2": (0.02633, 74.122), "P1": (0.02253, 62.829),
        "MN4": (0.02070, 61.185), "EP2": (0.02065, 87.702), "2MN6": (0.02009, 326.489),
        "Q1": (0.01979, 281.929), "M3": (0.01972, 14.989), "2SM2": (0.01776, 313.707),
        "2MS6": (0.01653, 60.208), "MSN2": (0.01391, 288.764), "X(2,1,0,1)": (0.01326, 228.963),
        "MB2": (0.01145, 104.467), "2MK2": (0.01122, 159.578), "ETA2": (0.01077, 188.266),
        "MK4": (0.01044, 179.788), "MSK2": (0.00991, 187.385), "MA2": (0.00947, 29.609),
        "MF": (0.00878, 192.680), "MM": (0.00878, 230.777), "SKM2": (0.00878, 298.412),
        "3L2": (0.00812, 11.943), "MKS2": (0.00756, 179.660), "X(1,0,0,0)": (0.00708, 265.988),
        "S1": (0.00692, 12.126), "2NM6": (0.00681, 299.144), "2MNU6": (0.00673, 319.250),
        "MSF": (0.00668, 38.661), "X(2,0,-2,2)": (0.00619, 275.003),
        "X(2,-3,0,3)": (0.00618, 48.413), "X(2,-3,4,-1)": (0.00602, 101.504),
        "4MS6": (0.00562, 96.449), "MSM": (0.00550, 13.134), "M2+N": (0.00545, 20.011),
        "MNU4": (0.00541, 64.558), "X(6,-3,2,1)": (0.00540, 74.475),
        "X(2,-2,3,0)": (0.00534, 21.090), "T3": (0.00522, 344.534),
        "X(2,-1,-1,1)": (0.00516, 7.480), "R3": (0.00500, 73.944), "R2": (0.00496, 160.834),
        "X(6,1,0,-1)": (0.00496, 183.970), "ML4": (0.00490, 128.303),
        "X(3,-1,0,1)": (0.00485, 5.839), "2MK6": (0.00452, 62.211), "3M2S2": (0.00443, 285.366),
        "X(2,-1,1,1)": (0.00426, 29.945), "MSN6": (0.00418, 38.586), "M2-N": (0.00413, 3.376),
        "MSL6": (0.00412, 278.148), "2MSN4": (0.00408, 38.179), "3MS8": (0.00404, 324.258),
        "SN4": (0.00403, 149.158), "SGM": (0.00390, 250.391), "X(2,-1,3,-1)": (0.00374, 23.638),
        "2Q1": (0.00372, 237.585), "OQ2": (0.00371, 323.640), "4MS10": (0.00369, 271.202),
        "RHO1": (0.00343, 288.614), "2NS2": (0.00319, 74.226), "N4": (0.00314, 13.849),
        "K2+N": (0.00292, 79.200), "X(3,2,0,0)": (0.00292, 53.276),
        "X(1,-2,0,0)": (0.00288, 239.955), "X(2,-4,2,0)": (0.00277, 274.122),
        "X(6,0,-2,0)": (0.00271, 183.418), "2SM6": (0.00266, 313.986),
        "X(2,1,2,-1)": (0.00265, 135.474), "X(2,4,-5,0)": (0.00250, 22.806),
        "X(6,1,0,1)": (0.00246, 34.835), "X(2,3,-2,1)": (0.00242, 157.647),
        "MTM": (0.00236, 187.200), "X(2,-1,1,0)": (0.00235, 277.379), "J1": (0.00234, 122.400),
        "M10": (0.00227, 225.505), "M8": (0.00226, 239.644), "X(10,1,-2,1)": (0.00225, 245.073),
        "3N2": (0.00223, 295.622), "X(1,1,-3,0)": (0.00223, 25.482), "4MN10": (0.00211, 201.295),
        "X(2,2,1,0)": (0.00204, 62.094), "3MS4": (0.00199, 240.056), "2MK5": (0.00198, 42.690),
        "X(2,-2,4,-2)": (0.00196, 60.121), "MSK6": (0.00181, 294.418), "S4": (0.00180, 270.773),
        "MSQM": (0.00177, 251.492), "MK3": (0.00176, 279.169), "3M2S10": (0.00170, 308.901),
        "2MS8": (0.00154, 54.250), "3MN8": (0.00153, 198.481), "N2-N": (0.00152, 339.476),
        "N2+N": (0.00150, 337.565), "OO1": (0.00148, 206.269), "S3": (0.00141, 155.507),
        "MO3": (0.00139, 57.267), "K2-N": (0.00128, 49.671), "2MSN8": (0.00127, 279.805),
        "2MNS4": (0.00126, 180.950), "3MK8": (0.00106, 321.566), "2MSK8": (0.00096, 53.453),
        "2MO5": (0.00095, 335.446), "SK4": (0.00093, 271.659), "SO3": (0.00092, 154.653),
        "2(MN)8": (0.00073, 160.701), "NO3": (0.00038, 331.872), "2SK2": (0.00034, 237.531),
        "5MS12": (0.00031, 77.809), "3MK7": (0.00021, 110.424), "4M2S12": (0.00017, 115.327),
        "M12": (0.00014, 16.271), "M1": (0.00012, 33.313),
    },
    inferer=False, exclues=frozenset(),
)

# Constantes TICON-4 (Hart-Davis, Dettmering & Seitz, 2025, SEANOE,
# https://doi.org/10.17882/109129, CC BY 4.0), via slackwater-database.
# H en mètres, G = déphasage de Greenwich en degrés (UTC).
#
# Brest, TICON-4 (marégraphe REFMAR n° 3, observations 1846-2025 ; Z0 observé
# sur 05/2006-04/2009) : point de départ des constantes du coefficient
# (BREST_COEFFICIENT) et liste des ondes de base de outils/ajuster.py.
BREST_TICON = Port(
    "Brest (TICON-4)", 48.3829, -4.4950, 4.163,
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

# Niveau moyen observé attendu en 2027 (m), d'après la tendance des moyennes
# annuelles REFMAR des 21 dernières années (outils/ajuster.py).
NIVEAU_MOYEN_OBSERVE_2027 = {"CHERBOURG": 3.938, "BREST": 4.278}

# Ondes semi-diurnes de Brest servant au coefficient de marée. Base TICON-4 ;
# M2, S2, N2, K2, L2 et NU2 sont des valeurs effectives ajustées pour reproduire
# les coefficients officiels (232 coefficients de novembre 2026 à février 2027,
# et la liste ouverte des grandes marées de Saint-Malo 2014-2028, Licence
# Ouverte) : le SHOM calcule le coefficient avec ses propres constantes de
# Brest, dont N2 en particulier diffère de l'analyse TICON-4. La phase de M2,
# indéterminée par cet ajustement, est celle de TICON-4.
BREST_COEFFICIENT = Port(
    "Brest (coefficient)", 48.3829, -4.4950, 0.0,
    {
        "M2": (2.03766, 108.993), "S2": (0.74067, 147.916), "N2": (0.39801, 84.741),
        "K2": (0.21583, 140.968), "MU2": (0.08572, 105.072), "NU2": (0.08291, 94.815),
        "2N2": (0.05703, 72.807), "L2": (0.05077, 133.552), "T2": (0.04169, 138.562),
        "LAMBDA2": (0.02627, 75.818), "EP2": (0.01969, 89.435), "2SM2": (0.01731, 312.051),
        "MKS2": (0.00760, 174.815), "R2": (0.00534, 157.872),
    },
    inferer=False, exclues=frozenset(),
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
    de marée, déduit de la différence de situation de M2 (−2,07 à 10,35 h :
    la marée remonte la Manche depuis Brest ; Brest elle-même est à 0 h)."""
    dg = (port.constantes["M2"][1] - BREST_COEFFICIENT.constantes["M2"][1] + 60.0) % 360.0 - 60.0
    return dg / vitesse("M2")


def pleines_mers_semi_diurnes_brest(debut: datetime, fin: datetime):
    """[(instant UTC, amplitude m)] des PM de la marée semi-diurne de Brest.

    Calculées avec BREST_COEFFICIENT, dont les constantes reproduisent les
    coefficients officiels.
    """
    bruts = _extrema_bruts(BREST_COEFFICIENT, en_utc(debut), en_utc(fin), _semi_diurne_coefficient)
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
    tz = fuseau(tz_nom)
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
    tz = fuseau(tz_nom)
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
    jour = date.fromisoformat(args.date) if args.date else datetime.now(fuseau()).date()
    for i in range(args.jours):
        print()
        afficher_annuaire(jour + timedelta(days=i), PORTS[args.port])


if __name__ == "__main__":
    main()
