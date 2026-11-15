#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculateur harmonique de marées pour Cherbourg-en-Cotentin (100% hors ligne).
Calcule avec précision les Pleines Mers, Basses Mers, hauteurs et coefficients
par le calcul astronomique pur (Méthode de Doodson / Schureman).
"""

import math
from datetime import datetime, date, time, timedelta, timezone
from zoneinfo import ZoneInfo

# ==============================================================================
# CONSTANTES HARMONIQUES OFFICIELLES DE CHERBOURG
# ==============================================================================
Z0_CHERBOURG = 3.86              # Niveau moyen au-dessus du Zéro Hydrographique (m)
UNITE_HAUTEUR_CHERBOURG = 2.62   # Unité de hauteur de vive-eau pour les coefficients

CONSTITUANTS_CHERBOURG = [
    # Ondes semi-diurnes
    ("M2",  1.860, 228.5, 2, -2,  2,  0,   0.0, "M2"),
    ("S2",  0.680, 268.0, 2,  0,  0,  0,   0.0, "S2"),
    ("N2",  0.360, 210.0, 2, -3,  2,  1,   0.0, "M2"),
    ("K2",  0.190, 266.0, 2,  0,  2,  0,   0.0, "K2"),
    ("2N2", 0.050, 191.0, 2, -4,  2,  2,   0.0, "M2"),
    ("MU2", 0.070, 260.0, 2, -4,  4,  0,   0.0, "M2"),
    ("NU2", 0.070, 213.0, 2, -3,  4, -1,   0.0, "M2"),
    ("L2",  0.060, 245.0, 2, -1,  2, -1, 180.0, "M2"),
    
    # Ondes diurnes
    ("K1",  0.090, 108.0, 1,  0,  1,  0,  90.0, "K1"),
    ("O1",  0.070, 355.0, 1, -2,  1,  0, -90.0, "O1"),
    ("P1",  0.030, 105.0, 1,  0, -1,  0, -90.0, "S2"),
    ("Q1",  0.020, 310.0, 1, -3,  1,  1, -90.0, "O1"),
    
    # Ondes de petits fonds (Manche / Asymétrie flot-jusant)
    ("M4",  0.140, 345.0, 4, -4,  4,  0,   0.0, "M4"),
    ("MS4", 0.090,  38.0, 4, -2,  2,  0,   0.0, "MS4"),
    ("MN4", 0.050, 325.0, 4, -5,  4,  1,   0.0, "M4"),
    ("M6",  0.040, 115.0, 6, -6,  6,  0,   0.0, "M6"),
]


# ==============================================================================
# MOTEUR ASTRONOMIQUE ET NODAL
# ==============================================================================
def julian_centuries_j2000(dt_utc: datetime) -> float:
    epoch_j2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=timezone.utc)
    return (dt_utc - epoch_j2000).total_seconds() / 86400.0 / 36525.0


def parametres_astronomiques(dt_utc: datetime) -> dict:
    T = julian_centuries_j2000(dt_utc)
    s = (218.3164477 + 481267.88123421 * T - 0.0015786 * T**2) % 360.0
    h = (280.4664600 + 36000.76982779 * T + 0.0003032 * T**2) % 360.0
    p = (83.3532465  + 4069.01372870 * T - 0.0103200 * T**2) % 360.0
    N = (125.0445200 - 1934.13626100 * T + 0.0020708 * T**2) % 360.0
    
    heures_utc = dt_utc.hour + dt_utc.minute / 60.0 + dt_utc.second / 3600.0
    theta_0 = (15.0 * heures_utc) % 360.0
    
    return {"s": s, "h": h, "p": p, "N": N, "theta_0": theta_0}


def corrections_nodales(N_deg: float) -> dict:
    N_rad = math.radians(N_deg)
    cos_N, sin_N = math.cos(N_rad), math.sin(N_rad)
    f_M2 = 1.000 - 0.037 * cos_N
    u_M2 = -2.1 * sin_N
    return {
        "S2":  (1.0, 0.0),
        "M2":  (f_M2, u_M2),
        "K2":  (1.024 + 0.286 * cos_N, -17.7 * sin_N),
        "K1":  (1.006 + 0.115 * cos_N, -8.9 * sin_N),
        "O1":  (1.009 + 0.187 * cos_N, +10.8 * sin_N),
        "M4":  (f_M2**2, 2.0 * u_M2),
        "MS4": (f_M2, u_M2),
        "M6":  (f_M2**3, 3.0 * u_M2),
    }


def hauteur_maree(dt: datetime) -> float:
    dt_utc = dt.astimezone(timezone.utc) if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
    astro = parametres_astronomiques(dt_utc)
    nodal = corrections_nodales(astro["N"])
    
    hauteur = Z0_CHERBOURG
    for _, H, G, c_theta, c_s, c_h, c_p, offset, fam in CONSTITUANTS_CHERBOURG:
        f, u = nodal[fam]
        V = c_theta * astro["theta_0"] + c_s * astro["s"] + c_h * astro["h"] + c_p * astro["p"] + offset
        hauteur += f * H * math.cos(math.radians((V + u - G) % 360.0))
    return hauteur


def hauteur_semi_diurne_pure(dt: datetime) -> float:
    dt_utc = dt.astimezone(timezone.utc) if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
    astro = parametres_astronomiques(dt_utc)
    nodal = corrections_nodales(astro["N"])
    h_sd = 0.0
    for _, H, G, c_theta, c_s, c_h, c_p, offset, fam in CONSTITUANTS_CHERBOURG:
        if c_theta == 2:
            f, u = nodal[fam]
            V = c_theta * astro["theta_0"] + c_s * astro["s"] + c_h * astro["h"] + c_p * astro["p"] + offset
            h_sd += f * H * math.cos(math.radians((V + u - G) % 360.0))
    return h_sd


# ==============================================================================
# RECHERCHE DES EXTREMA (PM & BM)
# ==============================================================================
def trouver_extrema_journee(jour: date, tz_nom: str = "Europe/Paris") -> list:
    tz = ZoneInfo(tz_nom)
    debut = datetime.combine(jour, time(0, 0), tzinfo=tz) - timedelta(minutes=45)
    pas_minutes = 5
    nb_points = int((26 * 60) / pas_minutes) + 1
    
    instants = [debut + timedelta(minutes=i * pas_minutes) for i in range(nb_points)]
    hauteurs = [hauteur_maree(t) for t in instants]
    
    extrema = []
    for i in range(1, len(hauteurs) - 1):
        y_prev, y_curr, y_next = hauteurs[i-1], hauteurs[i], hauteurs[i+1]
        est_pm = (y_curr >= y_prev) and (y_curr >= y_next) and (y_curr != y_prev or y_curr != y_next)
        est_bm = (y_curr <= y_prev) and (y_curr <= y_next) and (y_curr != y_prev or y_curr != y_next)
        
        if est_pm or est_bm:
            denom = 2.0 * (y_prev - 2.0 * y_curr + y_next)
            delta = (y_prev - y_next) / denom if abs(denom) > 1e-12 else 0.0
            t_exact = instants[i] + timedelta(minutes=delta * pas_minutes)
            
            if t_exact.date() == jour:
                h_exact = hauteur_maree(t_exact)
                type_m = "PM" if est_pm else "BM"
                coeff = None
                if est_pm:
                    amp_sd = abs(hauteur_semi_diurne_pure(t_exact))
                    coeff = int(round(100.0 * amp_sd / UNITE_HAUTEUR_CHERBOURG))
                    coeff = max(20, min(120, coeff))
                
                extrema.append({
                    "type": type_m,
                    "heure": t_exact,
                    "hauteur": h_exact,
                    "coeff": coeff
                })
    return extrema


def afficher_annuaire(jour: date):
    extrema = trouver_extrema_journee(jour)
    print("=" * 65)
    print(f" ANNUAIRE DES MARÉES (OFFLINE) — CHERBOURG : {jour.strftime('%d/%m/%Y')}")
    print("=" * 65)
    print(f"{'Marée':<14} | {'Heure (CET)':<12} | {'Hauteur':<12} | {'Coefficient'}")
    print("-" * 65)
    for e in extrema:
        libelle = "Pleine Mer" if e["type"] == "PM" else "Basse Mer"
        h_str = e["heure"].strftime("%Hh%M")
        hauteur_str = f"{e['hauteur']:.2f} m"
        c_str = str(e["coeff"]) if e["coeff"] is not None else "-"
        print(f"{libelle:<14} | {h_str:<12} | {hauteur_str:<12} | {c_str:^11}")
    print("=" * 65)


if __name__ == "__main__":
    auj = datetime.now(ZoneInfo("Europe/Paris")).date()
    demain = auj + timedelta(days=1)
    
    print("\n")
    afficher_annuaire(auj)
    print("\n")
    afficher_annuaire(demain)
