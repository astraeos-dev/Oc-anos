#!/usr/bin/env python3
"""Analyse harmonique complète d'un port à partir des observations REFMAR.

Outil de développement (numpy). Étapes :

1. jeu de base : ondes TICON-4 du port, ondes de petits fonds de tildes.ONDES
   (y compris celles que tildes infère quand elles manquent), et satellites
   nodaux libres (± N') de M2, N2 et K2 : à Cherbourg, la modulation nodale de
   M2 est plus faible que la théorie d'équilibre, comme ailleurs en Manche ;
2. ajustement par moindres carrés (niveau, tendance, H et G de chaque onde)
   puis retrait des mesures aberrantes (résidu hors surcote > 6 σ) ;
3. recherche dans le spectre du résidu des raies de marée encore absentes
   (≥ SEUIL_RAIE, bandes diurne à douzième-diurne). Chacune est identifiée par
   ses nombres de Doodson (τ, s, h, p) et reçoit la correction nodale de la
   combinaison d'ondes mères la plus plausible (DECOMPOSITION) ;
4. ajustement final, puis impression des définitions et constantes à reporter
   dans tildes.py.

    python3 outils/ajuster.py 13 1995 2025             # Cherbourg
    python3 outils/ajuster.py 13 1995 2023 --json f    # pour validation
"""

import argparse
import itertools
from datetime import datetime, timezone
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import harmonique as hm  # noqa: E402
import tildes  # noqa: E402

SEUIL_RAIE = 0.002           # m : amplitude minimale d'une raie ajoutée
ECART_MIN = 0.003            # °/h : écart minimal à une onde existante
ONDES_A_SATELLITES = ("M2", "N2", "K2")
MERES = ("M2", "S2", "N2", "K2", "K1", "O1")
# Amplitudes (m) servant à choisir la décomposition d'une raie : valeurs de
# Cherbourg, fixées pour que chaque raie ait la même correction nodale dans
# tous les ports.
AMPLITUDES_REFERENCE = {"M2": 1.87, "S2": 0.696, "N2": 0.371, "K2": 0.198,
                        "K1": 0.0931, "O1": 0.0630, "M4": 0.144}
NOMS_USUELS = {  # raies empiriques ayant un nom consacré
    (2, 4, -2, 0): "SKM2", (2, -4, 4, 0): "3M2S2", (2, -2, 0, 0): "2MK2",
    (1, -1, 2, 0): "TAU1", (1, 3, -2, 0): "SO1", (2, 3, 0, -1): "ETA2",
    (6, 3, -2, -1): "MSL6", (6, -1, 2, -1): "2MNU6",
}


def nom_raie(d):
    return NOMS_USUELS.get(tuple(d), "X({},{},{},{})".format(*d))


def satellites(meres):
    """Définit dans tildes.ONDES les satellites ± N' (sans correction nodale)."""
    noms = []
    for m in meres:
        d = list(tildes.ONDES[m][0])
        for sgn, suf in ((1, "+N"), (-1, "-N")):
            dd = d.copy()
            dd[4] += sgn
            tildes.ONDES[m + suf] = (tuple(dd), None)
            noms.append(m + suf)
    return noms


def decomposer(dood, amplitudes, ordre_max=5):
    """Combinaison entière de MERES reproduisant (c, s, h, p), la plus plausible.

    Une onde composée a une amplitude de l'ordre du produit de celles de ses
    ondes mères (rapportées à M2), divisé par un facteur ~ M2/M4 par ordre
    d'interaction : on maximise ce score. Retourne (règle nodale, constante de
    phase) ou (None, 0) si aucune combinaison ne convient.
    """
    D = {m: tildes.ONDES[m][0] for m in MERES}
    lg = {m: math.log(amplitudes[m] / amplitudes["M2"]) for m in MERES}
    penalite = math.log(amplitudes["M2"] / amplitudes["M4"])
    meilleur = None
    for k in itertools.product(range(-4, 5), repeat=len(MERES)):
        ordre = sum(abs(x) for x in k)
        if ordre == 0 or ordre > ordre_max:
            continue
        v = [sum(kk * D[m][i] for kk, m in zip(k, MERES)) for i in range(4)]
        if tuple(v) != tuple(dood):
            continue
        score = sum(abs(kk) * lg[m] for kk, m in zip(k, MERES)) - penalite * (ordre - 1)
        if meilleur is None or score > meilleur[0]:
            phase = sum(kk * D[m][6] for kk, m in zip(k, MERES)) % 4
            meilleur = (score, tuple((m, kk) for m, kk in zip(MERES, k) if kk), phase)
    return (meilleur[1], meilleur[2]) if meilleur else (None, 0)


def raies_du_residu(t, r, existantes):
    """[(c, s, h, p)] des raies du spectre du résidu ≥ SEUIL_RAIE, non modélisées."""
    g = np.round((t - t[0]) * 144).astype(int)  # grille de 10 min, trous à 0
    x = np.zeros(g[-1] + 1)
    x[g] = r - r.mean()
    amp = np.abs(np.fft.rfft(x * np.hanning(len(x)))) * 4.0 / len(x)
    freq = np.fft.rfftfreq(len(x), d=1.0 / 6.0) * 360.0
    W = tildes.VITESSES_ARGUMENTS
    wex = list(hm.vitesses(existantes))
    trouvees = []
    for i in np.argsort(-amp):
        if amp[i] < SEUIL_RAIE:
            break
        f = freq[i]
        if not 10.0 < f < 180.0 or amp[i] != amp[max(0, i - 2):i + 3].max():
            continue
        c = int(round(f / 15.0))
        cand = [(abs(s) + abs(h) + abs(p), (c, s, h, p))
                for s in range(-6, 7) for h in range(-6, 7) for p in range(-3, 4)
                if abs(c * W[0] + s * W[1] + h * W[2] + p * W[3] - f) < 0.0011]
        if not cand:
            continue
        d = min(cand)[1]
        w = sum(k * v for k, v in zip(d, W))
        if min(abs(np.array(wex) - w)) < ECART_MIN:
            continue
        trouvees.append(d)
        wex.append(w)
    return trouvees


def ajuster(station, annees, base, jour_ref, journal=print):
    t, h = hm.charger(station, 3, annees)
    journal(f"{len(t)} mesures {annees[0]}-{annees[-1]}")
    ondes = list(base) + satellites(ONDES_A_SATELLITES)
    z0, pente, cst, rms = hm.analyser(t, h, ondes, jour_ref=jour_ref)
    # mesures aberrantes : résidu hors surcote (moyenne glissante 25 h retirée) > 6 σ
    r = h - hm.predire(t, z0, cst) - pente / 365.25 * (t - jour_ref)
    hf = r - np.convolve(r, np.ones(151) / 151, mode="same")
    garde = np.abs(hf) < 6 * hf.std()
    t, h = t[garde], h[garde]
    journal(f"{np.sum(~garde)} mesures aberrantes retirées ; rms {rms * 100:.2f} cm")
    z0, pente, cst, rms = hm.analyser(t, h, ondes, jour_ref=jour_ref)
    r = h - hm.predire(t, z0, cst) - pente / 365.25 * (t - jour_ref)
    raies = []
    for d in raies_du_residu(t, r, ondes):
        regle, phase = decomposer(d, AMPLITUDES_REFERENCE)
        nom = nom_raie(d)
        tildes.ONDES[nom] = (tuple(d) + (0, 0, phase), regle)
        raies.append(nom)
    journal(f"{len(raies)} raies ajoutées : {' '.join(raies)}")
    ondes += raies
    z0, pente, cst, rms = hm.analyser(t, h, ondes, jour_ref=jour_ref)
    journal(f"{len(ondes)} ondes ; rms {rms * 100:.2f} cm ; Z0 {z0:.4f} m ; tendance {pente * 1000:+.2f} mm/an")
    return z0, pente, cst, raies


ANNEE_Z0 = 2027
ANS_TENDANCE = 21


def niveau_moyen(station, annees, cst):
    """Niveau moyen en ANNEE_Z0 : régression des moyennes annuelles observées
    (marée retirée) sur les ANS_TENDANCE dernières années. La hausse du niveau
    marin s'est accélérée : une tendance sur 30 ans sous-estime le niveau actuel.
    """
    x, y = [], []
    for a in list(annees)[-ANS_TENDANCE:]:
        t, h = hm.charger(station, 3, [a])
        if len(t) > 40000:  # année presque complète
            x.append(a + 0.5)
            y.append(float(np.mean(h - hm.predire(t, 0.0, cst))))
    pente, origine = np.polyfit(x, y, 1)
    return pente * ANNEE_Z0 + origine, pente


def code_python(cst, raies):
    """Texte des définitions d'ondes et des constantes, à reporter dans tildes.py."""
    lignes = ["# --- définitions (ONDES) ---"]
    sat = [o for o in cst if o.endswith(("+N", "-N"))]
    for o in sorted(raies + sat, key=tildes.vitesse):
        d, regle = tildes.ONDES[o]
        lignes.append(f'    "{o}": {" " * max(0, 8 - len(o))}({d}, {regle!r}),')
    lignes.append("# --- constantes {onde: (H m, G °)} ---")
    items = sorted(cst.items(), key=lambda kv: -kv[1][0])
    for i in range(0, len(items), 3):
        lignes.append("        " + " ".join(f'"{o}": ({H:.5f}, {G:.3f}),' for o, (H, G) in items[i:i + 3]))
    return "\n".join(lignes)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("station", type=int)
    p.add_argument("premiere", type=int)
    p.add_argument("derniere", type=int)
    p.add_argument("--json", help="enregistre le résultat (constantes et définitions)")
    a = p.parse_args()
    # Jeu de base : les 50 ondes de TICON-4 (mêmes pour tous les ports), et non
    # les constantes actuelles de tildes.py, pour que l'analyse soit reproductible.
    base = sorted(set(tildes.BREST_TICON.observees) | set(tildes.INFERENCES)
                  | set(PETITS_FONDS_SUPPLEMENTAIRES), key=tildes.vitesse)
    jour_ref = hm.jours_depuis_j2000([datetime(ANNEE_Z0, 1, 1, tzinfo=timezone.utc)])[0]
    annees = range(a.premiere, a.derniere + 1)
    z0, pente, cst, raies = ajuster(a.station, annees, base, jour_ref)
    z0_recent, pente_recente = niveau_moyen(a.station, annees, cst)
    print(f"Z0 en {ANNEE_Z0} : {z0_recent:.4f} m (tendance des {ANS_TENDANCE} dernières années "
          f"{pente_recente * 1000:+.2f} mm/an) ; tendance sur toute la période : {z0:.4f} m")
    if a.json:
        json.dump({"z0": z0_recent, "pente": pente_recente, "constantes": cst, "raies": raies,
                   "definitions": {o: tildes.ONDES[o] for o in cst}}, open(a.json, "w"))
    print(code_python(cst, raies))


# Ondes de petits fonds et astronomiques absentes de TICON-4 mais définies
# dans tildes.ONDES, ajustées en plus des ondes TICON et INFERENCES.
PETITS_FONDS_SUPPLEMENTAIRES = (
    "MSM", "OQ2", "2NS2", "MSK2", "2SK2", "MSN2", "MO3", "NO3", "SO3", "MK3", "2MNS4", "3MS4",
    "MNU4", "ML4", "2MSN4", "2NM6", "4MS6", "3MK7", "2(MN)8", "3MK8", "2MSK8", "4MN10", "M10",
    "4MS10", "3M2S10", "M12", "5MS12", "4M2S12",
)

if __name__ == "__main__":
    main()
