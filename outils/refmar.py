#!/usr/bin/env python3
"""Téléchargement des observations marégraphiques REFMAR (data.shom.fr).

Outil de développement : il ne sert qu'à recalculer les constantes
harmoniques et à mesurer la précision. tildes.py n'en dépend pas et reste
100 % hors ligne.

Les observations sont enregistrées dans donnees/refmar/<station>_s<source>_<année>.csv
(répertoire ignoré par git) : une ligne « AAAA-MM-JJTHH:MM:SSZ,hauteur » par
mesure, hauteur en mètres au-dessus du zéro hydrographique, temps UTC.

Sources REFMAR : 1 = brutes haute fréquence (1 min), 2 = brutes temps différé
(10 min), 3 = validées temps différé (10 min), 4 = validées horaires.
Le service limite chaque requête à 31 jours.

    python3 outils/refmar.py 13 3 2000 2025      # Cherbourg, validées, 2000-2025
    python3 outils/refmar.py 13 2 2026 2026      # Cherbourg, brutes, 2026
"""

import argparse
import json
import os
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

URL = ("https://services.data.shom.fr/maregraphie/observation/json/{station}"
       "?sources={source}&dtStart={debut}&dtEnd={fin}")
RACINE = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "donnees", "refmar")
STATIONS = {"CHERBOURG": 13, "BREST": 3}


def _mois(station, source, annee, mois):
    debut = datetime(annee, mois, 1, tzinfo=timezone.utc)
    fin = datetime(annee + (mois == 12), mois % 12 + 1, 1, tzinfo=timezone.utc)
    url = URL.format(station=station, source=source,
                     debut=debut.strftime("%Y-%m-%dT%H:%M:%SZ"),
                     fin=fin.strftime("%Y-%m-%dT%H:%M:%SZ"))
    for essai in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                donnees = json.load(r)["data"]
            break
        except Exception:
            if essai == 4:
                raise
            time.sleep(2 ** essai)
    lignes = []
    for d in donnees:
        t = datetime.strptime(d["timestamp"], "%Y/%m/%d %H:%M:%S")
        if d["value"] is not None and t < fin.replace(tzinfo=None):
            lignes.append(f"{t:%Y-%m-%dT%H:%M:%S}Z,{d['value']}")
    return lignes


def telecharger(station, source, annee):
    os.makedirs(RACINE, exist_ok=True)
    chemin = os.path.join(RACINE, f"{station}_s{source}_{annee}.csv")
    if os.path.exists(chemin):
        return chemin, None
    lignes = []
    for mois in range(1, 13):
        lignes += _mois(station, source, annee, mois)
    with open(chemin + ".part", "w") as f:
        f.write("\n".join(lignes) + ("\n" if lignes else ""))
    os.replace(chemin + ".part", chemin)
    return chemin, len(lignes)


def lire(station, source, annee):
    """[(datetime UTC aware, hauteur m)] d'un fichier téléchargé."""
    chemin = os.path.join(RACINE, f"{station}_s{source}_{annee}.csv")
    serie = []
    with open(chemin) as f:
        for ligne in f:
            t, h = ligne.rstrip().split(",")
            serie.append((datetime.strptime(t, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc), float(h)))
    return serie


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("station", help="numéro REFMAR (13 = Cherbourg, 3 = Brest) ou nom")
    p.add_argument("source", type=int)
    p.add_argument("premiere", type=int)
    p.add_argument("derniere", type=int)
    a = p.parse_args()
    station = STATIONS.get(a.station.upper(), None) or int(a.station)
    annees = range(a.premiere, a.derniere + 1)
    with ThreadPoolExecutor(max_workers=4) as pool:
        for annee, (chemin, n) in zip(annees, pool.map(lambda y: telecharger(station, a.source, y), annees)):
            print(f"{annee} : {'déjà présent' if n is None else f'{n} mesures'} → {os.path.relpath(chemin)}")


if __name__ == "__main__":
    main()
