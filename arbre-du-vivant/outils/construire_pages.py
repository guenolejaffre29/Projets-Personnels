#!/usr/bin/env python3
"""Assemble les pages : atlas.html (atlas commenté) et index.html (catalogue complet).

Usage : python3 construire_pages.py <dossier_sortie> [--fragment]
--fragment produit des pages sans <html>/<head>/<body>, pour une publication en Artifact
(seule la page principale y est enveloppée ; atlas.html reste un document complet).
"""
import pathlib, re, sys

ici = pathlib.Path(__file__).parent
CODES = set("c D R Rb P Pb Cs C Cb Ci Os O Ob Oi Fs F Fb T G S Sb x".split())


def donnees_atlas():
    texte = "\n".join(f.read_text().rstrip("\n") for f in sorted((ici / "atlas").glob("*.txt")))
    assert "</script" not in texte.lower()
    # Vérifie la structure (profondeurs, ancres, codes).
    connus, pile, erreurs = set(), [None], []
    for n, brut in enumerate(texte.split("\n"), 1):
        if not brut.strip() or brut.strip().startswith("#"):
            continue
        if brut.startswith("@@"):
            if brut[2:].strip() not in connus:
                erreurs.append(f"ancre inconnue : {brut}")
            pile = [None]
            continue
        prof = len(brut) - len(brut.lstrip(" "))
        if prof > len(pile) - 1:
            erreurs.append(f"saut de profondeur ligne {n} : {brut.strip()}")
        champs = brut.strip().split("|")
        if champs[0] not in CODES:
            erreurs.append(f"code inconnu ligne {n} : {brut.strip()}")
        connus.add(champs[1].strip())
        pile = pile[: prof + 1] + [champs[1]]
    if erreurs:
        sys.exit("\n".join(erreurs))
    return texte


def document(tete, corps, theme="#ECEFE9"):
    return ('<!doctype html>\n<html lang="fr">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            f'<meta name="theme-color" content="{theme}">\n' + tete.strip() + "\n</head>\n<body>\n" + corps.strip() + "\n</body>\n</html>\n")


def main(sortie, fragment=False):
    sortie = pathlib.Path(sortie)
    sortie.mkdir(parents=True, exist_ok=True)
    atlas = (ici / "atlas_modele.html").read_text().replace("/*__DATA__*/", donnees_atlas())
    tete_a, corps_a = atlas.split("<!--BODY-->")
    (sortie / "atlas.html").write_text(document(tete_a, corps_a))

    style = re.search(r"<style>.*?</style>", tete_a, re.S).group(0)
    style = style.replace("</style>", (ici / "catalogue_extra.css").read_text() + "</style>")
    liens = "\n".join(l for l in tete_a.splitlines() if l.startswith("<link"))
    tete_c = "<title>Catalogue du vivant</title>\n" + liens + "\n" + style
    corps_c = (ici / "catalogue_corps.html").read_text()
    if fragment:
        (sortie / "index.html").write_text(tete_c + "\n" + corps_c)
    else:
        (sortie / "index.html").write_text(document(tete_c, corps_c))


if __name__ == "__main__":
    main(sys.argv[1], "--fragment" in sys.argv)
