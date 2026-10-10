# Projets personnels

- [Arbre du vivant](arbre-du-vivant/) : explorateur mobile du vivant, sans les plantes.
  - `index.html` : **catalogue complet**, ~1,88 million d'espèces décrites d'après le Catalogue of Life, chargées à la demande depuis `d/`.
  - `atlas.html` : atlas commenté, arbre phylogénétique d'environ 5 500 groupes avec noms français et anecdotes.
  - `outils/` : sources de l'atlas et scripts de construction. Le workflow `.github/workflows/catalogue.yml` télécharge le Catalogue of Life et régénère `d/`. Les pages se servent depuis un serveur web (pas en `file://`).
