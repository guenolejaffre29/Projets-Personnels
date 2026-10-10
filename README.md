# Projets personnels

- [Arbre du vivant](arbre-du-vivant/) : explorateur mobile du vivant, sans les plantes.
  - `index.html` : **catalogue complet**, ~2,1 millions d'espèces décrites d'après la version étendue du Catalogue of Life, chargées à la demande depuis `d/`.
  - `atlas.html` : atlas commenté, arbre phylogénétique d'environ 5 500 groupes avec noms français et anecdotes.
  - `outils/` : sources de l'atlas et scripts de construction. Le workflow `.github/workflows/catalogue.yml` télécharge le Catalogue of Life et régénère `d/`. Les pages se servent depuis un serveur web (pas en `file://`).
