#!/bin/bash

# Se déplacer dans le dossier latex généré par Doxygen
cd docs/latex || exit

echo "--- Nettoyage des anciens fichiers ---"
rm -f *.aux *.log *.out *.toc *.idx *.ilg *.ind *.pdf

echo "--- Passe 1 : Génération de la structure ---"
pdflatex -interaction=nonstopmode refman.tex

echo "--- Génération de l'index ---"
makeindex refman.idx

echo "--- Passe 2 : Insertion de l'index et des liens ---"
pdflatex -interaction=nonstopmode refman.tex

echo "--- Passe 3 : Finalisation du PDF ---"
pdflatex -interaction=nonstopmode refman.tex

echo "------------------------------------------------"
echo "Terminé ! Le PDF est sur votre bureau Windows."

