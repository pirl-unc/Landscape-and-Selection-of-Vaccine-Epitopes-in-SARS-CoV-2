#!/bin/bash

set -e

notebooks=(
    "/solvent-accessibility/Normalize-solvent-accessibility-values-for-6vxx.ipynb"
    "/bcell/Extract-SARS2-literature-B-cell-epitopes.ipynb"
    "/bcell/Choose-B-cell-epitopes.ipynb"
    "/bcell/Visualize-B-cell-epitopes.ipynb"
    "/iedb/IEDB-SARS-and-SARSCoV2-MHC-ligands-and-T-cell-epitopes.ipynb"
    "/iedb/Predict-immunogenicity.ipynb"
    "/polymorphism/Extract-spike-sequence-from-whole-proteome-entropy-data.ipynb"
    "/vaccine-peptides/Extract-minimal-epitope-predictions.ipynb"
    "/vaccine-peptides/Select-vaccine-peptides.ipynb"
    "/vaccine-peptides/Summarize-Vaccine-Peptides.ipynb"
    "/vaccine-peptides-all-proteins/Extract-minimal-epitope-predictions.ipynb"
    "/vaccine-peptides-all-proteins/Select-vaccine-peptides.ipynb"
    "/vaccine-peptides-all-proteins/Summarize-Vaccine-Peptides.ipynb"
    "/tcell-validation/Curate-Snyder-MIRA-data.ipynb"
    "/tcell-validation/Map-ORF1ab-coordinates-to-nsp-sub-proteins.ipynb"
    "/tcell-validation/Look-for-overlap-between-vaccine-peptides-and-published-T-cell-epitopes.ipynb"
    "/tcell-validation/Generate-figures-of-T-cell-validation-data.ipynb"
)

for notebook in "${notebooks[@]}"
do
    echo jupyter nbconvert --to notebook --inplace --execute $notebook
    jupyter nbconvert --to notebook --inplace --execute $notebook
done
