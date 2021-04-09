#!/bin/bash

set -ex

notebooks=(
    "bcell/Extract-SARS2-literature-B-cell-epitopes.ipynb"
    "bcell/Choose-B-cell-epitopes.ipynb"
    "bcell/Visualize-B-cell-epitopes.ipynb"
    "iedb/IEDB-SARS-and-SARSCoV2-MHC-ligands-and-T-cell-epitopes.ipynb"
    "iedb/Predict-immunogenicity.ipynb"
    "polymorphism/Extract-spike-sequence-from-whole-proteome-entropy-data.ipynb"
    "vaccine-peptides/Extract-minimal-epitope-predictions.ipynb"
    "vaccine-peptides/Select-vaccine-peptides.ipynb"
    "vaccine-peptides/Summarize-Vaccine-Peptides.ipynb"
    "vaccine-peptides-all-proteins/Extract-minimal-epitope-predictions.ipynb"
    "vaccine-peptides-all-proteins/Select-vaccine-peptides.ipynb"
    "vaccine-peptides-all-proteins/Summarize-Vaccine-Peptides.ipynb"
    "tcell-validation/Curate-Snyder-MIRA-data.ipynb"
    "tcell-validation/Map-ORF1ab-coordinates-to-nsp-sub-proteins.ipynb"
    "tcell-validation/Look-for-overlap-between-vaccine-peptides-and-published-T-cell-epitopes.ipynb"
    "tcell-validation/Generate-figures-of-T-cell-validation-data.ipynb"
)

for notebook in "${notebooks[@]}"
do
    if [ ! -f "$notebook" ]; then
        echo "$notebook does not exist!"
        exit 0;
    else
        echo "OK exists: $notebook"
    fi
done

for notebook in "${notebooks[@]}"
do
    jupyter nbconvert --ExecutePreprocessor.timeout=1200 --to notebook --inplace --execute $notebook

done
