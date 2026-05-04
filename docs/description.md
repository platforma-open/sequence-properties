# Overview

Computes physico-chemical properties for peptide and antibody/TCR sequences and emits them as standardized PColumns for Lead Selection ranking. The block detects modality from the input axis automatically — peptide or antibody/TCR — and degrades gracefully with sequencing coverage.

Properties: net charge (pH 7), hydrophobicity (GRAVY), molecular weight, isoelectric point, extinction coefficients (oxidized and reduced), instability index, aliphatic index, aromaticity, and amino acid composition. Peptide mode computes them on the full sequence. Antibody/TCR mode computes them per CDR3 (CDR-H3/L3, or α3/β3 and γ3/δ3 for TCR), per full chain (VH/VL), and at the Fv level for paired antibody chains. Full-chain and Fv columns require all seven VDJ regions (FR1, CDR1, FR2, CDR2, FR3, CDR3, FR4); CDR3-only inputs receive CDR3 properties only.
