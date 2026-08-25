# Sequence Properties

Compute the physico-chemical properties of your sequences. This Platforma block calculates net charge, isoelectric point, hydrophobicity, molecular weight, extinction coefficients, instability and aliphatic indices, aromaticity, and amino acid composition for peptides and antibody or TCR sequences — emitted as columns you can rank, filter, and plot on.

Open-source analysis block for Platforma, the biologics discovery platform by MiLaboratories. For the full no-code workflow, see [platforma.bio](https://platforma.bio/).

## What it does

Whether a candidate can actually be made and used depends partly on properties you can read straight off the sequence. A peptide's charge at physiological pH affects solubility and formulation; an antibody's isoelectric point affects purification behavior and viscosity; hydrophobicity relates to aggregation propensity. These are cheap to compute and useful to have on every candidate before spending anything on them.

The block computes the standard set: **net charge** at pH 7, **hydrophobicity** as GRAVY, **molecular weight**, **isoelectric point**, **extinction coefficients** both oxidized and reduced, **instability index**, **aliphatic index**, **aromaticity**, and **amino acid composition**.

Modality is detected from the input, and the scope adapts to it. Peptides are computed on the full sequence. Antibody and TCR input is computed at three levels: per CDR3 — CDR-H3 and CDR-L3 for antibodies, α3/β3 and γ3/δ3 for TCRs — per full chain (VH, VL), and at the Fv level for paired antibody chains.

Coverage is handled gracefully rather than by failing. Full-chain and Fv properties need all seven V(D)J regions present (FR1 through FR4 and CDR1 through CDR3); a CDR3-only input simply gets CDR3 properties, without an error. You get what the data supports.

Distribution and relationship views let you see how properties spread across a library and how they covary before using them as filters.

## Inputs & outputs

* **Input:** peptide sequences from [Peptide Profiling](https://github.com/platforma-open/peptide-extraction), or antibody and TCR clonotypes from any Platforma clonotyping or import block. Modality is detected automatically.
* **Output:** standardized property columns per sequence — at CDR3, full-chain, and Fv level where the input supports it — available to any downstream block for ranking, filtering, plotting, and modeling.

## Specifications

| | |
|---|---|
| Block title in app | Sequence Properties |
| Properties | Net charge (pH 7), hydrophobicity (GRAVY), molecular weight, isoelectric point, extinction coefficient (oxidized and reduced), instability index, aliphatic index, aromaticity, amino acid composition |
| Peptide scope | Full sequence |
| Antibody / TCR scope | Per CDR3 (CDR-H3/L3; α3/β3, γ3/δ3), per full chain (VH, VL), and Fv level for paired chains |
| Coverage requirement | Full-chain and Fv properties need all seven V(D)J regions; CDR3-only inputs get CDR3 properties |
| Modality detection | Automatic from the input |
| Implementation | [Biopython](https://biopython.org/) |
| Views | Main table, property distributions, property relationships |

## Use cases

* **Formulation triage:** flag candidates whose charge or isoelectric point will complicate purification or formulation.
* **Aggregation risk:** use hydrophobicity to identify candidates more likely to aggregate.
* **Peptide developability:** compute charge, hydrophobicity, and molecular weight across a peptide library before synthesis.
* **Lead ranking:** include property columns as criteria in [Lead Selection](https://github.com/platforma-open/antibody-tcr-lead-selection).
* **Quantification setup:** use extinction coefficients to plan concentration measurement for expressed candidates.
* **Library characterization:** read property distributions to see whether a library is biased in charge or hydrophobicity.
* **Feature input:** supply numeric per-sequence features to downstream modeling or plotting.

## FAQ

### Which properties are computed?

Net charge at pH 7, GRAVY hydrophobicity, molecular weight, isoelectric point, extinction coefficients (oxidized and reduced), instability index, aliphatic index, aromaticity, and amino acid composition.

### At what level are they computed for antibodies?

Three levels: per CDR3, per full chain, and at the Fv level for paired chains. That means you can compare candidates on the binding loop alone, on a whole domain, or on the assembled Fv, depending on what the property means for your question.

### What if my clonotypes only have CDR3?

You get CDR3 properties. Full-chain and Fv properties require all seven V(D)J regions to be present; when they are not, the block reports what it can rather than failing.

### Does it work on peptides?

Yes. Modality is detected from the input, and peptides are computed on the full sequence — no region decomposition, since there are no regions.

### Why does isoelectric point matter?

It governs how a molecule behaves during purification and how it charges at a given pH, which relates to solubility and, for antibodies, viscosity at high concentration. It is one of the cheapest early signals of downstream manufacturing difficulty.

### What is GRAVY?

The grand average of hydropathy — mean hydrophobicity across the sequence. Higher values indicate a more hydrophobic sequence, which correlates with aggregation propensity and can affect expression and solubility.

### How does this relate to the liabilities blocks?

Properties are continuous physico-chemical measures of the whole sequence. [Sequence Liabilities](https://github.com/platforma-open/antibody-sequence-liabilities) instead flags specific problematic motifs at specific positions. Together they cover bulk behavior and localized chemical risk.

## Citation

Properties are computed with Biopython. If you use this block in your research, please cite:

> Cock, P. J. A., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., & de Hoon, M. J. L. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* **25**(11), 1422–1423. [https://doi.org/10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)

## Part of the Platforma ecosystem

This block is part of [Platforma](https://platforma.bio/) by [MiLaboratories](https://github.com/milaboratory), built on [Biopython](https://biopython.org/). Explore the other open-source blocks at [github.com/platforma-open](https://github.com/platforma-open) and the docs for antibody discovery at [docs.platforma.bio/biology-guides/antibody-discovery](https://docs.platforma.bio/biology-guides/antibody-discovery/).
