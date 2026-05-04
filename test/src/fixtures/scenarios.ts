/**
 * Synthetic-corpus scenario builders for sequence-properties block tests.
 *
 * Each scenario has:
 *   - a TSV file in `../assets/<name>.tsv` (committed, diffable)
 *   - an xsv-import `BlockArgs.spec` builder that publishes axes + columns
 *     with the exact name / domain / annotations seqprops's bb.addMulti
 *     queries expect to discover.
 *
 * This is the "synthetic-PColumn-publisher" route: xsv-import exposes full
 * control over axis + column specs, so we can mint VDJ-shaped or peptide-
 * shaped anchored data without running MiXCR or peptide-extraction. Tests
 * stay deterministic and fast (xsv-import is a thin TSV→PFrame import).
 *
 * The MiXCR canary test (`mixcrClonotyping` upstream) lives outside this
 * synthetic path — it runs the real fastq pipeline, guarding the contract
 * that the synthetic axis specs here actually match what MiXCR emits.
 */

import type { Spec } from '@platforma-open/milaboratories.xsv-import.model';

type Feature = 'FR1' | 'CDR1' | 'FR2' | 'CDR2' | 'FR3' | 'CDR3' | 'FR4';

/** Anchor abundance column annotations matching seqprops's inputAnchorSpecs. */
const ANCHOR_ANNOTATIONS = {
  'pl7.app/isAnchor': 'true',
  'pl7.app/isAbundance': 'true',
  'pl7.app/abundance/isPrimary': 'true',
  'pl7.app/abundance/normalized': 'false',
  'pl7.app/abundance/unit': 'cells',
  'pl7.app/label': 'Number of Cells',
};

/** Build a per-region VDJ sequence column spec (matches MiXCR's emission shape). */
function vdjSequenceColumn(args: {
  column: string;
  chain: 'A' | 'B';
  feature: Feature | 'FR4InFrame';
  receptor: 'IG' | 'TCRAB' | 'TCRGD';
}): Spec['columns'][number] {
  return {
    column: args.column,
    spec: {
      name: 'pl7.app/vdj/sequence',
      valueType: 'String',
      domain: {
        'pl7.app/alphabet': 'aminoacid',
        'pl7.app/vdj/feature': args.feature,
        'pl7.app/vdj/scClonotypeChain': args.chain,
        'pl7.app/vdj/scClonotypeChain/index': 'primary',
        'pl7.app/vdj/receptor': args.receptor,
      },
      annotations: {
        'pl7.app/label': `${args.chain} ${args.feature} aa Primary`,
      },
    },
  };
}

/** Peptide mode: peptide.tsv. Drives detection via variantKey + extractionRunId. */
export function peptideScenario() {
  return {
    tsv: './assets/peptide.tsv',
    fileExt: 'tsv' as const,
    spec: {
      axes: [
        {
          column: 'sampleId',
          spec: { name: 'pl7.app/sampleId', type: 'String' as const },
        },
        {
          column: 'variantKey',
          spec: {
            name: 'pl7.app/variantKey',
            type: 'String' as const,
            domain: { 'pl7.app/peptide/extractionRunId': 'synthetic-pep-run' },
            annotations: { 'pl7.app/label': 'Variant ID' },
          },
        },
      ],
      columns: [
        {
          column: 'cellCount',
          spec: {
            name: 'pl7.app/abundance',
            valueType: 'Long' as const,
            annotations: ANCHOR_ANNOTATIONS,
          },
        },
        {
          column: 'sequence',
          spec: {
            name: 'pl7.app/sequence',
            valueType: 'String' as const,
            domain: {
              'pl7.app/feature': 'peptide',
              'pl7.app/alphabet': 'aminoacid',
            },
            annotations: { 'pl7.app/label': 'Peptide sequence' },
          },
        },
      ],
      storageFormat: 'Binary' as const,
      partitionKeyLength: 0,
    } satisfies Spec,
  };
}

/** Build the column spec list for a VDJ scenario (paired chain) with a chosen receptor. */
function vdjColumnsFull(receptor: 'IG' | 'TCRAB' | 'TCRGD'): Spec['columns'] {
  // FR4 is emitted by MiXCR as 'FR4InFrame'. seqprops normalises both back
  // to 'FR4' in REQUIRED_FEATURES. We use 'FR4InFrame' here to mirror the
  // real producer and keep this deviation covered by the synthetic path.
  const featuresOnDisk = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4InFrame'] as const;
  const cols: Spec['columns'] = [
    {
      column: 'cellCount',
      spec: {
        name: 'pl7.app/vdj/uniqueCellCount',
        valueType: 'Long' as const,
        annotations: ANCHOR_ANNOTATIONS,
      },
    },
  ];
  for (const chain of ['A', 'B'] as const) {
    for (const f of featuresOnDisk) {
      cols.push(
        vdjSequenceColumn({
          column: `${chain}_${f.replace('InFrame', '')}`,
          chain,
          feature: f,
          receptor,
        }),
      );
    }
  }
  return cols;
}

/** Common axes for sc clonotype-keyed scenarios. */
function scAxes(receptor: 'IG' | 'TCRAB' | 'TCRGD'): Spec['axes'] {
  return [
    {
      column: 'sampleId',
      spec: { name: 'pl7.app/sampleId', type: 'String' as const },
    },
    {
      column: 'scClonotypeKey',
      spec: {
        name: 'pl7.app/vdj/scClonotypeKey',
        type: 'String' as const,
        domain: {
          'pl7.app/vdj/clonotypingRunId': 'synthetic-vdj-run',
          'pl7.app/vdj/receptor': receptor,
        },
        annotations: { 'pl7.app/label': 'Clonotype ID' },
      },
    },
  ];
}

/** sc IG, full 7-region paired coverage. Drives full_chain + Fv path. */
export function scIgFullScenario() {
  return {
    tsv: './assets/sc-ig-full.tsv',
    fileExt: 'tsv' as const,
    spec: {
      axes: scAxes('IG'),
      columns: vdjColumnsFull('IG'),
      storageFormat: 'Binary' as const,
      partitionKeyLength: 0,
    } satisfies Spec,
  };
}

/** sc IG, CDR3-only. Drives the R11a info-banner path. */
export function scIgCdr3OnlyScenario() {
  return {
    tsv: './assets/sc-ig-cdr3-only.tsv',
    fileExt: 'tsv' as const,
    spec: {
      axes: scAxes('IG'),
      columns: [
        {
          column: 'cellCount',
          spec: {
            name: 'pl7.app/vdj/uniqueCellCount',
            valueType: 'Long' as const,
            annotations: ANCHOR_ANNOTATIONS,
          },
        },
        vdjSequenceColumn({ column: 'A_CDR3', chain: 'A', feature: 'CDR3', receptor: 'IG' }),
        vdjSequenceColumn({ column: 'B_CDR3', chain: 'B', feature: 'CDR3', receptor: 'IG' }),
      ],
      storageFormat: 'Binary' as const,
      partitionKeyLength: 0,
    } satisfies Spec,
  };
}

/** sc TCRαβ, full coverage. Asserts no Fv emission for TCR receptor. */
export function scTcrAbFullScenario() {
  return {
    tsv: './assets/sc-tcrab-full.tsv',
    fileExt: 'tsv' as const,
    spec: {
      axes: scAxes('TCRAB'),
      columns: vdjColumnsFull('TCRAB'),
      storageFormat: 'Binary' as const,
      partitionKeyLength: 0,
    } satisfies Spec,
  };
}

/**
 * Per-feature VDJ sequence column matching bulk MiXCR's emission: no
 * `scClonotypeChain`, no `pl7.app/vdj/receptor`. Chain identity lives on
 * the axis-domain `pl7.app/vdj/chain` key — SD-008 derives the receptor
 * from there.
 */
function bulkVdjSequenceColumn(
  column: string,
  feature: Feature | 'FR4InFrame',
): Spec['columns'][number] {
  return {
    column,
    spec: {
      name: 'pl7.app/vdj/sequence',
      valueType: 'String',
      domain: {
        'pl7.app/alphabet': 'aminoacid',
        'pl7.app/vdj/feature': feature,
      },
      annotations: {
        'pl7.app/label': `${feature} aa`,
      },
    },
  };
}

/** Common axes for bulk clonotypeKey-keyed scenarios. */
function bulkAxes(chain: 'IGHeavy' | 'IGLight' | 'TCRAlpha' | 'TCRBeta' | 'TCRGamma' | 'TCRDelta'): Spec['axes'] {
  return [
    {
      column: 'sampleId',
      spec: { name: 'pl7.app/sampleId', type: 'String' as const },
    },
    {
      column: 'clonotypeKey',
      spec: {
        name: 'pl7.app/vdj/clonotypeKey',
        type: 'String' as const,
        // Mirrors bulk MiXCR's emission: chain on the axis, no receptor key.
        // SD-008 derives the receptor from this chain.
        domain: {
          'pl7.app/vdj/clonotypingRunId': 'synthetic-bulk-run',
          'pl7.app/vdj/chain': chain,
        },
        annotations: { 'pl7.app/label': 'Clonotype ID' },
      },
    },
  ];
}

/** Build the column spec list for a bulk scenario (single chain) with full coverage. */
function bulkVdjColumnsFull(): Spec['columns'] {
  // Same FR4InFrame normalisation as sc; bulk emits a single chain so no chain prefix.
  const featuresOnDisk = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4InFrame'] as const;
  const cols: Spec['columns'] = [
    {
      column: 'umiCount',
      spec: {
        name: 'pl7.app/vdj/uniqueMoleculeCount',
        valueType: 'Long' as const,
        annotations: {
          ...ANCHOR_ANNOTATIONS,
          'pl7.app/abundance/unit': 'umis',
          'pl7.app/label': 'Number of UMIs',
        },
      },
    },
  ];
  for (const f of featuresOnDisk) {
    cols.push(bulkVdjSequenceColumn(f.replace('InFrame', ''), f));
  }
  return cols;
}

/** Bulk MiXCR IGHeavy, full coverage. SD-008 derives receptor=IG from the axis chain. */
export function bulkIgHeavyScenario() {
  return {
    tsv: './assets/bulk-ig-heavy.tsv',
    fileExt: 'tsv' as const,
    spec: {
      axes: bulkAxes('IGHeavy'),
      columns: bulkVdjColumnsFull(),
      storageFormat: 'Binary' as const,
      partitionKeyLength: 0,
    } satisfies Spec,
  };
}

/** Bulk MiXCR TCRAlpha, full coverage. SD-008 derives receptor=TCRAB from the axis chain. */
export function bulkTcrAlphaScenario() {
  return {
    tsv: './assets/bulk-tcr-alpha.tsv',
    fileExt: 'tsv' as const,
    spec: {
      axes: bulkAxes('TCRAlpha'),
      columns: bulkVdjColumnsFull(),
      storageFormat: 'Binary' as const,
      partitionKeyLength: 0,
    } satisfies Spec,
  };
}

/** Bogus axis name — drives the R1a panic path. */
export function bogusAxisScenario() {
  return {
    tsv: './assets/bogus-axis.tsv',
    fileExt: 'tsv' as const,
    spec: {
      axes: [
        {
          column: 'sampleId',
          spec: { name: 'pl7.app/sampleId', type: 'String' as const },
        },
        {
          column: 'mysteryKey',
          spec: { name: 'pl7.app/notReal/mysteryKey', type: 'String' as const },
        },
      ],
      columns: [
        {
          column: 'cellCount',
          spec: {
            name: 'pl7.app/abundance',
            valueType: 'Long' as const,
            annotations: ANCHOR_ANNOTATIONS,
          },
        },
        {
          column: 'sequence',
          spec: {
            name: 'pl7.app/sequence',
            valueType: 'String' as const,
            domain: { 'pl7.app/alphabet': 'aminoacid' },
          },
        },
      ],
      storageFormat: 'Binary' as const,
      partitionKeyLength: 0,
    } satisfies Spec,
  };
}
