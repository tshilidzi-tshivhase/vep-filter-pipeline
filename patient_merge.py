"""Patient genotype functions (imported by run_pipeline.py)."""
import gzip
import pandas as pd

# (heading in the output file, column it comes from), in the supervisor's order
OUTPUT_COLUMNS = [
    ('Sample', 'patient'),
    ('Gene symbol', 'SYMBOL'),
    ('Position', 'Position'),
    ('Ref', 'ref'),
    ('Alt', 'alt'),
    ('Zygosity', 'zygosity'),
    ('Allele frequency', 'max_pop_af'),
    ('Transcript', 'Feature'),
    ('Consequence', 'Consequence'),
    ('Predictions', 'predictor_hits_count'),
    ('HGVSc', 'HGVSc'),
    ('HGVSp', 'HGVSp'),
    ('Amino Acid', 'Amino_acids'),
    ('Existing Variant', 'Existing_variation'),
    ('Clinvar Clinical Significance', 'clinvar_clnsig'),
]

# shown unless the user asks for other prediction columns
DEFAULT_PREDICTION_COLUMNS = [
    ('SIFT', 'SIFT'),
    ('PolyPhen', 'PolyPhen'),
    ('BayesDelNoAF Score', 'BayesDel_noAF_score'),
    ('CADDphred', 'CADD_PHRED'),
    ('Primate AI Score', 'PrimateAI_score'),
    ('REVEL Score', 'REVEL_score'),
    ('MutPred2 Score', 'MutPred2_score'),
    ('ESM1b Score', 'ESM1b_score'),
    ('VarityR Score', 'VARITY_R_score'),
]

INFO_COLS = ['chrom', 'vcf_pos', 'ref', 'alt', 'vep_pos', 'vep_allele']
CARRIER_STATES = ['Heterozygous', 'Homozygous', 'Hemizygous']


def add_vep_keys(vep):
    vep = vep.copy()
    parts = vep['Location'].str.split(':', expand=True)
    vep['chrom'] = parts[0].str.replace('chr', '')
    vep['vep_pos'] = parts[1].str.split('-', expand=True)[0].astype(int)
    vep['vep_allele'] = vep['Allele'].fillna('-')
    return vep


def vep_allele(ref, alt):
    if ref[0] == alt[0] and (len(ref) > 1 or len(alt) > 1):
        alt = alt[1:]
    return alt if alt != '' else '-'


def classify(cell, alt_index):
    gt = str(cell).split(':')[0].replace('|', '/')
    alleles = gt.split('/')
    if '.' in alleles:
        return 'Missing'
    copies = alleles.count(str(alt_index))
    if copies == 0:
        return 'Reference'
    if len(alleles) == 1:
        return 'Hemizygous'
    if copies == 2:
        return 'Homozygous'
    return 'Heterozygous'


def read_vcf_genotypes(vcf_file, vep_keys, log=print):
    wanted = set((c, p) for c, p, _ in vep_keys)
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for header_idx, line in enumerate(f):
            if line.startswith('#CHROM'):
                break
        else:
            raise ValueError("No #CHROM header line found - is this a VCF file?")

    rows = []
    dropped = []
    found_keys = set()
    samples = []
    scanned = 0

    for chunk in pd.read_csv(vcf_file, sep='\t', skiprows=header_idx, dtype=str, chunksize=50000):
        scanned += len(chunk)
        chunk = chunk.rename(columns={'#CHROM': 'CHROM'})
        samples = list(chunk.columns[9:])
        chunk['chrom'] = chunk['CHROM'].str.replace('chr', '')
        chunk['pos'] = chunk['POS'].astype(int)

        in_wanted = [(c, p) in wanted for c, p in zip(chunk['chrom'], chunk['pos'])]
        chunk = chunk[in_wanted]

        failed = chunk[chunk['FILTER'] != 'PASS']
        for _, r in failed.iterrows():
            dropped.append((r['chrom'], r['pos'], r['FILTER']))
        chunk = chunk[chunk['FILTER'] == 'PASS']

        for _, row in chunk.iterrows():
            for alt_index, alt in enumerate(row['ALT'].split(','), start=1):
                if alt == '*' or alt.startswith('<'):
                    continue
                key = (row['chrom'], row['pos'], vep_allele(row['REF'], alt))
                if key not in vep_keys:
                    continue
                found_keys.add(key)
                entry = {'chrom': row['chrom'], 'vcf_pos': row['pos'], 'ref': row['REF'],
                         'alt': alt, 'vep_pos': key[1], 'vep_allele': key[2]}
                for sample in samples:
                    entry[sample] = classify(row[sample], alt_index)
                rows.append(entry)

        log(f"    scanned {scanned:,} VCF records; {len(rows)} matching variants so far")

    if not rows:
        raise ValueError("No PASS variants in the VCF matched the filtered VEP variants")

    genotypes = pd.DataFrame(rows)
    genotypes['mutation'] = (genotypes['chrom'] + ':' + genotypes['vcf_pos'].astype(str)
                             + ' ' + genotypes['ref'] + '>' + genotypes['alt'])
    return genotypes, samples, found_keys, dropped


def genotype_matrix(genotypes, patient_cols):
    return genotypes.set_index('mutation')[patient_cols]


def build_patient_table(vep, genotypes, patient_cols, prediction_columns=None):
    if prediction_columns is None:
        prediction_columns = DEFAULT_PREDICTION_COLUMNS
    if 'predictor_hits_count' not in vep.columns:
        vep = vep.assign(predictor_hits_count=float('nan'))

    wanted_cols = OUTPUT_COLUMNS + list(prediction_columns)
    from_carriers = ['patient', 'Position', 'ref', 'alt', 'zygosity']
    vep_cols = []
    for _, source in wanted_cols:
        if source not in from_carriers and source not in vep_cols:
            vep_cols.append(source)
    missing = [c for c in vep_cols if c not in vep.columns]
    if missing:
        raise ValueError(f"Columns not found in the VEP file: {missing}")

    long_table = genotypes.melt(id_vars=['mutation'] + INFO_COLS, value_vars=patient_cols,
                                var_name='patient', value_name='zygosity')
    carriers = long_table[long_table['zygosity'].isin(CARRIER_STATES)].copy()
    carriers['Position'] = 'chr' + carriers['chrom'] + ':' + carriers['vcf_pos'].astype(str)

    annotation = vep[['chrom', 'vep_pos', 'vep_allele'] + vep_cols]
    annotated = carriers.merge(annotation, how='inner', on=['chrom', 'vep_pos', 'vep_allele'])

    annotated['chrom_num'] = pd.to_numeric(annotated['chrom'], errors='coerce').fillna(99)
    annotated = annotated.sort_values(['patient', 'chrom_num', 'vcf_pos'])

    output = annotated[[source for _, source in wanted_cols]].copy()
    output.columns = [heading for heading, _ in wanted_cols]
    return output


def mutations_of_patient(table, patient_id):
    return table[table['Sample'] == patient_id]


def patients_with_mutation(table, chrom, pos):
    return table[table['Position'] == f"chr{chrom}:{pos}"]
