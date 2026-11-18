"""Filter functions for the VEP annotation table (imported by run_pipeline.py)."""
import operator as op
import pandas as pd


def load_vep(path):
    df = pd.read_csv(path, sep="\t", na_values=["-", "."], low_memory=False,
                     float_precision='round_trip')    
    df.rename(columns={'#Uploaded_variation': 'Uploaded_variation'}, inplace=True)
    df['Allele'] = df['Allele'].fillna('-')
    return df


def filter_population_af(df, af_columns, threshold=0.01, keep_missing=True):
    missing = [c for c in af_columns if c not in df.columns]
    if missing:
        raise ValueError(f"Allele frequency columns not found: {missing}")
    df = df.copy()
    df[af_columns] = df[af_columns].apply(pd.to_numeric, errors='coerce')
    df['max_pop_af'] = df[af_columns].max(axis=1, skipna=True)
    keep = df['max_pop_af'] <= threshold
    if keep_missing:
        keep = keep | df['max_pop_af'].isna()
    return df[keep]


def filter_canonical(df, transcripts_file):
    with open(transcripts_file) as f:
        canonical = set(line.strip() for line in f if line.strip())
    base_id = df['Feature'].str.split('.').str[0]
    return df[base_id.isin(canonical)]


def filter_genes(df, genes_file, gene_column='SYMBOL'):
    with open(genes_file) as f:
        genes = set(line.strip() for line in f if line.strip())
    return df[df[gene_column].isin(genes)]


def filter_consequence(df, consequence_list, column='Consequence'):
    terms = df[column].str.split(',')
    keep = terms.apply(
        lambda t: any(x in consequence_list for x in t) if isinstance(t, list) else False
    )
    return df[keep]


def filter_clinsig(df, clinsig_terms=None, require_non_empty=True, column='clinvar_clnsig'):
    series = df[column]
    keep = series.notna() if require_non_empty else pd.Series(True, index=df.index)
    if clinsig_terms:
        terms = series.str.split('/')
        has_term = terms.apply(
            lambda t: any(x in clinsig_terms for x in t) if isinstance(t, list) else False
        )
        keep = keep & has_term
    return df[keep]


def filter_known_novel(df, which):
    if which == 'known':
        return df[df['Existing_variation'].notna()]
    if which == 'novel':
        return df[df['Existing_variation'].isna()]
    return df


def count_predictor_hits(df, criteria, new_column='predictor_hits_count'):
    ops = {'>=': op.ge, '>': op.gt, '<=': op.le, '<': op.lt, '==': op.eq}
    df = df.copy()
    hits = pd.Series(0, index=df.index)
    available = pd.Series(0, index=df.index)
    for column, (cutoff, operator_str) in criteria.items():
        if column not in df.columns:
            raise ValueError(f"Predictor column not found: {column}")
        if operator_str not in ops:
            raise ValueError(f"Unknown operator '{operator_str}' (use one of {list(ops)})")
        values = pd.to_numeric(df[column], errors='coerce')
        hits += ops[operator_str](values, cutoff).astype(int)
        available += values.notna().astype(int)
    df[new_column] = hits
    df[new_column + '_available'] = available
    return df
