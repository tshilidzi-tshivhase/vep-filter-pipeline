"""Run the VEP filtering and patient-merge pipeline, one step at a time."""
import argparse
import os
import sys
import time

import patient_merge as pm
import vep_filters as vf

TOTAL_STEPS = 11
OPERATORS = {'ge': '>=', 'gt': '>', 'le': '<=', 'lt': '<', 'eq': '=='}
DEFAULT_CONSEQUENCES = [
    'frameshift_variant', 'inframe_deletion', 'inframe_insertion', 'missense_variant',
    'start_lost', 'stop_gained', 'stop_lost', 'splice_acceptor_variant', 'splice_donor_variant',
]


def say(message):
    print(message, flush=True)


def parse_predictor(text):
    parts = text.split(':')
    if len(parts) != 3 or parts[1] not in OPERATORS:
        raise argparse.ArgumentTypeError(
            f"'{text}' is not valid. Use COLUMN:OPERATOR:CUTOFF with OPERATOR one of "
            f"{', '.join(sorted(OPERATORS))} (for example CADD_PHRED:ge:20)")
    try:
        cutoff = float(parts[2])
    except ValueError:
        raise argparse.ArgumentTypeError(f"cutoff '{parts[2]}' is not a number")
    return parts[0], OPERATORS[parts[1]], cutoff


def build_parser():
    p = argparse.ArgumentParser(
        description="Filter VEP annotations and list which patients carry which mutations.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    p.add_argument('--vep', required=True, help="VEP output file")
    p.add_argument('--vcf', required=True, help="patient VCF file (.vcf or .vcf.gz)")
    p.add_argument('--canonical', required=True, help="file of canonical transcript IDs, one per line")
    p.add_argument('--genes', help="file of candidate genes, one per line (leave out to keep all genes)")
    p.add_argument('--out', default='results', help="folder for the output files")

    p.add_argument('--af-columns', nargs='+', default=['AF'],
                   help="allele frequency columns to check (AF is 1000 Genomes)")
    p.add_argument('--af-threshold', type=float, default=0.01,
                   help="remove variants more common than this")
    p.add_argument('--drop-missing-af', action='store_true',
                   help="also remove variants with no frequency value (they are kept by default)")
    p.add_argument('--consequences', nargs='+', default=DEFAULT_CONSEQUENCES,
                   help="consequence types to keep")
    p.add_argument('--clinsig', nargs='+',
                   help="keep only these exact ClinVar terms, e.g. Pathogenic Likely_pathogenic")
    p.add_argument('--clinvar-has-value', action='store_true',
                   help="keep only variants that have a ClinVar entry")
    p.add_argument('--variant-type', choices=['all', 'known', 'novel'], default='all',
                   help="known = has an existing ID, novel = has none")
    p.add_argument('--predictor', action='append', type=parse_predictor,
                   metavar='COLUMN:OPERATOR:CUTOFF',
                   help="count how many predictors pass; repeat for each predictor. "
                        "OPERATOR is ge, gt, le, lt or eq")
    p.add_argument('--show', nargs='+',
                   help="prediction columns to show in the output (default: the standard set)")
    p.add_argument('--patients', nargs='+', help="only show these patient IDs")
    p.add_argument('--positions', nargs='+', help="only show these positions, e.g. chr6:32580779")
    return p


def run_filter(number, label, df, function, *args, **kwargs):
    before = len(df)
    df = function(df, *args, **kwargs)
    say(f"[{number}/{TOTAL_STEPS}] {label}: {before:,} -> {len(df):,} variants")
    return df


def skipped(number, label, reason):
    say(f"[{number}/{TOTAL_STEPS}] {label}: skipped ({reason})")


def main():
    args = build_parser().parse_args()
    started = time.time()

    for path in [args.vep, args.vcf, args.canonical, args.genes]:
        if path and not os.path.isfile(path):
            sys.exit(f"Error: file not found: {path}")
    os.makedirs(args.out, exist_ok=True)
    say(f"Pipeline started {time.strftime('%Y-%m-%d %H:%M:%S')}")

    vep = vf.load_vep(args.vep)
    say(f"[1/{TOTAL_STEPS}] Loaded {args.vep}: {len(vep):,} rows, {vep.shape[1]} columns")

    vep = run_filter(2, f"Population frequency filter ({', '.join(args.af_columns)} <= {args.af_threshold})",
                     vep, vf.filter_population_af, args.af_columns,
                     threshold=args.af_threshold, keep_missing=not args.drop_missing_af)
    vep = run_filter(3, "Canonical transcript filter", vep, vf.filter_canonical, args.canonical)

    if args.genes:
        vep = run_filter(4, "Candidate gene filter", vep, vf.filter_genes, args.genes)
    else:
        skipped(4, "Candidate gene filter", "no --genes file given")

    vep = run_filter(5, "Consequence filter", vep, vf.filter_consequence, args.consequences)

    if args.clinsig:
        vep = run_filter(6, f"ClinVar filter ({', '.join(args.clinsig)})", vep, vf.filter_clinsig, args.clinsig)
    elif args.clinvar_has_value:
        vep = run_filter(6, "ClinVar filter (has a ClinVar entry)", vep, vf.filter_clinsig, None)
    else:
        skipped(6, "ClinVar filter", "not requested")

    if args.variant_type != 'all':
        vep = run_filter(7, f"Known/novel filter ({args.variant_type} only)", vep,
                         vf.filter_known_novel, args.variant_type)
    else:
        skipped(7, "Known/novel filter", "keeping all variants")

    if len(vep) == 0:
        sys.exit("No variants are left after filtering. Relax a filter and try again.")

    if args.predictor:
        criteria = {column: (cutoff, operator) for column, operator, cutoff in args.predictor}
        vep = vf.count_predictor_hits(vep, criteria)
        say(f"[8/{TOTAL_STEPS}] Predictor hit count added using: {', '.join(criteria)}")
    else:
        skipped(8, "Predictor hit count", "no --predictor given, so the Predictions column will show -")

    vep = pm.add_vep_keys(vep)
    vep_keys = set(zip(vep['chrom'], vep['vep_pos'], vep['vep_allele']))
    say(f"[9/{TOTAL_STEPS}] Reading genotypes from {args.vcf} ...")
    genotypes, patients, found_keys, dropped = pm.read_vcf_genotypes(args.vcf, vep_keys, log=say)
    say(f"[9/{TOTAL_STEPS}] {len(found_keys)} of {len(vep_keys)} variants have a PASS record "
        f"({len(patients)} patients); {len(dropped)} VCF records at these positions did not say PASS")

    heading_for = {source: heading for heading, source in pm.DEFAULT_PREDICTION_COLUMNS}
    prediction_columns = None
    if args.show:
        prediction_columns = [(heading_for.get(column, column), column) for column in args.show]
    table = pm.build_patient_table(vep, genotypes, patients, prediction_columns)

    if args.patients:
        unknown = sorted(set(args.patients) - set(patients))
        if unknown:
            say(f"    Warning: these patient IDs are not in the VCF: {', '.join(unknown)}")
        table = table[table['Sample'].isin(args.patients)]
    if args.positions:
        table = table[table['Position'].isin(args.positions)]

    counts = ", ".join(f"{name} {count}" for name, count in table['Zygosity'].value_counts().items())
    say(f"[10/{TOTAL_STEPS}] Patient table: {len(table):,} rows, {table['Sample'].nunique()} patients ({counts})")

    table_path = os.path.join(args.out, 'patient_mutations.tsv')
    matrix_path = os.path.join(args.out, 'genotype_matrix.tsv')
    table.to_csv(table_path, sep='\t', index=False, na_rep='-')
    pm.genotype_matrix(genotypes, patients).to_csv(matrix_path, sep='\t')
    say(f"[11/{TOTAL_STEPS}] Saved {table_path} and {matrix_path}")
    say(f"Done in {time.time() - started:.0f} seconds.")


if __name__ == '__main__':
    try:
        main()
    except ValueError as error:
        sys.exit(f"Error: {error}")

