   VEP filtering and patient-merge pipeline

   What it does
     Filters VEP annotations (population frequency, canonical transcripts, candidate
     genes, consequence type, ClinVar, known/novel), counts how many prediction
     tools pass user-chosen cutoffs, then reads a patient VCF and lists which
     patients carry which mutations, with homozygous/heterozygous status.

   Files
     run_pipeline.py    the command you run (options + progress messages)
     vep_filters.py     the filter functions
     patient_merge.py   reads the VCF and builds the patient table
     run_pipeline.pbs   example script for submitting the pipeline to the queue

   Requirements
     python3 with pandas

   Inputs
     --vep        VEP output file (tab-delimited; the header line starts with #)
     --vcf        patient VCF (.vcf or .vcf.gz). Only records with FILTER = PASS
                  are used, and only the genotype (GT) field.
     --canonical  canonical transcript IDs, one per line (version numbers ignored)
     --genes      candidate genes, one per line, exact match (optional)

   Example
     python3 run_pipeline.py --vep samples.vepout --vcf patients.vcf \
         --canonical canonical_transcripts.txt --genes candidate_genes.txt --out results

   Options (python3 run_pipeline.py --help lists them all with their defaults)
     --af-columns COL ...      allele frequency columns to check (default AF = 1000 Genomes)
     --af-threshold X          remove variants more common than X (default 0.01)
     --drop-missing-af         also remove variants with no frequency value
                               (they are kept by default)
     --consequences TYPE ...   consequence types to keep (default: the nine agreed types)
     --clinsig TERM ...        keep only these exact ClinVar terms,
                               e.g. Pathogenic Likely_pathogenic
     --clinvar-has-value       keep only variants that have a ClinVar entry
     --variant-type all|known|novel
                               known = has an existing ID, novel = has none
     --predictor COLUMN:OP:CUTOFF
                               count how many predictors pass; repeat the option once
                               per predictor. OP is ge, gt, le, lt or eq. ge and le
                               include the cutoff. Numeric score columns only
                               (not the SIFT/PolyPhen text columns). Choosing the
                               right direction is the user's responsibility: CADD is
                               higher = more damaging (ge), PROVEAN is lower = more
                               damaging (le).
     --show COL ...            prediction columns to show in the output
     --patients ID ...         show only these patients
     --positions POS ...       show only these positions, e.g. chr6:32580779
     --out FOLDER              output folder (default: results)

   Outputs (in the output folder)
     patient_mutations.tsv   one row per patient per mutation they carry. Missing
                             values are written as -
     genotype_matrix.tsv     mutations as rows, patients as columns. Each cell is
                             Reference, Heterozygous, Homozygous, Hemizygous or Missing

   Notes
     - The Predictions column is the number of --predictor criteria a variant meets.
       Newer prediction tools have many missing values, and a missing score counts
       as "not met". With no --predictor, the column shows -.
     - Allele frequency is the highest value across the --af-columns used for the
       filter.
     - Heterozygous = one copy of the variant, Homozygous = two copies,
       Hemizygous = a single allele (for example the X chromosome in males).
