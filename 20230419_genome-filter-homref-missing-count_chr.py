import hail as hl
chr = "chr22"


# Initialize Hail
TMP_DIR = 'gs://gsp-ccdg-f3-tmp-4day/'
REF = 'GRCh38'
hl.init(tmp_dir = TMP_DIR, default_reference = REF)

# Read in annotation tables
VARIANT_QC_HT = "gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/ccdg_genomes_variant_qc_vqsr_alleleSpecificTrans_split.ht"
variant_qc_ht = hl.read_table(VARIANT_QC_HT)

VARIANT_INFO_HT = "gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/ccdg_genomes_variant_info_split.ht"
variant_info_ht = hl.read_table(VARIANT_INFO_HT)

# Read in genomes VDS
VDS = "gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/ccdg_genome_136k.vds"
vds = hl.vds.read_vds(VDS)

# Filter to chr
print(f"Filtering to {chr}...")
vds = hl.vds.filter_chromosomes(vds, keep = chr)

# Compute GT
print("Computing GT...")
vds.variant_data = vds.variant_data.annotate_entries(GT = hl.vds.lgt_to_gt(vds.variant_data.LGT, vds.variant_data.LA))
vds.variant_data = vds.variant_data.drop(vds.variant_data.LGT)

# Densify VDS
print("Densifying...")
mt = hl.vds.to_dense_mt(vds)

# Split multiallelic
print("Splitting multi...")
mt = hl.split_multi_hts(mt)

# Compute entry filter stats
print("Computing entry filter stats...")
mt = mt.compute_entry_filter_stats()

# Count bad entries for each row (use to filter out later)
## IF hom_ref or NA (due to entry filtering)
print("Counting number of bad entries per row...")
mt = mt.annotate_rows(num_bad = hl.agg.count_where(mt.GT.is_hom_ref()) + mt.entry_stats_row.n_filtered)

# Focus on only rows
mt_rows = mt.rows()

# Annotate info.AS_VQSLOD, AS_lowqual, and is_snv
mt_rows = mt_rows.annotate(AS_VQSLOD = variant_qc_ht[mt_rows.locus, mt_rows.alleles].info.AS_VQSLOD)
mt_rows = mt_rows.annotate(AS_lowqual = variant_info_ht[mt_rows.locus, mt_rows.alleles].AS_lowqual)
mt_rows = mt_rows.annotate(is_snv = hl.is_snp(mt_rows.alleles[0], mt_rows.alleles[1]))

data_type = "genomes"
snv_cutoff = -2.4464 if data_type == "genomes" else -7.4038
indel_cutoff = -0.3533 if data_type == "genomes" else -2.8062

# Compute PASS/FAIL status
mt_rows = mt_rows.annotate(filters =
        hl.if_else((~hl.is_missing(mt_rows.AS_VQSLOD)) & (~mt_rows.AS_lowqual) & 
                    (((mt_rows.is_snv) & (mt_rows.AS_VQSLOD >=  snv_cutoff))|
                    ((~mt_rows.is_snv) & (mt_rows.AS_VQSLOD >=  indel_cutoff))),
                    {"PASS"}, {"FAIL"}))

# Checkpoint
mt_rows.write(f"gs://bgen-temp/genome-variant-data_rows-annotated-bad_{chr}.ht", overwrite = True)




