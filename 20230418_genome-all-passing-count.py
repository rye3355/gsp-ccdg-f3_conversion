import hail as hl

# Initialize Hail
TMP_DIR = 'gs://gsp-ccdg-f3-tmp-4day/'
REF = 'GRCh38'
hl.init(tmp_dir = TMP_DIR, default_reference = REF)

# Read in annotation tables
VARIANT_QC_HT = "gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/ccdg_genomes_variant_qc_vqsr_alleleSpecificTrans_split.ht"
variant_qc_ht = hl.read_table(VARIANT_QC_HT)

VARIANT_INFO_HT = "gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/ccdg_genomes_variant_info_split.ht"
variant_info_ht = hl.read_table(VARIANT_INFO_HT)

# Read in MT of variant data
mt = "gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/ccdg_genome_136k.vds/variant_data"
mt = hl.read_matrix_table(mt)

# Grab only variant data rows
mt_rows = mt.rows()

# Split multi
mt_rows = hl.split_multi_hts(mt_rows)

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
mt_rows.write("gs://bgen-temp/genome-variant-data_rows-annotated.ht", overwrite = True)

