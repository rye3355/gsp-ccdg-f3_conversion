#!/usr/bin/env python3

__author__ = "rye"

import hail as hl
import logging
import argparse
from datetime import date

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)




def get_bucket(
    data_type: str
) -> str:
    """
    Return path prefix to desired CCDG bucket.
    :param data_type: Whether data is from CCDG exomes or genomes
    :return: Path prefix to CCDG bucket
    """
    bucket = "gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/" if data_type == "genomes" else "gs://fc-secure-7e69c896-d6c0-4a4e-8490-42cb2d4fdebf/"
    return bucket




def read_in_vds(
    VDS: str
) -> hl.vds.VariantDataset:
    """
    Read in and return specified VDS. Tries to handle legacy VDS representation.
    :param VDS: Path to VDS of interest
    :return: VDS of interest
    """
    logger.info(f"Starting to read VDS...")
    # Attempt the usual reading strategy
    try:
        vds = hl.vds.read_vds(VDS)

    # If it failed, try Tim's solution for old VDS format
    except:
        logger.info(f"Usual VDS read function failed. Trying legacy read...")

        def read_old_vds(path):
            rd = hl.read_matrix_table(path + '/reference_data')
            vd = hl.read_matrix_table(path + '/variant_data')
            rd = rd.transmute_rows(ref_allele=rd.alleles[0])
            return hl.vds.VariantDataset(rd, vd)

        vds = read_old_vds(VDS)

    return vds



def prepare_vds_to_mt(
    vds: hl.vds.VariantDataset,
    chr: str
) -> hl.MatrixTable:
    """
    Filter VDS to chr, densify to MT, compute GT from LGT+LA and finally split
    multiallelic sites.
    :param vds: VDS to be used
    :param chr: String representation of chromosome to be filtered to
    :return: Densified MT filtered to chr with computed GT and split multi sites
    """
    logger.info(f"Starting to prepare VDS (filtering to {chr})...")

    # Filter VDS to chr
    vds = hl.vds.filter_chromosomes(vds, keep = chr)

    # Compute GT
    vds.variant_data = vds.variant_data.annotate_entries(GT = hl.vds.lgt_to_gt(vds.variant_data.LGT, vds.variant_data.LA))

    # Densify VDS
    mt = hl.vds.to_dense_mt(vds)

    # Split multiallelic sites
    mt = hl.split_multi_hts(mt)

    return mt




def prepare_mt(
    mt: hl.MatrixTable,
    data_type: str,
    bucket: str
) -> hl.MatrixTable:
    """
    Read in variant_qc and variant_info HTs to annotate MT with AS_VQSLOD and
    AS_lowqual fields (respectively). Also annotate is_snv.
    :param mt: MT to be annotated
    :param data_type: Whether data is from CCDG exomes or genomes
    :param bucket: Path prefix to either exomes or genomes buckets
    :return: MT annotated with AS_VQSLOD, AS_lowqual, is_snv variant fields
    """
    logger.info("Starting to annotate AS_VQSLOD, AS_lowqual, and is_snv...")

    # Read in annotation tables
    VARIANT_QC_HT = bucket + "ccdg_" + data_type + "_variant_qc_vqsr_alleleSpecificTrans_split.ht"
    variant_qc_ht = hl.read_table(VARIANT_QC_HT)
    VARIANT_INFO_HT = bucket + "ccdg_" + data_type + "_variant_info_split.ht"
    variant_info_ht = hl.read_table(VARIANT_INFO_HT)

    # Annotate AS_VQSLOD
    mt = mt.annotate_rows(AS_VQSLOD = variant_qc_ht[mt.locus, mt.alleles].info.AS_VQSLOD)

    # Annotate AS_lowqual
    mt = mt.annotate_rows(AS_lowqual = variant_info_ht[mt.locus, mt.alleles].AS_lowqual)

    # Annotate is_snv
    mt = mt.annotate_rows(is_snv = hl.is_snp(mt.alleles[0], mt.alleles[1]))

    return mt




def compute_passing(
    mt: hl.MatrixTable,
    data_type: str
) -> hl.MatrixTable:
    """
    Compute variant passing status based on Wenhan's default recommendations:
    If WGS: (AS_lowqual not True) and (AS_VQSLOD < -2.4464 | SNV) and (AS_VQSLOD < -0.3533 | INDEL)
    If WES: (AS_lowqual not True) and (AS_VQSLOD < -7.4038 | SNV) and (AS_VQSLOD < -2.8062 | INDEL)
    :param mt: MT to be annotated
    :param data_type: Whether data is from CCDG exomes or genomes
    :return: MT annotated with filter variant field
    """
    logger.info("Starting to annotate passing filter...")

    # Set different cutoffs for each data type
    snv_cutoff = -2.4464 if data_type == "genomes" else -7.4038
    indel_cutoff = -0.3533 if data_type == "genomes" else -2.8062

    # Compute PASS/FAIL status
    mt = mt.annotate_rows(filters =
            hl.if_else((~mt.AS_lowqual) &
                     (((mt.is_snv) & (mt.AS_VQSLOD <  snv_cutoff))|
                      ((~mt.is_snv) & (mt.AS_VQSLOD <  indel_cutoff))),
                      {"PASS"}, {"FAIL"}))

    return mt




def format_mt_for_vcf(
    mt: hl.MatrixTable
) -> hl.MatrixTable:
    """
    Make final formatting adjustments to mt before writing to VCF
    Columns: [only 's' is taken by export_vcf]
    info becomes INFO field
    filters becomes FILTERS field
    rsid becomes ID field
    qual becomes QUAL field
    No other row fields exported
    :param mt: MT to be annotated
    :return: MT annotated with only GT entry field and info field
    """
    logger.info("Starting to format MT for VCF...")

    # Keep only the GT entry field
    mt = mt.select_entries(mt.GT)

    # Construct the info field
    mt = mt.annotate_rows(info = hl.struct(was_split = mt.was_split,
                                                    AS_VQSLOD = mt.AS_VQSLOD,
                                                    AS_lowqual = mt.AS_lowqual,
                                                    is_snv = mt.is_snv))

    return mt



def export_VDS_to_VCF(
    vds: hl.vds.VariantDataset,
    data_type: str,
    bucket: str,
    chr: str,
    out_path: str
):
    """
    Driver function for annotating and exporting VDS to VCF.
    Filters to provided chr, lifts variant annotations over, computes PASS/FAIL,
    and finally writes to VCF.
    :param vds: VDS to be used
    :param data_type: Whether data is from CCDG exomes or genomes
    :param bucket: Path prefix to either exomes or genomes buckets
    :param chr: String representation of chromosome to be filtered to
    :param out_path: VCF output path
    """
    logger.info(f"Starting to export VDS ({chr}) to VCF...")

    # Filter VDS to chr, convert to MT, compute GT, and split multi
    mt = prepare_vds_to_mt(vds, chr)

    # Annotate MT with AS_VQSLOD, AS_lowqual, and is_snv
    mt = prepare_mt(mt, data_type, bucket)

    # Compute passing status using Wenhan's default suggestions and annotate
    mt = compute_passing(mt, data_type)

    # Last formatting for mt before writing to VCF
    mt = format_mt_for_vcf(mt)

    # Format metadata
    snv_cutoff = -2.4464 if data_type == "genomes" else -7.4038
    indel_cutoff = -0.3533 if data_type == "genomes" else -2.8062
    metadata = {'filter': {'PASS': {'Description': 'Variant passing QC'}, 'FAIL': {'Description': 'Variant failing QC'}},
                'info': {'was_split': {'Description': 'True if this variant was originally multiallelic, otherwise False.',
                                'Number': '0',
                                'Type': 'Flag'},
                        'AS_VQSLOD': {'Description': f"AS_VQSLOD score. Passing variants have score < {snv_cutoff} if SNV and < {indel_cutoff} otherwise (indel).",
                                'Number': '1',
                                'Type': 'Float'},
                        'AS_lowqual': {'Description': 'AS_lowqual classification. Passing variants have value FALSE.',
                                'Number': '1',
                                'Type': 'Flag'},
                        'is_snv': {'Description': 'Whether a variant is an SNV',
                                'Number': '0',
                                'Type': 'Flag'}
                    }
           }

    # Write to VCF
    logger.info(f"Starting to write to VCF...")
    hl.export_vcf(mt, out_path, metadata = metadata)



def main(args):
    data_type = "exomes" if args.exomes else "genomes"
    hl.init(default_reference="GRCh38",
            tmp_dir="gs://bgen-temp/tmp_dir")

    # Grab workspace bucket prefix (for either exomes or genomes)
    bucket = get_bucket(data_type)

    # Read in VDS
    vds = read_in_vds(args.vds)

    # Convert VDS
    # If no specific chromosome specified, cycle through chr1 to chr21
    if (not args.chr):
        chrs = list(range(1, 23)) + ['M', 'X', 'Y']
        for c in chrs:
            path = f"{args.out}_chr{c}.vcf.bgz"
            export_VDS_to_VCF(vds, data_type, bucket, f"chr{c}", path)
    else:
        path = f"{args.out}_chr{args.chr}.vcf.bgz"
        export_VDS_to_VCF(vds, data_type, bucket, f"chr{args.chr}", path)


    logger.info("Copying log to logging bucket...")
    hl.copy_log(f"gs://bgen-temp/log/{data_type}_export_{date.today()}_chr{args.chr}.log")





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--exomes",
        help = "Subset CCDG exomes. Otherwise, default is CCDG genomes.",
        action = "store_true",
    )
    parser.add_argument(
        "--vds",
        help = "Full path to input VDS for conversion",
        type = str,
    )
    parser.add_argument(
        "--chr",
        help="Chromosome to subset down to and export as VCF. If no input, run on chr1-chr21",
        default=0,
        type=int,
    )
    parser.add_argument(
        "--out",
        help="Output path prefix for final VCF WITHOUT file extension (.vcf.bgz will be used)",
        type=str,
    )

    args = parser.parse_args()

    main(args)
