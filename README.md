# gsp-ccdg-f3_conversion

Script to annotate and convert gsp-ccdg-f3 VDS files to standard VCF format. 

Variant filter status annotated with Wenhan's suggested variant qc metrics, and additional [was_split, AS_VQSLOD, AS_lowqual, is_snv] fields also included.

A variant is labelled PASS if:
* (info.AS_lowqual != True) & (info.AS_VQSLOD >= [cutoff] | is.snv) & (info.AS_VQSLOD >= [cutoff] | ~is.snv), depending on genomes or exomes data
* See here for more detailed information: https://docs.google.com/document/d/1kG7Hf2QhDKmhhRfmmhnBkC1HLX0mj2jZ38Z4Jp427qs/edit?usp=sharing


Annotations for AS_lowqual and AS_VQSLOD are pulled from Wenhan's annotation HailTables:
* gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/ccdg_genomes_variant_qc_vqsr_alleleSpecificTrans_split.ht
* gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/ccdg_genomes_variant_info_split.ht
* gs://fc-secure-7e69c896-d6c0-4a4e-8490-42cb2d4fdebf/ccdg_exomes_variant_qc_vqsr_alleleSpecificTrans_split.ht
* gs://fc-secure-7e69c896-d6c0-4a4e-8490-42cb2d4fdebf/ccdg_exomes_variant_info_split.ht


Example of how to run:
* First, spin up workers on dataproc using Hail:
  * hailctl dataproc start [cluster-name] --num-workers 6 --num-secondary-workers 120 
  * Can increase number of preemptibles (--num-secondary-workers) and non-preemptibles (--num-workers)
  * From Tim: "It's not bad practice to make sure your ratio of preemptible to non-preemptible isn't larger than around 20:1"
* Next, submit script with relevant flags
  * Running on ONLY CHROMOSOME 21 of exome dataset:
     * hailctl dataproc submit [cluster-name] gsp-ccdg-f3_vds-to-vcf.py --exomes --vds gs://fc-secure-7e69c896-d6c0-4a4e-8490-42cb2d4fdebf/ccdg_exome_203k.vds/ --chr 21 --out gs://bgen-temp/ccdg_exome_203k/ccdg_exome_203k
  * Running on ALL CHROMOSOMES of genome dataset:
     * hailctl dataproc submit [cluster-name] gsp-ccdg-f3_vds-to-vcf.py --vds gs://fc-secure-9e3357c0-389c-41d7-94ee-56673db6b75f/ccdg_genome_136k.vds/ --out gs://bgen-temp/ccdg_genome_136k/ccdg_genome_136k
  * Running on ALL CHROMOSOMES of exome dataset:
     * hailctl dataproc submit [cluster-name] gsp-ccdg-f3_vds-to-vcf.py --exomes --vds gs://fc-secure-7e69c896-d6c0-4a4e-8490-42cb2d4fdebf/ccdg_exome_203k.vds/ --out gs://bgen-temp/ccdg_exome_203k/ccdg_exome_203k
  * Use --help for more info
