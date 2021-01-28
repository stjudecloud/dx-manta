<!-- dx-header -->
# Manta SV Caller

DNAnexus platform implementation of Illumina's Manta structural variant caller.
Can perform single sample or paired analysis.

See the Manta User Guide (https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md)
for details on Manta.

<!-- /dx-header -->

## Inputs

  * `normal_bam` (Required): One or more germline BAM files. Required for
    either germline only or paired analysis. If paired, only the first file
    will be used.

  * `normal_bam_index`: Index file(s) for germline BAM files. If none are
    provided, the app will generate them.

  * `tumor_bam`: Tumor BAM file. If provided, app will assume tumor normal
    analysis.

  * `tumor_bam_index`: Tumor BAM index file. If one is not provided, the app
    will generate it.

  * `reference_fasta` (Required): Reference FASTA file. If not provided, a
    default hg38 reference FASTA will be used. If running hg19 samples, this
    input is required or the app will error.

  * `reference_fasta_index`: Reference FASTA file index. If not provided, an
    index will be generated.

  * `output_contig`: Boolean flag for the 'outputContig' option on manta.
    [default: true]

  * `call_regions`: Boolean flag for the 'callRegions' option on manta. Only
    works on hg38 samples. [default: true]

  * `output_prefix`: Prefix to add to output files. 


## Outputs

  * `variant_outputs`: Outputs under the 'variant' folder. It includes:
                        candidateSV.vcf.gz
                        candidateSV.vcf.gz.tbi
                        diploidSV.vcf.gz
                        diploidSV.vcf.gz.tbi
                        candidateSmallIndels.vcf.gz
                        candidateSmallIndels.vcf.gz.tbi
                        somaticSV.vcf.gz
                        somaticSV.vcf.gz.tbi


  * `statistics_outputs`: Outputs under the 'statistics' folder. 
                          alignmentStatsSummary.txt
                          svLocusGraphStats.tsv
                          svCandidateGenerationStats.xml
                          svCandidateGenerationStats.tsv 

## Tips

* Instance Types - Manta can parallelize across multiple cores, therefore it
  is preferable to choose an instance type with many cores. The other key
  factor in determining instance type is the storage required to store all the
  input and reference files. We recommend `mem1_ssd1_x16` or bigger instances.

* Runtime - To reduce runtime and costs, make sure to include a reference FASTA
  file and all index files. The app will need the input files to be indexed
  before running.

* Hg19 samples - To run hg19 samples, specify the `reference_fasta` and
  `reference_fasta_index` inputs as the default references are GRCh38.
  Call Regions option should be set to `false`.

## Troubleshooting

### Storage

Ensure the instance type has enough storage space for all the inputs
and reference files. The reference fasta is about 3GB in size.

Typically for single WGS germline analysis, we recommend `azure:mem1_ssd1_x8`
or `azure:mem1_ssd1_x16` instances. For tumor-normal paired analysis, we
recommend `azure:mem2_ssd1_x16`.

### Input mismatches

If the index files are not correct or if the germline, tumor samples aren't a
pair, then the app will error possibly with a vague error message. If the app
has errors where it cannot read the input bams, ensure that all the indices are
correct. It is also recommended to input the correct reference FASTA to the app.
