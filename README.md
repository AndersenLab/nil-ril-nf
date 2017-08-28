# nils-nf

Align, call variants, and generate datasets for NIL sequence data

## Usage

```
# cd to directory of fastqs
nextflow run main.nf -resume --fq_folder data/test_fq<fq_directory>
```

The __NIL__ workflow works a little differently than the RIL and WI workflows. NIL sequence sets generally only need to be run once, and the resultant datasets are not mixed together as they are with wild isolate and RIAIL sequence data.

In order to process NIL data, you need to move the sequence data to a folder and create a `fq_sheet.tsv`. This file defines the fastqs that should be processed. __Note__: Unlike other workflows the fastq path provided is *relative* and not *absolute* for each fastq. This is because the NIL workflow processes fastqs located within a folder rather.

## fq_sheet.tsv

The `fq_sheet.tsv` defines the fastqs to be processed as part of the NIL workflow. It comes in the following format:

| strain   | fastq_pair_id   | library   | fastq-1-path   | fastq-2-path   |
|:-------|:-----------------------|:------------------|:-------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------|
| QX98   | QX98_GCTACGCTGCTCGAA   | GCTACGCTGCTCGAA   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX98_GCTACGCT-GCTCGAA_L003_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX98_GCTACGCT-GCTCGAA_L003_R2_001.fq.gz   |
| QX99   | QX99_AGGCAGAAGCTAATC   | AGGCAGAAGCTAATC   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L005_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L005_R2_001.fq.gz   |
| QX99   | QX99_AGGCAGAAGCTAATC   | AGGCAGAAGCTAATC   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L006_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_AGGCAGAA-GCTAATC_L006_R2_001.fq.gz   |
| QX99   | QX99_CGAGGCTGTGGCAAT   | CGAGGCTGTGGCAAT   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L003_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L003_R2_001.fq.gz   |
| QX99   | QX99_CGAGGCTGTGGCAAT   | CGAGGCTGTGGCAAT   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L004_R1_001.fq.gz   | /projects/b1059/data/fastq/RIL/dna/processed/151009_D00422_0262_BC7NJ0ANXX-ECA/QX99_CGAGGCTG-TGGCAAT_L004_R2_001.fq.gz   |

This file can be generated using a script that *looks* like this, but not that you may have to modify it depending on the naming convention used for the provided fastqs.

```
ls -1 *.gz | grep 'R1' | awk -v pwd=`pwd` '{
                            split($0, a, "_");
                            SM = a[1]; 
                            gsub("-","",a[2]); // LB
                            ID = a[1] "_" a[2];
                            fq2 = $1; gsub("_R1", "_R2", fq2);
                            print SM "\t" ID "\t" a[2] "\t" pwd "/" $1 "\t" pwd "/" fq2;
                            }' > fq_sheet.tsv
```

## Important Notes

* NIL sequencing uses low coverage by design. No pre-processing (ie trimming) takes place because variants will be called at specific positions and low quality/adapter contamination are unlikely to be problematic.

## Generate CB4856 Sitelist

The `CB4856.20160408.sitelist.tsv.gz` file is probably fine forever. It was generated with the following command:

```
bcftools view --samples CB4856,N2 -m 2 -M 2 WI.20160408.filtered.vcf.gz | \
vk filter ALT --min=1 - | \
vk filter REF --min=1 - | \
bcftools view --samples CB4856 - | \
bcftools query --include 'FORMAT/GT == "1/1"' -f '%CHROM\t%POS\t%REF,%ALT\n' | \
bgzip -c > CB4856.20160408.sitelist.tsv.gz && tabix -s 1 -b 2 -e 2 CB4856.20160408.sitelist.tsv.gz
```

## Output

The output directory looks like this:

```
├── concordance
│   └── fq_concordance.tsv
├── cross_object
│   ├── breakpoint_sites.tsv.gz
│   ├── cross_obj.Rdata
│   ├── cross_obj_geno.tsv
│   ├── cross_obj_pheno.tsv
│   └── cross_obj_strains.tsv
├── duplicates
│   └── bam_duplicates.tsv
├── fq
│   ├── fq_bam_idxstats.tsv
│   ├── fq_bam_stats.tsv
│   ├── fq_coverage.full.tsv
│   └── fq_coverage.tsv
├── fq_ril_sheet.tsv
└── hmm
    ├── gt_hmm.png
    ├── gt_hmm.svg
    └── gt_hmm.tsv
```

### concordance/
 
__fq_concordance.tsv__ - Contains the concordance among fastqs belonging to the same strain. 

* `a` - First fastq
* `b` - Second fastq
* `concordant_sites` - Number of sites concordant between `a` and `b`
* `total_sites` - The total number of sites called in both `a` and `b`
* `concordance` - frequency concordant.
* `SM` - Strain

| a                     | b                     |   concordant_sites |   total_sites |   concordance | SM    |
|:----------------------|:----------------------|-------------------:|--------------:|--------------:|:------|
| QX204_CAGAGAGGCCGGATA | EA-G02_TTAACTC_L001   |             148794 |        151385 |      0.982885 | QX204 |
| EA-G02_TTAACTC_L001   | QX204_CAGAGAGGCCGGATA |             148794 |        151385 |      0.982885 | QX204 |

### cross_object/

* __breakpoint_sites.tsv.gz__ - A gzipped file of sites flanking loci where recombination has occured.
* __cross_obj.Rdata__ - A cross object usable with `qtl`
* __cross_obj_geno.tsv__ - Genotypes file used by the R `qtl` package.
* __cross_obj_pheno.tsv__ - Phenotypes file used by the R `qtl` package.
* __cross_obj_strains.tsv__ - List of strains within the cross object.

The script used to generate the cross object is called `generate_cross_object.R` and is located in the root of this repo.

### duplicates/

__bam_duplicates.tsv__ - A summary of duplicate reads from aligned bams.

### fq/

* __fq_bam_idxstats.tsv__
* __fq_bam_stats.tsv__
* __fq_coverage.full.tsv__
* __fq_coverage.tsv__

### SM/

* __SM_bam_idxstats.tsv__
* __SM_bam_stats.tsv__
* __SM_coverage.full.tsv__
* __SM_coverage.tsv__

### hmm/

* __gt_hmm.(png/svg)__ - Image depicting RILs.
* __gt_hmm.tsv__ - Long form genotypes file.

### plots/

__coverage_comparison.(png/svg)__

__duplicates.(png/svg)__

__unmapped_reads.(png/svg)__

### vcf/

* __gt_hmm.tsv__ - Haplotypes defined by region with associated information. 
* __gt_hmm_fill.tsv__ - Same as above, but using `--infill` and `--endfill` with VCF-Kit. For more information, see [VCF-Kit Documentation](http://vcf-kit.readthedocs.io/en/latest/)
* __RIL.filter.vcf.gz__ - A VCF of filtered genotypes. 
* __RIL.filtered.stats.txt__ - Summary of filtered genotypes. Generated by `bcftools stats RIL.filter.vcf.gz`
* __RIL.hmm.vcf.gz__ - The RIL VCF as output by VCF-Kit; HMM applied to determine genotypes.
* __union_vcfs.txt__ - A list of VCFs that were merged to generate RIL.filter.vcf.gz
