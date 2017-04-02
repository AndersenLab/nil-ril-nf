#!/usr/bin/bash
# Generate fq_sheet for 150225_D00422_0165_BHGMN3ADXX
prefix=/projects/b1059/projects/kammenga-nils-nf/fastq
ls -1 ${prefix} | grep 'R1' | awk -v prefix=${prefix} '{
                            split($0, a, "_");
                            SM = a[1]; 
                            gsub("-","",a[2]); // LB
                            ID = a[1] "_" a[2];
                            fq2 = $1; gsub("_R1_", "_R2_", fq2);
                            print SM "\t" ID "\t" a[2] "\t" prefix "/" $1 "\t" prefix "/" fq2;
                            }' | sort -k1,1 > fq_nil_sheet.tsv
