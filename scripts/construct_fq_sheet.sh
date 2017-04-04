#!/usr/bin/bash
# Generate fq_sheet for CB/N2 Strains
prefix=/projects/b1059/projects/kammenga-nils-nf/fastq
ls -1 ${prefix} | grep 'R1' | awk -v prefix=${prefix} '{
                            split($0, a, "_");
                            SM = a[1]; 
                            gsub("-","",a[2]); // LB
                            ID = a[1] "_" a[2];
                            fq2 = $1; gsub("_R1", "_R2", fq2);
                            print SM "\t" ID "\t" a[2] "\t" prefix "/" $1 "\t" prefix "/" fq2;
                            }' > fq_temp.tsv

# Generate fq_sheet for Alternative Introgressions
ls -1 ${prefix} | grep 'ECA23*' | grep '_1P' | awk -v prefix=${prefix} '{
                            split($0, a, "-");
                            split($0, b, "_");
                            SM = a[2]; 
                            gsub("-","",a[2]); // LB
                            ID = a[1] "_" a[2] "_" b[2];
                            LB = b[5] "_" b[6];
                            gsub("index","",LB); gsub("_", "", LB);
                            fq2 = $1; gsub("_1P", "_2P", fq2);
                            print SM "\t" ID "\t" LB "\t" prefix "/" $1 "\t" prefix "/" fq2;
                            }' >> fq_temp.tsv

cat fq_temp.tsv | sort -k1,1n > ../fq_nil_sheet.tsv
rm fq_temp.tsv