#!/usr/bin/env bash

script_path="../main.nf"
data_path="../test_data"

curl --version >/dev/null 2>&1 || { echo >&2 "I require curl, but it's not installed. Aborting."; exit 1; }
tar --version >/dev/null 2>&1 || { echo >&2 "I require tar, but it's not installed. Aborting."; exit 1; }
docker -v >/dev/null 2>&1 || { echo >&2 "I require docker, but it's not installed. Visit https://www.docker.com/products/overview#/install_the_platform  ."; exit 1; }
nextflow -v >/dev/null 2>&1 || { echo >&2 "I require nextflow, but it's not installed. If you hava Java, run 'curl -fsSL get.nextflow.io | bash'. If not, install Java."; exit 1; }

run_name="Test NIL-nf: "$(date +%s)

# Download reference genome
curl andersen/genome/c_elegans/WS245 > WS245.tar.gz
tar -xvzf WS245.tar.gz


cmd="""nextflow run ../main.nf -resume \
                               -name \"$run_name\" \
                               --reference WS245.fa.gz
                               --fqs ../test_data/fq_sheet.tsv \
                               --with-docker andersenlab/nil-ril \
                               --VCF ../test_data/N2_CB.simple.vcf.gz \
                               --cores 1 """

echo "Starting nextflow... Command:"
echo $cmd
echo "-----"
eval $cmd
