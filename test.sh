#!/usr/bin/env bash

curl --version >/dev/null 2>&1 || { echo >&2 "I require curl, but it's not installed. Aborting."; exit 1; }
tar --version >/dev/null 2>&1 || { echo >&2 "I require tar, but it's not installed. Aborting."; exit 1; }
docker -v >/dev/null 2>&1 || { echo >&2 "I require docker, but it's not installed. Visit https://www.docker.com/products/overview#/install_the_platform  ."; exit 1; }
nextflow -v >/dev/null 2>&1 || { echo >&2 "I require nextflow, but it's not installed. If you hava Java, run 'curl -fsSL get.nextflow.io | bash'. If not, install Java."; exit 1; }

run_name="Test NIL-nf: "$(date +%s)

# Download reference genome
curl https://storage.googleapis.com/elegansvariation.org/genome/WS245/WS245.tar.gz > WS245.tar.gz
tar -xvzf WS245.tar.gz

cmd="""nextflow run main.nf  \
                   -with-docker andersenlab/nil-ril-nf \
                   -profile travis \
                   --reference WS245/WS245.fa.gz
                   --fqs ${TRAVIS_BUILD_DIR}/test_data/fq_sheet.tsv \
                   --vcf ${TRAVIS_BUILD_DIR}/test_data/N2_CB.simple.vcf.gz \
                   --cores 1 \
                   --tmpdir ${TRAVIS_BUILD_DIR}/tmp \
                   -resume """

echo "Starting nextflow... Command:"
echo $cmd
echo "-----"
eval $cmd
