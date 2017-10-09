FROM linuxbrew/linuxbrew:1.3.1

# Install homebrew files
RUN brew install gcc
RUN brew tap homebrew/science \
    && brew install \
            bwa \
            samtools \
            bcftools \
            bedtools \
            fastqc \
            nextflow \
            sambamba \
            vcflib \
            vcftools \
            python2 \
            R

RUN pip2 install numpy cython
RUN pip2 install vcf-kit
RUN pip2 install https://github.com/AndersenLab/bam-toolbox/archive/0.0.3.tar.gz

ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

# Install R Packages v2
RUN Rscript -e 'install.packages(c("tidyverse", "cowplot"))'