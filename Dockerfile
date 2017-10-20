FROM linuxbrew/linuxbrew:1.3.1

# Install homebrew files
RUN brew install gcc
RUN brew install https://raw.githubusercontent.com/Linuxbrew/homebrew-core/043fb1f50af078db481b971d36c605f0dcf72ccd/Formula/jdk.rb
RUN brew tap homebrew/science \
    && brew install \
            bwa \
            samtools \
            bcftools \
            bedtools \
            nextflow \
            sambamba \
            vcflib \
            vcftools \
            python2 \
            R \
            pigz

RUN brew install fastqc --ignore-dependencies

RUN pip2 install numpy cython
RUN pip2 install vcf-kit
RUN pip2 install https://github.com/AndersenLab/bam-toolbox/archive/0.0.3.tar.gz

ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

# Install R Packages v2
RUN Rscript -e 'install.packages(c("tidyverse", "cowplot"))'

# Link python
RUN ln /home/linuxbrew/.linuxbrew/bin/python2 /home/linuxbrew/.linuxbrew/bin/python

# Build Olson TZ database
RUN mkdir -p /home/linuxbrew/zone_info \
       && cd /home/linuxbrew/zone_info \
       && curl https://github.com/danielecook/danielecook.github.io/raw/master/downloads/zone_info.zip > /home/linuxbrew/zone_info/zone_info.zip \
       && echo "export TZDIR=/home/linuxbrew/zone_info" >> ~/.bash_profile