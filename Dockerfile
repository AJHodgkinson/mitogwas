# Use CentOS as the base image
FROM centos:latest

RUN cd /etc/yum.repos.d/
RUN sed -i 's/mirrorlist/#mirrorlist/g' /etc/yum.repos.d/CentOS-*
RUN sed -i 's|#baseurl=http://mirror.centos.org|baseurl=http://vault.centos.org|g' /etc/yum.repos.d/CentOS-*

#Load dependencies

RUN dnf install -y redhat-rpm-config
RUN yum -y update
RUN  yum -y --enablerepo=extras install epel-release

RUN yum install -y \
git \
python2 \
wget \
epel-release \
python2-pip \
gcc \ 
python2-devel \
make \
zlib-devel \
gcc-c++ \
bzip2 \
bzip2-devel \
ncurses-devel \
xz-devel \
perl-Env \
java-devel \
perl-core \
gsl \
gsl-devel

RUN yum -y install dnf-plugins-core && \
    yum config-manager --set-enabled powertools

RUN yum -y install R

RUN wget -O- https://cpanmin.us | perl - App::cpanminus
RUN cpanm --mirror-only --mirror https://cpan.metacpan.org/ Statistics::R

RUN echo 'options(repos = list(CRAN = "https://cloud.r-project.org/"))' > /usr/lib64/R/etc/Rprofile.site
RUN Rscript -e "install.packages('calibrate')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('gridExtra')"

RUN Rscript -e "install.packages(c('BiocManager', 'remotes'))"
RUN Rscript -e "BiocManager::install(version = '3.17')"

RUN yum -y install curl-devel
RUN yum -y install libxml2-devel
RUN yum -y install openssl-devel

ENV R_LIBS_USER /usr/lib64/R/library

RUN Rscript -e "BiocManager::install('beadarray', dependencies = TRUE)"
RUN Rscript -e "BiocManager::install('preprocessCore', dependencies = TRUE, configure.args='--disable-threading', force = TRUE)"


#Get peer:
RUN yum install -y cmake
RUN yum install -y swig

RUN git clone https://github.com/PMBio/peer.git /peer
WORKDIR peer
RUN mkdir build && cd build && \
    cmake -D BUILD_PEERTOOL=1 ./.. && \
    make && \
    make install
WORKDIR /

#Get plink
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
RUN unzip plink_linux_x86_64_20230116.zip
RUN mv plink /bin/

#Get Eigensoft:
RUN git clone https://github.com/DReichLab/EIG.git
WORKDIR EIG/src
RUN make
RUN make install
WORKDIR /
ENV PATH="/EIG/bin:${PATH}"


RUN Rscript -e "install.packages('qqman')"
RUN yum -y install pigz



