FROM rocker/verse
RUN R -e "install.packages('DynTxRegime')"
RUN apt-get update -y && apt-get install -y ca-certificates curl apt-transport-https lsb-release gnupg
COPY microsoft.asc /etc/apt/trusted.gpg.d/microsoft.gpg 
RUN apt-get update -y && apt-get install -y azure-cli
RUN apt-get update -y && apt-get install -y jq emacs
RUN R -e "devtools::install_github('Azure/rAzureBatch')"
RUN R -e "devtools::install_github('Azure/doAzureParallel')"
RUN R -e "devtools::install_github('pROC')"
