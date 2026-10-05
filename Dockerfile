# OmicOne — dependency-free distribution.
# Everything (R + Seurat + Bioconductor + scop + the app) is baked into this image, so a
# user only needs Docker. Build once, then anyone can run:
#
#   docker run --rm -p 3838:3838 -m 16g <image>
#   # then open http://localhost:3838
#
# The Bioconductor base image ships R plus a prebuilt Bioconductor toolchain,
# which makes the Bioc dependencies (scDblFinder, SingleR, scater/scran) install
# quickly and reliably.

# Bioc 3.22 / R 4.5: the release the pipelines were checked against (2026-10).
FROM bioconductor/bioconductor_docker:RELEASE_3_22

LABEL org.opencontainers.image.title="OmicOne" \
      org.opencontainers.image.description="Local interactive multi-omics analysis app (single-cell and WES complete; bulk/spatial/integration planned)" \
      org.opencontainers.image.source="https://github.com/Tianqi-Ma/OmicOne"

# System libs occasionally needed by leiden/igraph/plotly stacks are already in
# the Bioconductor base. Install R package dependencies in a cached layer.
RUN R -e "install.packages(c( \
      'shiny','bslib','ggplot2','Matrix','plotly','DT','shinyWidgets', \
      'promises','future','progressr','remotes','Seurat','SeuratObject', \
      'harmony','survival','patchwork','leidenbase','igraph','RANN','maxstat','pROC', \
      'hdf5r','data.table'), \
      repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c( \
      'SingleCellExperiment','SummarizedExperiment','scater','scran', \
      'scDblFinder','glmGamPoi','SingleR','celldex','UCell','clusterProfiler', \
      'ComplexHeatmap','slingshot','maftools','MAST','edgeR','limma','DESeq2'), update=FALSE, ask=FALSE)"

# Mutational-signature extraction needs a reference genome and NMF. This layer
# is large (~700 MB for the BSgenome); drop it if you never run that step.
RUN R -e "install.packages(c('NMF','mclust'), repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg19', update=FALSE, ask=FALSE)"

# The analysis + plotting engine (scop) and its ecosystem. scop pulls a large
# tree; give it its own layer. LIANA/mascarade/copykat for cell-cell comm & CNV.
# scop is pinned to the commit the wrappers in fct_scop.R were checked against
# (0.9.2, 2026-09-12); bump it together with tests/testthat/test-scop.R.
RUN R -e "remotes::install_github('mengxu98/scop@bea13be5', upgrade='never')" && \
    R -e "remotes::install_github(c('saezlab/liana','alserglab/mascarade','navinlabcode/copykat', \
      'immunogenomics/presto','jinworks/CellChat'), upgrade='never')"

# Pre-bake the Python/conda environment for the Python-backed analyses
# (scVelo, PAGA, Palantir, scVI, scanorama, BBKNN). Doing this at BUILD time
# means users never wait for PrepareEnv() at runtime.
RUN R -e "reticulate::install_miniconda()" && \
    R -e "tryCatch(scop::PrepareEnv(), error=function(e) message('PrepareEnv at build: ', conditionMessage(e)))"

# Install the app itself (copy source and install from local path).
WORKDIR /opt/OmicOne
COPY . /opt/OmicOne
RUN R -e "remotes::install_local('/opt/OmicOne', dependencies = FALSE, upgrade = 'never')"

EXPOSE 3838

# Bind to 0.0.0.0 so the host browser can reach the container; do NOT auto-open a
# browser inside the container. Users open http://localhost:3838 themselves.
CMD ["R", "-e", "OmicOne::run_app(host='0.0.0.0', port=3838, launch.browser=FALSE)"]
