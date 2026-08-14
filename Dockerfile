# bac2food — reproducible environment, Docker variant.
#
# The SUPPORTED container is the Singularity/Apptainer image built from bac2food.def:
# it needs no daemon and no root, so it runs on HPC systems where Docker is unavailable.
# This Dockerfile is kept for anyone who does have Docker.
#
# The image carries the CODE and its dependencies. It deliberately does NOT carry the
# reference tables: they are 2.1 GB, not every source's licence permits redistribution
# inside an image, and they are versioned separately on Zenodo. Mount them at run time:
#
#   docker build -t bac2food .
#   docker run --rm \
#       -v /path/to/zenodo/deposit:/data/bac2food/exports:ro \
#       -v "$PWD/out":/out \
#       bac2food \
#       python 4_predict/bac2food_predict.py \
#           --mag_tsv /out/my_metagenome.tsv --out_prefix /out/run1 --jobs 4
#
# Python is pinned to 3.10.8 because that is the interpreter the published results were
# produced with, and pyarrow 11 (see requirements.txt) is built against it.
FROM python:3.10.8-slim-bookworm

# procps is for the memory reporting the predictor does at startup; nothing else is
# needed at run time. Cleaned in the same layer so the list is not baked into the image.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/bac2food

# Dependencies first, in their own layer: requirements.txt changes far less often than
# the code, so a code edit does not re-run the install.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# The predictor peaks around 5.7 GB on a full metagenome, driven by loading the
# reference layer rather than by the community size. Give the container at least 8 GB
# (docker run --memory=8g); below that the Linux OOM killer takes it, which surfaces as
# an unexplained exit 137 rather than a Python error.
ENV PYTHONUNBUFFERED=1 \
    BAC2FOOD_CACHE=/data/bac2food/cache

# Derived caches (the 1.1 GB modeled index, the parquet reference cache) are rebuilt on
# first run and belong on a mounted volume, not in the image layers.
VOLUME ["/data/bac2food"]

CMD ["python", "4_predict/bac2food_predict.py", "--help"]
