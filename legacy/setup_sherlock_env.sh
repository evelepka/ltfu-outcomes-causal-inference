#!/bin/bash
set -e

echo "Starting Sherlock GIS Environment Setup (CentOS 7 Compatible)..."

rm -rf ~/miniconda3
mkdir -p ~/miniconda3
cd ~/miniconda3

echo "Downloading Minconda Py39 (GLIBC 2.17 compatible)..."
wget -q https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh -O miniconda.sh
echo "Installing Miniconda silently..."
bash miniconda.sh -b -f -p ~/miniconda3
rm -rf miniconda.sh

# Initialize conda in this shell session
source ~/miniconda3/bin/activate

echo "Creating tb_geo_env with conda-forge python=3.10 geopandas..."
conda create -y -n tb_geo_env -c conda-forge python=3.10 geopandas earthengine-api

echo "Activating tb_geo_env..."
conda activate tb_geo_env

echo "Installing geobr via pip..."
pip install geobr

echo "Environment setup complete!"

echo "Now running the data acquisition scripts using the new environment..."
cd ~/tb_geocoding/
python3 02b_download_historical_ibge.py
python3 02c_extract_mapbiomas.py

echo "All tasks finished entirely automatically!"
