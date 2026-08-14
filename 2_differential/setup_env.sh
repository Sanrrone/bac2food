#!/bin/bash

echo "Creating Python virtual environment..."
python3 -m venv ec2food_env

echo "Activating environment..."
# Note: If you are on Windows using Git Bash, this path might be ec2food_env/Scripts/activate
source ec2food_env/bin/activate

echo "Upgrading pip..."
pip install --upgrade pip

echo "Installing required packages for the Differential Prescription Engine..."
# Pinning versions to ensure stability
pip install pandas==2.2.1 numpy==1.26.4 pyarrow==15.0.0

echo "=========================================================="
echo "Environment setup complete!"
echo "To activate this environment in the future, run:"
echo "source ec2food_env/bin/activate"
echo "=========================================================="
