#!/bin/bash
cd "$(dirname "$0")/../.."
PROJECT_DIR="$(realpath .)"

if [ ! -d ${PROJECT_DIR}/venv ]; then
    echo "Installing env in ${PROJECT_DIR}/venv"
    python3.10 -m venv venv
fi

source venv/bin/activate
pip install -r ${PROJECT_DIR}/scripts/setup/requirements.txt
