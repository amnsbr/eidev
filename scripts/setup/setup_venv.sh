#!/bin/bash
cd "$(dirname "$0")/../.."
PROJECT_DIR="$(realpath .)"

if [ ! -d ${PROJECT_DIR}/venv ]; then
    echo "Installing env in ${PROJECT_DIR}/venv"
    python3.10 -m venv venv
fi

source venv/bin/activate
pip install -r ${PROJECT_DIR}/scripts/setup/requirements.txt

# install cubnm
export CUBNM_NOISE_WHOLE=1
pip install git+https://github.com/amnsbr/cuBNM@be473d471ff4aa8bf92b8971b075a2d265adb1b5 -vvv