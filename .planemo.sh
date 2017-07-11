#!/bin/bash

pip install planemo
planemo database_create galaxy
planemo test --install_galaxy --galaxy_branch release_17.01 --skip_venv --no_conda_auto_install --no_conda_auto_init --postgres galaxy/wrapper