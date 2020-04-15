#!/bin/bash

################################################################################
# HCV genotyping pipeline preparation script.
################################################################################
# Objective: this script prepares the system installing all required software
# 	and dependencies.
#
# IMPORTANT: run this script as root user (sudo). Tested on Ubuntu 17.10.
#
# HCVgenotool
# Copyright (C) 2018  David Piñeyro
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
# Author: David Piñeyro <igtp.genbioinfo@gmail.com>
# License: GPL-3.0-or-later
# Version: 0.0.1
# Date: February 21st, 2018
################################################################################

################################################################################
# Setup hcv_genotyping virtualenv.
################################################################################
# Update ubuntu packages.
apt-get update
# Install virtualenv.
apt-get install virtualenv
# Install python3 pip.
apt-get install python3-pip
# Install python2 pip.
#apt-get install python-pip
# Turn on bash autocomplete for pip (optional, but recommended).
pip3 completion --bash >> ~/.bashrc
#pip completion --bash >> ~/.bashrc
source ~/.bashrc
# Install virtualenvwrapper locally via pip3.
pip3 install --user virtualenvwrapper
#pip install --user virtualenvwrapper
echo "export VIRTUALENVWRAPPER_PYTHON=/usr/bin/python3" >> ~/.bashrc
#echo "export VIRTUALENVWRAPPER_PYTHON=/usr/bin/python" >> ~/.bashrc
echo "source ~/.local/bin/virtualenvwrapper.sh" >> ~/.bashrc
# Export the WORKON_HOME variable which contains the directory in which our
# virtual environments are to be stored (~/.virtualenvs).
export WORKON_HOME=~/.virtualenvs
# Create ~/.virtualenvs
mkdir $WORKON_HOME
# Put this variable in ~/.bashrc
echo "export WORKON_HOME=$WORKON_HOME" >> ~/.bashrc
# Makes sure that if pip creates an extra virtual environment, it is also placed
# in our WORKON_HOME directory.
echo "export PIP_VIRTUALENV_BASE=$WORKON_HOME" >> ~/.bashrc
# Load the changes.
source ~/.bashrc
# Create our virtualenv (and activate it).
mkvirtualenv -p python3 hcv_genotyping

################################################################################
# Installing Cutadapt.
################################################################################
pip3 install cutadapt

################################################################################
# Installing python3 required modules.
################################################################################
pip3 install numpy
pip3 install biopython
pip3 install pandas

################################################################################
# Installing ea-utils/fastq-join.
################################################################################
# Downloaded from:
#https://github.com/ExpressionAnalysis/ea-utils/zipball/master
# to compile:
#make

#Requeriments on UBUNTU :
#apt-get install libgsl0-dev
# Change the require lines in the test files t/multx.t, t/join.t, t/mcf.t from
#require (dirname(__FILE__) . "/test-prep.pl");
#to
#require "/path/to/test-prep.pl"

#Requeriments on CENTOS/REDHAT :
#rpm -i gsl-devel
