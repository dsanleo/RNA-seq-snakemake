#!/bin/bash
# **************************************************************************
# * (C)opyright <+$YEAR$;R+> by David San Leon Granado++++++++++++++++++++>*
# *                                                                        *
# * This program is free software; you can redistribute it and/or modify   *
# *  it under the terms of the GNU General Public License as published by  *
# *  the Free Software Foundation version 2 of the License.                *
# *                                                                        *
# * This program is distributed in the hope that it will be useful,        *
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
# *  GNU General Public License for more details.                          *
# *                                                                        *
# * You should have received a copy of the GNU General Public License      *
# *  along with this program; if not, write to the                         *
# *  Free Software Foundation, Inc.,                                       *
# *  59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
# *                                                                        *
# * Contact author at:                                                     *
# *  David San Leon Granado++++++++++++++++++++++++++++++++++++++++++++++> *
# *  d.san_leon_granado@lumc.nl++++++++++++++++++++++++++++++++++++++++++> *
# **************************************************************************


# Date:
date
# Version: 1.0.0
# Usage: sh 

#path to drmaa in Shark
export DRMAA_LIBRARY_PATH=/share/isilon/system/local/OpenGridScheduler/gridengine/lib/linux-x64/libdrmaa.so

userDir=`whoami`

echo Start time : `date`

if ! conda_loc="$(type -p conda)" ; then
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -p /exports/humgen/$userDir/miniconda3
	echo "PATH=\$PATH:/exports/"$userDir"/miniconda/bin" >> $HOME/.bashrc
	export PATH=$PATH:/exports/humgen/$userDir/miniconda3/bin
	conda config --add channels defaults
	conda config --add channels conda-forge
	conda config --add channels bioconda
	conda install -y snakemake
	conda install -y drmaa
fi


snakemake -p \
	--snakefile Snakefile \
	--latency-wait 90 \
	--use-conda \
	--cluster-config $(pwd)/cluster_config.yml \
	--drmaa " -N rnaseqAnalysis -pe BWA {threads} -l h_vmem=20G -q all.q -cwd -V" \
	--drmaa-log-dir $(pwd)/cluster_logs \
	--jobs 100 \
	--wait-for-files \
	--restart-times 2 \
	--rerun-incomplete


echo End time : `date`
