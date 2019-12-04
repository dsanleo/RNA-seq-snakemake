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


echo Start time : `date`
snakemake -p \
	--snakefile Snakefile \
	--cluster-config cluster_config.yml \
	--cluster "qsub -N rnaseqAnalysis -pe BWA 10  -q all.q -cwd -V" \
	--jobs 100 \
	--latency-wait 300 \
	--max-jobs-per-second 10 \


echo End time : `date`
