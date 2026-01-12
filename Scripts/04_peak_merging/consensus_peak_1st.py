import os
import pandas as pd
import pyranges as pr
from pycisTopic.iterative_peak_calling import get_consensus_peaks

# Chromosome sizes
chromsizes = pd.read_table(
    "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
chromsizes.head()

# Narrow peaks
def load_narrow_peak(path):
    """Load MACS2 narrow peak files as :class:`pr.PyRanges`."""
    narrow_peak = pd.read_csv(path, sep="\t", header=None)
    narrow_peak.columns = [
        "Chromosome", 
        "Start", 
        "End", 
        "Name", 
        "Score", 
        "Strand", 
        "FC_summit", 
        "-log10_pval", 
        "-log10_qval", 
        "Summit",]
    narrow_peak_pr = pr.PyRanges(narrow_peak)
    return narrow_peak_pr

narrow_peak_dict = {}
with open(os.environ['PEAKPATH'], "r") as f:
    for line in f:
        k, v = line.strip().split(": ")
        narrow_peak_dict.update({k: load_narrow_peak(v)})

# Other param
peak_half_width=250
path_to_blacklist="/scratch/user/s4799891/pycisTopic/blacklist/hg38-blacklist.v2.bed"

# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist)

consensus_peaks.to_bed(
    path = os.environ['OUTDIR'] + "/consensus_regions.bed",
    keep =True,
    compression = 'infer',
    chain = False)
