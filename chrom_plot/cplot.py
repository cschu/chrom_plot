import os
import sys
import csv
import argparse
import math

from collections import Counter

from chrom_plot import *
import svg
from data_io import read_genemap, read_expression_map, read_genes, read_regions, read_proteins


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("region_file", type=str, help="A tab-separated file with region data to be visualised")
    ap.add_argument("--gene_coords", type=str, help="A gff file with genes to be visualised on the regions", default="")
    ap.add_argument("--gene_groups", type=str, help="A file mapping genes to groups", default="")
    ap.add_argument("--expression_map", type=str, help="A file with expression data for each gene", default="")
    ap.add_argument("--proteins", "-p", type=str, help="A list of protein coordinates", default="")
    ap.add_argument("--orientation", choices=("landscape", "portrait", "l", "p"), default="l")
    ap.add_argument("--outfile", "-o", type=str, default="chrom_plot.tmp.svg")
    args = ap.parse_args()

    plot_config = {
        "genes": read_genes(args.gene_coords) if args.gene_coords else None,
        "gene_groups": read_genemap(args.gene_groups) if args.gene_groups else None,
        "expression_map": read_expression_map(args.expression_map) if args.expression_map else None,
        "proteins": read_proteins(args.proteins.split(","))
    }

    print(plot_config)

    regions = read_regions(args.region_file)
    print(regions)
    cplot = (ChromPlotLandscape if args.orientation.startswith("l") else ChromPlotPortrait)(regions, **plot_config)
    with open(args.outfile, "w") as svg_out:
        cplot.draw(out=svg_out)
    
    

if __name__ == "__main__":
    main()
