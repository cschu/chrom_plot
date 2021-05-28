import os
import csv

def read_genemap(geneset_file):
    gene_map = dict()
    for line in open(geneset_file):
        if line.startswith("#"):
            timepoints = set(map(int, line.replace("#", "").replace("dpi", "").strip().split(" ")[:-1]))
            col = "#cc0000"
            if 7 in timepoints:
                col = "#0000cc"
                if len(timepoints) > 1:
                    col = "#00cc00"
        else:
            gene_map[line.strip().strip("*")] = col
    return gene_map

def read_expression_map(fn, threshold=0.25):
    expression_map = dict()
    for gene, *values, _, _ in csv.reader(open(fn), delimiter="\t"):
        values = list(map(float, values))
        early = sum(values[:-3]) / 9
        late = sum(values[-3:]) / 3
        if early > late:
            expression_map[gene] = "#cc0000" if late / early < threshold else "#00cc00"
        elif late > early:
            expression_map[gene] = "#0000cc" if early / late < threshold else "#00cc00"
        else:
            expression_map[gene] = "#00cc00"

    return expression_map

def read_genes(gene_file):
    genes = dict()                                                                                        
    for row in csv.reader(open(gene_file), delimiter="\t"):                                                   
        genes.setdefault(row[0], list()).append((row[8].split(";")[0].split("=")[1], int(row[3])))    
    return genes

def read_regions(region_file):
    regions = dict()

    y_last = None
    for i, (region, chrom, color) in enumerate(csv.reader(open(region_file), delimiter="\t")):
        region, coords = region.split("=")[1].split(":")
        chrom = int(chrom.split("=")[1])
    
        length = int(coords.split("..")[1])
        if False: #chrom < 11: #length >= 3000000:
            y_offset = i + 1
            x_offset = 0
            last = y_offset
        else:
            y_offset = i + 1 #last
            x_offset = (1 if chrom % 2 == 1 else 0)#i % 2
            #x_offset = #i % 2
            #last += (1 if chrom % 2 else 0)
            
        regions[region] = {
            "chrom": chrom,
            "length": length,
            #"y_offset": (i + 1) if length >= 3000000 else i,
            "y_offset": y_offset,
            #i+1,
            #i // 2 + 1,
            #"x_offset": 0 if length >= 3000000 else i % 2
            "x_offset": x_offset,
            #0,
            #i % 2
            "chr_color": color
        }

    return regions

def read_proteins(protein_files):
    print("PF", protein_files)
    proteins = dict()
    for protein_file in protein_files:
        if protein_file:
            proteins[os.path.basename(protein_file)] = read_genes(protein_file)
    return proteins



