import os
import sys
import csv
import argparse
import math

from collections import Counter

"""
<svg version="1.1"
     baseProfile="full" 
     width="1000" height="500"
     xmlns="http://www.w3.org/2000/svg">

  <!-- <rect width="100%" height="100%" fill="red" /> -->

  <!-- <circle cx="150" cy="100" r="80" fill="green" /> -->

  <!-- <text x="150" y="125" font-size="60" text-anchor="middle" fill="white">SVG</text> -->

  <line x1="100" y1="100" x2="200" y2="100" fill="black" stroke="black" />

</svg>


ID=NC_030986.1:1..6854980    chromosome=1
ID=NC_030987.1:1..5577357    chromosome=2
"""

def write_svg(body, width=1000, height=500, stream=sys.stdout):
    print('<svg version="1.1" '
          'baseProfile="full" '
          'width="{width}" height="{height}" '
          'xmlns="http://www.w3.org/2000/svg">'.format(width=width, height=height), file=stream, flush=True)
    print('<defs>\n'
          '<!-- arrowhead marker definition -->\n'
          '<marker id="arrow" viewBox="0 0 10 10" refX="5" refY="5"'
          ' markerWidth="6" markerHeight="6"'
          ' orient="auto-start-reverse">'
          '<path d="M 0 0 L 10 5 L 0 10 z" />'
          '</marker></defs>', file=stream, flush=True)
    print(body, file=stream, flush=True)
    print('</svg>', file=stream, flush=True)

def write_legend_rect(**args):
    return '<rect x="{x}" y="{y}" width="{width}" height="{height}" stroke="#000000" fill="#ffffff" stroke-width="1" />'.format(**args)
def write_rect(**args):
    return '<rect x="{x}" y="{y}" width="{width}" height="{height}" stroke="#000000" fill="{fill}" stroke-width="{stroke_width}" rx="{rx}" />'.format(**args)
def write_outer_rect(**args):
    return '<rect x="{x}" y="{y}" width="{width}" height="{height}" stroke="#000000" fill="none" stroke-width="1.5" rx="3" />'.format(**args)

def write_line(**args):
    return '<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{stroke_width}"/>\n'.format(**args)
def write_dashed_line(**args):
    return '<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{stroke_width}" stroke-dasharray="2,2" />\n'.format(**args)

def write_arrow(**args):
    #return '<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{stroke_width}" marker-start="url(#arrow)" marker-end="url(#arrow)" />'.format(**args)
    return '<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{stroke_width}" marker-start="url(#arrow)" />'.format(**args)
    

def write_text(**args):
    #return '<text x="{x}" y="{y}" font-size="10" font-family="sans-serif" text-anchor="middle">{text}</text>\n'.format(**args)
    return '<text x="{x}" y="{y}" font-size="{font_size}" font-family="sans-serif" text-anchor="middle" dominant-baseline="central">{text}</text>\n'.format(**args)
def write_text_right(**args):
    return '<text x="{x}" y="{y}" font-size="8" font-family="sans-serif" text-anchor="end">{text}</text>\n'.format(**args)
def write_text_left(**args):
    return '<text x="{x}" y="{y}" font-size="8" font-family="sans-serif" dy="0.25em">{text}</text>\n'.format(**args)

def write_circle(**args):
    return '<circle cx="{x}" cy="{y}" r="{radius}" fill="{fill}"/>'.format(**args)


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
        if chrom < 11: #length >= 3000000:
            y_offset = i + 1
            x_offset = 0
            last = y_offset
        else:
            y_offset = last
            x_offset = (1 if chrom % 2 == 1 else 0)#i % 2
            #x_offset = #i % 2
            last += (1 if chrom % 2 else 0)
            
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("region_file", type=str, help="A tab-separated file with region data to be visualised")
    #ap.add_argument("gene_coords", type=str, help="A gff file with genes to be visualised on the regions")
    #ap.add_argument("gene_groups", type=str, help="A file mapping genes to groups")
    #ap.add_argument("expression_map", type=str, help="A file with expression data for each gene")
    ap.add_argument("--proteins", "-p", type=str, help="A list of protein coordinates", default="")
    args = ap.parse_args()

    regions = read_regions(args.region_file)
    #genes = read_genes(args.gene_coords)
    #gene_map = read_genemap(args.gene_groups)
    #expression_map = read_expression_map(args.expression_map)

    proteins = dict()
    for protein_file in args.proteins.split(","):
        if protein_file:
            proteins[os.path.basename(protein_file)] = read_genes(protein_file)

    #print("GENE_MAP", len(gene_map))
    #print("EXPRESSION_MAP", len(expression_map))
    #print("GENES", sum(len(v) for v in genes.values()))
    body = ""
    y_start = 550
    x_offsets = 50, 525
    max_axis_length = 425
    
    with open("chrom_plot_ls.proteins.svg", "wt") as svg_out:
        max_region_length = max(region["length"] for region in regions.values())
        scale_x_length = 1e6 / max_region_length * max_axis_length
        body += write_line(x1=50, y1=50, x2=50 + scale_x_length, y2=50, color="#000000", stroke_width=1.5)
        body += write_line(x1=50, y1=45, x2=50, y2=55, color="#000000", stroke_width=1.5)
        body += write_line(x1=50 + scale_x_length, y1=45, x2=50 + scale_x_length, y2=55, color="#000000", stroke_width=1.5)
        body += write_text(x=50 + (scale_x_length)/2, y=65, text="1 Mbp", font_size=10)
    
        for i, (region_id, region) in enumerate(regions.items()):
            y = y_start * region["y_offset"]
            
            region_length = int(region["length"] / max_region_length * max_axis_length)
            if region["x_offset"]:
                x = 300 #last_region_end + 100
            else:
                last_region_end = region_length + 50
                x = 50
            #x = x_offsets[region["x_offset"]]
    
            #body += write_text(x=x+(max_axis_length//2), y=y-30, text="Chr{chrom} ({length_mb:.2f} Mbp)".format(**region, length_mb=region["length"]/1e6))
            body += '<g id="{gid}">\n'.format(gid=region["chrom"])
            #body += write_text(x=x, y=y-30, text="Chr{chrom} ({length_mb:.2f} Mbp)".format(**region, length_mb=region["length"]/1e6))
    
            x = 35 + x_offsets[0] * (i + 1)
            #\n({length_mb:.2f} Mbp)".format(**region, length_mb=region["length"]/1e6))
            body += write_rect(x=x, y=y_start-region_length, width=10, height=region_length, fill=region["chr_color"], stroke_width=1.0, rx=3)
            body += write_text(x=x + 5, y=y_start + 15, text="Chr{chrom}".format(**region), font_size=10)
    
            # deal with accessory regions on Chr1, Chr2
            if i < 2:
                acc_height_factor = 7/8 if i == 0 else 6/8
                acc_height = (1e6 * acc_height_factor) / max_region_length * max_axis_length
                acc_y = y_start-region_length + acc_height
                body += write_rect(x=x, y=y_start-region_length, width=10, height=acc_height, fill="#bbbbbb", stroke_width=0.0, rx=3)
                body += write_dashed_line(x1=x, x2=x+10, y1=acc_y, y2=acc_y, color="#000000", stroke_width="1.0")
    
            #for gene, pos in genes.get(region_id, list()):
            #    #x_pos = x + int(gene / region["length"] * region_length)
            #    y_pos = y_start - int(pos / region["length"] * region_length)
            #    col = gene_map.get(gene, None)
            #    if not col:
            #        print("xxx", gene)
            #    else:
            #        if col == "#00cc00":
            #            col = expression_map.get(gene, col)
            #            print(gene, col)
            #        if True: #col == "#0000cc":
            #            body += write_line(x1=x, y1=y_pos, x2=x+10, y2=y_pos, color=col, stroke_width=1.5)

            bins = Counter()
            binsize = 500000
            for pid, (protein_set, protein_coords) in enumerate(proteins.items()):
                for protein, pos in protein_coords.get(region_id, list()):
                    bins[pos // binsize] += 1
                    #y_pos = y_start - int(pos / region["length"] * region_length)
                    #col = "#00cccc" if pid == 0 else "#cc0000"
                    #col = "#003333"
                    #body += write_dashed_line(x1=x+11, y1=y_pos, x2=x+21, y2=y_pos, color=col, stroke_width=1)
            print(bins)
            for pos, count in bins.items():
                y_pos = y_start - int((pos * binsize) / region["length"] * region_length)
                col = "#003333"
                col = "#000"
                # body += write_dashed_line(x1=x+13, y1=y_pos, x2=x+15, y2=y_pos, color=col, stroke_width=1)
                body += write_arrow(x1=x+13, y1=y_pos, x2=x+15, y2=y_pos, color=col, stroke_width=1)
                radius = math.sqrt(4 * count/math.pi)
        
                #body += write_circle(x=x+21, y=y_pos, radius=radius, fill=col) 
                body += write_text(x=x+20, y=y_pos, text=count, font_size=8)
                # return '<text x="{x}" y="{y}" font-size="10" font-family="sans-serif" text-anchor="middle">{text}</text>\n'.format(**args)

            # overpaint chromosome borders
            body += write_rect(x=x, y=y_start-region_length, width=10, height=region_length, fill="none", stroke_width=1.5, rx=3)
            
            # 100kbp yellow marker-band for validating coordinate order 
            #body += write_line(x1=x, y1=y_start - int(100000 / region["length"] * region_length), x2=x+10, y2=y_start - int(100000 / region["length"] * region_length), color="#ffff00", stroke_width=3)
    
            #body += write_line(x1=x, y1=y_start, x2=x, y2=y_start - region_length, color="#000000", stroke_width=3.0) 
            #body += write_line(x1=x+10, y1=y_start, x2=x+10, y2=y_start - region_length, color="#000000", stroke_width=1.0) 
            #return '<rect x="{x}" y="{y}" width="{width}" height="{height}" stroke="#000000" fill="#111111" rx="15"/>'.format(**args)
            body += "</g>\n"
            
    
            #x = x_1 if i % 2 == 0 else x_2
        
        #print(regions)    
        #write_svg('<line x1="100" y1="100" x2="200" y2="100" fill="black" stroke="black" />')
        #write_svg(body, width=600, height=1500, stream=svg_out)
        #body += write_line(x1=75, y1=y_start, x2=75, y2=y_start - max_axis_length, stroke_width="1.0", color="#000000")
        body += write_text_right(x=73, y=(y_start - max_axis_length - 3), text="{:.2f} Mbp".format(max_region_length/1e6))
        body += write_text_right(x=73, y=(y_start - max_axis_length//2 - 3), text="{:.2f} Mbp".format((max_region_length/2)/1e6))
        #body += write_text_right(x=73, y=y_start, text="1 bp")
        #return '<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{stroke_width}"/>\n'.format(**args)
        #return '<text x="{x}" y="{y}" font-size="10" font-family="sans-serif" text-anchor="middle">{text}</text>\n'.format(**args)
    
        legend_x = 915//2 - 150
        legend_width = 200
    
        body += write_legend_rect(x=legend_x, y=y_start+25, width=legend_width, height=40)
    
        #body += write_line(x1=legend_x + 10, y1=y_start+25+12, x2=legend_x + 20, y2=y_start+25+12, color="#cc0000", stroke_width=1.5)
        #body += write_text_left(x=legend_x + 20 + 2, y=y_start+25+12, text="early (1,2,3 dpi)")
    
        #body += write_line(x1=legend_x + 20 + 2 + 100, y1=y_start+25+12, x2=legend_x + 20 + 2 + 110, y2=y_start+25+12, color="#0000cc", stroke_width=1.5)
        #body += write_text_left(x=legend_x + 20 + 2 + 110 + 2, y=y_start+25+12, text="late (7 dpi)")
    
        #body += write_line(x1=legend_x + 20 + 2 + 200, y1=y_start+25+12, x2=legend_x + 20 + 2 + 210, y2=y_start+25+12, color="#00cc00", stroke_width=1.5)
        #body += write_text_left(x=legend_x + 20 + 2 + 210 + 2, y=y_start+25+12, text="all")
    
    
        legend_y = y_start+25+12#+20
    
        body += write_line(x1=legend_x+ 10, y1=legend_y, x2=legend_x + 20, y2=legend_y, color="#eeeeee", stroke_width=1.5)
        body += write_text_left(x=legend_x + 20 + 2, y=legend_y, text="core")
        
        body += write_line(x1=legend_x + 20 + 2 + 100, y1=legend_y, x2=legend_x + 20 + 2 + 110, y2=legend_y, color="#bbbbbb", stroke_width=1.5)
        # body += write_text_left(x=legend_x + 20 + 2 + 110 + 2, y=legend_y, text="core (fast-evolving)")
        body += write_text_left(x=legend_x + 20 + 2 + 110 + 2, y=legend_y, text="accessory")

        body += write_arrow(x1=legend_x+12, y1=legend_y+15, x2=legend_x+12, y2=legend_y+15, color=col, stroke_width=1)
        #body += write_text(x=legend_x+15, y=legend_y+15, text="n", font_size=8)
        body += write_text_left(x=legend_x + 20, y=legend_y+15, text="n proteins encoded per 500kbp window")
        #for np in range(1, 7): #min(bins.values()), max(bins.values()) + 1):
        #    body += write_circle(x=legend_x + 10 + 15 * np, y=legend_y+30, radius=math.sqrt(4 * np/math.pi), fill="#000")
        #    body += write_text(x=legend_x + 10 + 15 * np + 8, y=legend_y+30, font_size=6, text=np)
    
        #body += write_line(x1=legend_x + 20 + 2 + 200, y1=legend_y, x2=legend_x + 20 + 2 + 210, y2=legend_y, color="#888888", stroke_width=1.5)    
        #body += write_text_left(x=legend_x + 20 + 2 + 210 + 2, y=legend_y, text="accessory")
        
         #return '<rect x="{x}" y="{y}" width="{width}" height="{height}" stroke="#000000" fill="#cccccc" stroke-width="1" rx="3" />'.format(**args)
        write_svg(body, width=1000, height=650, stream=svg_out)


if __name__ == "__main__":
    main()

"""
NC_030986.1    RefSeq    gene    298241    299509    .    -    .    ID=gene-FOXG_10922;Dbxref=GeneID:28952358;Name=FOXG_10922;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10922
NC_030986.1    RefSeq    gene    298701    301474    .    +    .    ID=gene-FOXG_10923;Dbxref=GeneID:28952359;Name=FOXG_10923;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10923
NC_030986.1    RefSeq    gene    345744    347414    .    +    .    ID=gene-FOXG_10940;Dbxref=GeneID:28952374;Name=FOXG_10940;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10940
NC_030986.1    RefSeq    gene    371296    372054    .    -    .    ID=gene-FOXG_10949;Dbxref=GeneID:28952383;Name=FOXG_10949;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10949
NC_030986.1    RefSeq    gene    373454    374362    .    +    .    ID=gene-FOXG_10950;Dbxref=GeneID:28952384;Name=FOXG_10950;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10950
NC_030986.1    RefSeq    gene    437332    438668    .    +    .    ID=gene-FOXG_10974;Dbxref=GeneID:28952404;Name=FOXG_10974;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10974
NC_030986.1    RefSeq    gene    595549    597284    .    -    .    ID=gene-FOXG_11033;Dbxref=GeneID:28952461;Name=FOXG_11033;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_11033
NC_030986.1    RefSeq    gene    835811    838567    .    -    .    ID=gene-FOXG_11103;Dbxref=GeneID:28952530;Name=FOXG_11103;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_11103
NC_030986.1    RefSeq    gene    890007    892816    .    +    .    ID=gene-FOXG_11116;Dbxref=GeneID:28952541;Name=FOXG_11116;end_range=892816,.;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_11116;partial=true
NC_030986.1    RefSeq    gene    1926815    1930446    .    -    .    ID=gene-FOXG_00102;Dbxref=GeneID:28942425;Name=FOXG_00102;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_00102
"""
