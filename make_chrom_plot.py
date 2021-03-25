import sys
import csv


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


ID=NC_030986.1:1..6854980	chromosome=1
ID=NC_030987.1:1..5577357	chromosome=2
"""

def write_svg(body, width=1000, height=500, stream=sys.stdout):
	print('<svg version="1.1" '
		  'baseProfile="full" '
		  'width="{width}" height="{height}" '
          'xmlns="http://www.w3.org/2000/svg">'.format(width=width, height=height), file=stream, flush=True)
	print(body, file=stream, flush=True)
	print('</svg>', file=stream, flush=True)


def write_line(**args):
	return '<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="1.5"/>\n'.format(**args)
def write_text(**args):
	#return '<text x="{x}" y="{y}" font-size="10" font-family="sans-serif" text-anchor="middle">{text}</text>\n'.format(**args)
	return '<text x="{x}" y="{y}" font-size="10" font-family="sans-serif">{text}</text>\n'.format(**args)


def read_regions(region_file):
	regions = dict()


	y_last = None
	for i, (region, chrom) in enumerate(csv.reader(open(region_file), delimiter="\t")):
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
		}

	return regions

def read_genes(gene_file):
	genes = dict()
	for row in csv.reader(open(gene_file), delimiter="\t"):
		genes.setdefault(row[0], list()).append(int(row[3]))
	return genes


body = ""
regions = read_regions(sys.argv[1])
genes = read_genes(sys.argv[2])
y_start = 50
x_offsets = 50, 525
max_axis_length = 425

with open("chrom_plot.svg", "wt") as svg_out:
	max_region_length = max(region["length"] for region in regions.values())
	for region_id, region in regions.items():
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
		body += write_text(x=x, y=y-30, text="Chr{chrom} ({length_mb:.2f} Mbp)".format(**region, length_mb=region["length"]/1e6))
		body += "</g>\n"

		for gene in genes.get(region_id, list()):
			x_pos = x + int(gene / region["length"] * region_length)
			body += write_line(x1=x_pos, y1=y, x2=x_pos, y2=y-10, color="#006600")

		body += write_line(x1=x, y1=y, x2=x+region_length, y2=y, color="#000000") 
		

		#x = x_1 if i % 2 == 0 else x_2
	
	#print(regions)	
	#write_svg('<line x1="100" y1="100" x2="200" y2="100" fill="black" stroke="black" />')
	write_svg(body, width=600, height=1500, stream=svg_out)

"""
NC_030986.1	RefSeq	gene	298241	299509	.	-	.	ID=gene-FOXG_10922;Dbxref=GeneID:28952358;Name=FOXG_10922;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10922
NC_030986.1	RefSeq	gene	298701	301474	.	+	.	ID=gene-FOXG_10923;Dbxref=GeneID:28952359;Name=FOXG_10923;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10923
NC_030986.1	RefSeq	gene	345744	347414	.	+	.	ID=gene-FOXG_10940;Dbxref=GeneID:28952374;Name=FOXG_10940;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10940
NC_030986.1	RefSeq	gene	371296	372054	.	-	.	ID=gene-FOXG_10949;Dbxref=GeneID:28952383;Name=FOXG_10949;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10949
NC_030986.1	RefSeq	gene	373454	374362	.	+	.	ID=gene-FOXG_10950;Dbxref=GeneID:28952384;Name=FOXG_10950;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10950
NC_030986.1	RefSeq	gene	437332	438668	.	+	.	ID=gene-FOXG_10974;Dbxref=GeneID:28952404;Name=FOXG_10974;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_10974
NC_030986.1	RefSeq	gene	595549	597284	.	-	.	ID=gene-FOXG_11033;Dbxref=GeneID:28952461;Name=FOXG_11033;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_11033
NC_030986.1	RefSeq	gene	835811	838567	.	-	.	ID=gene-FOXG_11103;Dbxref=GeneID:28952530;Name=FOXG_11103;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_11103
NC_030986.1	RefSeq	gene	890007	892816	.	+	.	ID=gene-FOXG_11116;Dbxref=GeneID:28952541;Name=FOXG_11116;end_range=892816,.;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_11116;partial=true
NC_030986.1	RefSeq	gene	1926815	1930446	.	-	.	ID=gene-FOXG_00102;Dbxref=GeneID:28942425;Name=FOXG_00102;gbkey=Gene;gene_biotype=protein_coding;locus_tag=FOXG_00102
"""
