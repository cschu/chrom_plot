import sys
from collections import Counter

import svg

class ChromPlot:
    x_offsets = 50, 525
    max_axis_length = 425
    y_start = 50
    width, height = 1000, 650

    def __init__(self, regions, genes=None, gene_groups=None, expression_map=None, proteins=None):
        self.body = list()
        #return '<rect x="{x}" y="{y}" width="{width}" height="{height}" stroke="#000000" fill="{fill}" stroke-width="{stroke_width}" rx="{rx}" />'.format(**args)
        self.body.append(svg.write_rect(x=0, y=0, width=self.width, height=self.height, fill="#fff", stroke_width=0, rx=0))
        self.regions = regions
        self.genes = genes if genes else dict()
        self.gene_groups = gene_groups if gene_groups else dict()
        self.expression_map = expression_map if expression_map else dict()
        self.proteins = proteins if proteins else dict()

        self.max_region_length = max(region["length"] for region in regions.values())
        self.scale_x_length = 1e6 / self.max_region_length * ChromPlot.max_axis_length

    def calc_protein_bins(self, region_id, binsize=500000):
        bins = Counter()
        print("SELFP", self.proteins)
        for pid, (protein_set, protein_coords) in enumerate(self.proteins.items()):
            for protein, pos in protein_coords.get(region_id, list()):
                bins[pos // binsize] += 1
        return bins

    def draw(self, out=sys.stdout):
        self.draw_regions()
        self.draw_legend()
        #print(self.body)
        svg.write_svg("".join(self.body), width=self.width, height=self.height, stream=out)


class ChromPlotLandscape(ChromPlot):
    y_start = 550
    width = 1000
    height = 650

    def __init__(self, regions, **kwargs):
        super().__init__(regions, **kwargs)
        print("L")
    def calc_accessory_region(self, chrom, region_length, offset):
         acc_height_factor = 7/8 if chrom == 0 else 6/8
         acc_height = (1e6 * acc_height_factor) / self.max_region_length * ChromPlot.max_axis_length        
         acc_y = offset - region_length + acc_height
         return acc_height, acc_y

        
    def draw_regions(self):
        for i, (region_id, region) in enumerate(self.regions.items()):
            print(region_id, region)        
            region_length = int(region["length"] / self.max_region_length * ChromPlot.max_axis_length)
            x = 35 + ChromPlot.x_offsets[0] * (i + 1)
   
            self.body.extend([ 
                '<g id="{gid}">\n'.format(gid=region["chrom"]),
                svg.write_rect(x=x, y=self.y_start-region_length, width=10, height=region_length, fill=region["chr_color"], stroke_width=1.0, rx=3),
                svg.write_text(x=x + 5, y=self.y_start + 15, text="Chr{chrom}".format(**region), font_size=10)
            ])
    
            # deal with accessory regions on Chr1, Chr2
            if i < 2:
                acc_height, acc_y = self.calc_accessory_region(i, region_length, self.y_start)
                self.body.extend([
                    svg.write_rect(x=x, y=self.y_start-region_length, width=10, height=acc_height, fill="#bbbbbb", stroke_width=0.0, rx=3),
                    svg.write_dashed_line(x1=x, x2=x+10, y1=acc_y, y2=acc_y, color="#000000", stroke_width="1.0")
                ])

            for gene, pos in self.genes.get(region_id, list()):
                #x_pos = x + int(gene / region["length"] * region_length)
                y_pos = self.y_start - int(pos / region["length"] * region_length)
                col = self.gene_groups.get(gene, None)
                if not col:
                    print("xxx", gene)
                else:
                    if col == "#00cc00":
                        col = self.expression_map.get(gene, col)
                        print(gene, col)
                    self.body.append(svg.write_line(x1=x, y1=y_pos, x2=x+10, y2=y_pos, color=col, stroke_width=1.5))

            binsize = 500000
            bins = self.calc_protein_bins(region_id, binsize=binsize)
            for pos, count in bins.items():
                y_pos = self.y_start - int((pos * binsize) / region["length"] * region_length)
                print("FOSI")
                self.body.extend([
                    svg.write_arrow(x1=x+13, y1=y_pos, x2=x+15, y2=y_pos, color="#000", stroke_width=1),
                    svg.write_text(x=x+20, y=y_pos, text=count, font_size=10)
                ])

            # overpaint chromosome borders
            self.body.extend([
                svg.write_rect(x=x, y=self.y_start-region_length, width=10, height=region_length, fill="none", stroke_width=1.5, rx=3),
                "</g>\n"
            ])

        x = 35 + ChromPlot.x_offsets[0] * (len(self.regions) + 1)

        self.body.extend([
            svg.write_line(x1=x, y1=self.y_start, x2=x, y2=self.y_start - self.scale_x_length, color="#000000", stroke_width=1.5),
            svg.write_line(x1=x-5, y1=self.y_start, x2=x+5, y2=self.y_start, color="#000000", stroke_width=1.5),
            svg.write_line(x1=x-5, y1=self.y_start - self.scale_x_length, x2=x+5, y2=self.y_start - self.scale_x_length, color="#000000", stroke_width=1.5),
            svg.write_text_left(x=x+5, y=self.y_start - self.scale_x_length/2, text="1 Mbp", font_size=8),
        ])



    def draw_legend(self):
        legend_x, legend_width = 915//2 - 150, 200
        legend_y = self.y_start+25+12#+20

        if self.proteins:
            self.body.extend([
                svg.write_legend_rect(x=legend_x, y=self.y_start+25, width=legend_width, height=40),
                svg.write_line(x1=legend_x+ 10, y1=legend_y, x2=legend_x + 20, y2=legend_y, color="#eeeeee", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2, y=legend_y, text="core", font_size=8),
                svg.write_line(x1=legend_x + 20 + 2 + 100, y1=legend_y, x2=legend_x + 20 + 2 + 110, y2=legend_y, color="#bbbbbb", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2 + 110 + 2, y=legend_y, text="accessory", font_size=8),
                svg.write_arrow(x1=legend_x+12, y1=legend_y+15, x2=legend_x+12, y2=legend_y+15, color="#000", stroke_width=1),
                svg.write_text_left(x=legend_x + 20, y=legend_y+15, text="n proteins encoded per 500kbp window", font_size=8)
            ])
        else:
            legend_width = 300
            self.body.extend([
                svg.write_legend_rect(x=legend_x, y=self.y_start+25, width=legend_width, height=40),
                #
                svg.write_line(x1=legend_x+ 10, y1=legend_y, x2=legend_x + 20, y2=legend_y, color="#eeeeee", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2, y=legend_y, text="core", font_size=8),
                svg.write_line(x1=legend_x + 20 + 2 + 100, y1=legend_y, x2=legend_x + 20 + 2 + 110, y2=legend_y, color="#bbbbbb", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2 + 110 + 2, y=legend_y, text="accessory", font_size=8),
                #
                svg.write_line(x1=legend_x + 10, y1=legend_y + 15, x2=legend_x + 20, y2=legend_y + 15, color="#cc0000", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2, y=legend_y + 15, text="early (1,2,3 dpi)", font_size=8),
                #
                svg.write_line(x1=legend_x + 20 + 2 + 100, y1=legend_y + 15, x2=legend_x + 20 + 2 + 110, y2=legend_y + 15, color="#0000cc", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2 + 110 + 2, y=legend_y + 15, text="late (7 dpi)"),
                #
                svg.write_line(x1=legend_x + 20 + 2 + 200, y1=legend_y + 15, x2=legend_x + 20 + 2 + 210, y2=legend_y + 15, color="#00cc00", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2 + 210 + 2, y=legend_y + 15, text="all", font_size=8)
            ])
    

class ChromPlotPortrait(ChromPlot):
    y_start = 50
    width = 650
    height = 650

    def __init__(self, regions, **kwargs):
        super().__init__(regions, **kwargs)
        print("P")

    def calc_accessory_region(self, chrom, region_length, offset):
         acc_height_factor = 7/8 if chrom == 0 else 6/8
         acc_height = (1e6 * acc_height_factor) / self.max_region_length * ChromPlot.max_axis_length        
         acc_y = offset + region_length - acc_height
         return acc_height, acc_y

    def draw_regions(self):
        y_offset = 35
        for i, (region_id, region) in enumerate(self.regions.items()):
            print(region_id, region)        
            region_length = int(region["length"] / self.max_region_length * ChromPlot.max_axis_length)

            y = self.y_start + y_offset * i
   
            self.body.extend([ 
                '<g id="{gid}">\n'.format(gid=region["chrom"]),
                svg.write_rect(x=50, y=y, height=10, width=region_length, fill=region["chr_color"], stroke_width=1.0, rx=3),
                svg.write_text_right(x=40, y=y+5, text="Chr{chrom}".format(**region), font_size=10)
            ])

            # deal with accessory regions on Chr1, Chr2
            if i < 2:
                acc_height, acc_y = self.calc_accessory_region(i, region_length, 50)
                self.body.extend([
                    svg.write_rect(x=acc_y, y=y, width=acc_height, height=10, fill="#bbbbbb", stroke_width=0.0, rx=3),
                    svg.write_dashed_line(x1=acc_y, x2=acc_y, y1=y, y2=y+10, color="#000000", stroke_width="1.0")
                ])

            for gene, pos in self.genes.get(region_id, list()):
                print(region)
                x_pos = 50 + int(pos / region["length"] * region_length)
                # y_pos = y_start - int(pos / region["length"] * region_length)
                col = self.gene_groups.get(gene, None)
                if not col:
                    print("xxx", gene)
                else:
                    if col == "#00cc00":
                        col = self.expression_map.get(gene, col)
                        print(gene, col)
                    #body += write_line(x1=x, y1=y_pos, x2=x+10, y2=y_pos, color=col, stroke_width=1.5)
                    self.body.append(svg.write_line(x1=x_pos, y1=y, x2=x_pos, y2=y+10, color=col, stroke_width=1.5))

            binsize = 500000
            bins = self.calc_protein_bins(region_id, binsize=binsize)
            for pos, count in bins.items():
                x_pos = 50 + int((pos * binsize) / region["length"] * region_length)
                self.body.extend([
                    svg.write_arrow(x1=x_pos, y1=y-3, x2=x_pos, y2=y-5, color="#000", stroke_width=1),
                    svg.write_text(x=x_pos, y=y-10, text=count, font_size=10)
                ])

            # overpaint chromosome borders
            self.body.extend([
                svg.write_rect(x=50, y=y, width=region_length, height=10, fill="none", stroke_width=1.5, rx=3),
                "</g>\n"
            ])

    def draw_legend(self):
        legend_x = 300
        legend_y = self.y_start + 35 * (14) - 30 - 20
        legend_width = 200

        if self.proteins:
            self.body.extend([
                svg.write_legend_rect(x=legend_x, y=legend_y, width=legend_width, height=60),
                svg.write_line(x1=legend_x+ 10, y1=legend_y +10, x2=legend_x + 20, y2=legend_y +10, color="#eeeeee", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2, y=legend_y+10, text="core", font_size=8),
                svg.write_line(x1=legend_x + 20 + 2 + 100, y1=legend_y+10, x2=legend_x + 20 + 2 + 110, y2=legend_y+10, color="#bbbbbb", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2 + 110 + 2, y=legend_y+10, text="accessory", font_size=8),
                svg.write_arrow(x1=legend_x+12, y1=legend_y+16 + 8, x2=legend_x+12, y2=legend_y+14 + 8, color="#000", stroke_width=1),
                svg.write_text_left(x=legend_x + 20, y=legend_y+15 + 10, text="n proteins encoded per 500kbp window", font_size=8),
                svg.write_line(x1=legend_x + legend_width - self.scale_x_length - 10, y1=legend_y + 40, x2=legend_x + legend_width - 10, y2=legend_y + 40, color="#000000", stroke_width=1),
                svg.write_line(x1=legend_x + legend_width - self.scale_x_length - 10, y1=legend_y + 35, x2=legend_x + legend_width - self.scale_x_length - 10, y2=legend_y+45, color="#000000", stroke_width=1),
                svg.write_line(x1=legend_x + legend_width - 10, y1=legend_y + 35, x2=legend_x + legend_width - 10, y2=legend_y+45, color="#000000", stroke_width=1),
                svg.write_text(x=(legend_x + legend_width - self.scale_x_length - 10) + self.scale_x_length/2, y=legend_y + 50, text="1 Mbp", font_size=8),
            ])
        else:
            self.body.extend([
                svg.write_legend_rect(x=legend_x, y=legend_y, width=legend_width, height=60),
                svg.write_line(x1=legend_x+ 10, y1=legend_y +10, x2=legend_x + 20, y2=legend_y +10, color="#eeeeee", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2, y=legend_y+10, text="core", font_size=8),
                svg.write_line(x1=legend_x + 20 + 2 + 100, y1=legend_y+10, x2=legend_x + 20 + 2 + 110, y2=legend_y+10, color="#bbbbbb", stroke_width=1.5),                
                svg.write_text_left(x=legend_x + 20 + 2 + 110 + 2, y=legend_y+10, text="accessory", font_size=8),
                #
                svg.write_line(x1=legend_x + 10, y1=legend_y + 25, x2=legend_x + 20, y2=legend_y + 25, color="#cc0000", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2, y=legend_y + 25, text="early (1,2,3 dpi)", font_size=8),
                #
                svg.write_line(x1=legend_x + 10, y1=legend_y + 35, x2=legend_x + 20, y2=legend_y + 35, color="#0000cc", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2, y=legend_y + 35, text="late (7 dpi)", font_size=8),
                #
                svg.write_line(x1=legend_x + 10, y1=legend_y + 45, x2=legend_x + 20, y2=legend_y + 45, color="#00cc00", stroke_width=1.5),
                svg.write_text_left(x=legend_x + 20 + 2, y=legend_y + 45, text="all", font_size=8),
                #
                svg.write_line(x1=legend_x + legend_width - self.scale_x_length - 10, y1=legend_y + 40, x2=legend_x + legend_width - 10, y2=legend_y + 40, color="#000000", stroke_width=1),
                svg.write_line(x1=legend_x + legend_width - self.scale_x_length - 10, y1=legend_y + 35, x2=legend_x + legend_width - self.scale_x_length - 10, y2=legend_y+45, color="#000000", stroke_width=1),
                svg.write_line(x1=legend_x + legend_width - 10, y1=legend_y + 35, x2=legend_x + legend_width - 10, y2=legend_y+45, color="#000000", stroke_width=1),
                svg.write_text(x=(legend_x + legend_width - self.scale_x_length - 10) + self.scale_x_length/2, y=legend_y + 50, text="1 Mbp", font_size=8),

            ])
