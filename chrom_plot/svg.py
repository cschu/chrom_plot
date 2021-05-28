import sys
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
    #return '<text x="{x}" y="{y}" font-size="10" font-family="Arial, sans-serif" text-anchor="middle">{text}</text>\n'.format(**args)
    print("FONT_SIZE: ", args["font_size"])
    return '<text x="{x}" y="{y}" font-size="{font_size}" font-family="Arial, sans-serif" text-anchor="middle" dominant-baseline="central">{text}</text>\n'.format(**args)
def write_text_right(**args):
    return '<text x="{x}" y="{y}" font-size="{font_size}" font-family="Arial, sans-serif" text-anchor="end" dominant-baseline="central">{text}</text>\n'.format(**args)
def write_text_left(**args):
    return '<text x="{x}" y="{y}" font-size="{font_size}" font-family="Arial, sans-serif" dy="0.25em">{text}</text>\n'.format(**args)
def write_text_rotated(**args):
    #return '<text x="{x}" y="{y}" font-size="10" font-family="Arial, sans-serif" text-anchor="middle">{text}</text>\n'.format(**args)
    return '<text x="{x}" y="{y}" font-size="{font_size}" font-family="Arial, sans-serif" text-anchor="middle" dominant-baseline="central" transform="rotate(90)">{text}</text>\n'.format(**args)

def write_circle(**args):
    return '<circle cx="{x}" cy="{y}" r="{radius}" fill="{fill}"/>'.format(**args)

