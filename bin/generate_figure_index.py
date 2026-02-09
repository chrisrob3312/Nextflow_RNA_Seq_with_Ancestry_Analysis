#!/usr/bin/env python3
"""Generate HTML index page for all pipeline figures."""

import argparse
import os
import glob

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--figure-dir', default='figures')
    parser.add_argument('--output', default='figure_index.html')
    args = parser.parse_args()

    figures = sorted(glob.glob(os.path.join(args.figure_dir, '**/*.png'), recursive=True))

    html = """<!DOCTYPE html>
<html><head><title>Pipeline Figures</title>
<style>
body { font-family: sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }
.gallery { display: grid; grid-template-columns: repeat(auto-fill, minmax(400px, 1fr)); gap: 20px; }
.figure { border: 1px solid #ddd; padding: 10px; border-radius: 5px; }
.figure img { width: 100%; }
.figure p { text-align: center; font-size: 13px; color: #555; }
h2 { color: #2c3e50; border-bottom: 1px solid #ccc; }
</style></head><body>
<h1>Pipeline Figure Gallery</h1>
"""

    # Group by subdirectory
    groups = {}
    for fig in figures:
        rel = os.path.relpath(fig, args.figure_dir)
        parts = rel.split(os.sep)
        group = parts[0] if len(parts) > 1 else 'General'
        groups.setdefault(group, []).append(fig)

    for group, figs in sorted(groups.items()):
        html += f'<h2>{group.replace("_", " ").title()}</h2>\n<div class="gallery">\n'
        for fig in figs:
            name = os.path.basename(fig).replace('.png', '').replace('_', ' ')
            html += f'<div class="figure"><img src="{fig}" alt="{name}"><p>{name}</p></div>\n'
        html += '</div>\n'

    html += '</body></html>'

    with open(args.output, 'w') as f:
        f.write(html)

if __name__ == '__main__':
    main()
