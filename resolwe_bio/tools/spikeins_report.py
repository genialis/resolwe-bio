#!/usr/bin/env python3
"""Create spikeins-qc report."""
import argparse
import datetime
import os

from jinja2 import Environment, FileSystemLoader

parser = argparse.ArgumentParser(description="Fill data into html template file.")
parser.add_argument('--sample_names', help="Sample names.", nargs='+')
parser.add_argument('--images_dir', help="List of files to include in report.", required=True)
parser.add_argument('--template', help="Report template (*.html file).", required=True)


if __name__ == '__main__':
    args = parser.parse_args()

    ercc_imgs = []
    images = [os.path.join(args.images_dir, name) for name in os.listdir(args.images_dir)]
    for sample_name in args.sample_names:
        # Find the corresponding image
        expected_name = "{} (ERCC spike-in's).png".format(sample_name)
        image = next((img for img in images if os.path.basename(img) == expected_name), None)
        if image:
            ercc_imgs.append({
                'path': os.path.relpath(os.path.abspath(image), start=os.getcwd()),
            })
        else:
            ercc_imgs.append({
                'sample_name': sample_name,
            })

    content = {
        'sample_names': args.sample_names,
        'today': datetime.datetime.now().date().isoformat(),
        'ercc_imgs': ercc_imgs,
    }

    # Prepare template and fill with content.
    template_dir = os.path.dirname(os.path.abspath(args.template))
    template_basename = os.path.basename(args.template)
    env = Environment(loader=FileSystemLoader(template_dir))
    with open('report.html', 'wt', encoding='utf-8') as ofile:
        template = env.get_template(template_basename)
        ofile.write(template.render(**content))
