import argparse
from pathlib import Path

from src.plink import create_plink_bcfile


def parse_arguments():
    desc = "Create PLINK format file for GWAS"
    parser = argparse.ArgumentParser(description=desc)
    help_calling = "(Required) VCF File."
    parser.add_argument("--vcf", 
                        "-v", type=str,
                        help=help_calling,
                        required=True)
    help_output = "(Required) Output dir"
    parser.add_argument("--out", "-o",
                        type=str, help=help_output,
                        required=True)
    return parser


def get_options():
    parser = parse_arguments()
    options = parser.parse_args()
    vcf_fpath = Path(options.vcf)
    output_dir = Path(options.out)
    return {'vcf_fpath': vcf_fpath,
            'out_dir': output_dir}


if __name__ == '__main__':
    options = get_options()
    vcf_path = options["vcf_fpath"]
    base_path = options["vcf_fpath"].parent / options["vcf_fpath"].stem.replace(".vcf", "")
    create_plink_bcfile(vcf_path, base_path)