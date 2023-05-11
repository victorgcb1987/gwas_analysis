import argparse
import subprocess
import sys
from pathlib import Path

from config import EXECUTABLES_REQUIREMENTS as exec_recs
from dependencies import get_executables

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

#Parse and return values given to options when running this program
def get_options():
    parser = parse_arguments()
    options = parser.parse_args()
    vcf_fpath = Path(options.vcf)
    output_dir = Path(options.out)
    return {'vcf_fpath': vcf_fpath,
            'out_dir': output_dir}


def create_plink_bcfile(vcf_path, base_path):

    cmd = [get_executables(exec_recs["plink2"])]
    cmd.extend(['--vcf', str(vcf_path)])
    cmd.extend(['--out', str(base_path)])
    cmd.extend(['--allow-extra-chr', '--double-id', '--vcf-half-call', 'missing',
                '--set-missing-var-ids', '@:# ', '--make-bed'])
    stderr_path = Path(str(base_path) + '.bfiles.stderr')
    stdout_path = Path(str(base_path) + '.bfiles.stderr')

    print('Running: ', ' '.join(cmd))
    subprocess.run(cmd, stdout=stdout_path.open('wt'),
                   stderr=stderr_path.open('wt'), check=True)

if __name__ == '__main__':
    options = get_options()
    vcf_path = options["vcf_fpath"]
    base_path = options["vcf_fpath"].parent / options["vcf_fpath"].stem
    create_plink_bcfile(vcf_path, base_path)