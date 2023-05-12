import argparse
from pathlib import Path

from src.plink import create_ld_indep_variants_file, VariantFilters





def parse_arguments():
    desc = "Create PLINK format file for GWAS"
    parser = argparse.ArgumentParser(description=desc)
    help_plink_input = "(Required) PLINK basename path"
    parser.add_argument("--plink", 
                        "-i", type=str,
                        help=help_plink_input,
                        required=True)
    help_pruned_vars = "(Required) Pruned variations path"
    parser.add_argument("--pruned_vars", "-p",
                        type=str, help=help_pruned_vars,
                        required=True)
    help_plink_output = "(Required) Pruned PLINK output path"
    parser.add_argument("--pruned_plink", "-o",
                        type=str, help=help_plink_output,
                        required=True)
    help_max_missing_rate = "(Optional) maximun variant missing rate"
    parser.add_argument("--max_missing_rate", "-m",
                        type=float, help=help_max_missing_rate,
                        default=0.1)
    help_windows_size = "(Optional), windows size for removing LD"
    parser.add_argument("--windows_size", "-w",
                        type=int, help=help_windows_size,
                        default=50)
    help_step_size = "(Optional) windows steps"
    parser.add_argument("--step_size", "-s",
                        type=int, help=help_step_size,
                        default=5)
    help_bad_ld = "(Optional) allow --bad-ld in plink."
    parser.add_argument("--bad-ld", "-b", type=bool, 
                        action="store_true", help=help_bad_ld)
    return parser


def get_options():
    parser = parse_arguments()
    options = parser.parse_args()
    base_plink_path = Path(options.plink)
    pruned_vars = Path(options.pruned_vars)
    pruned_plink_path = Path(options.pruned_plink)
    max_missing_rate = options.max_missing_rate
    windows_size = options.windows_size
    step_size = options.step_size
    bad_ld = options.bad-ld
    return {'base_plink_path': base_plink_path,
            'pruned_vars': pruned_vars,
            "pruned_plink_path": pruned_plink_path,
            "max_missing_rate": max_missing_rate,
            "windows_size": windows_size,
            "step_size": step_size,
            "bad_ld": bad_ld}


if __name__ == '__main__':
    options = get_options()
    create_ld_indep_variants_file(bfiles_base_path=options["base_plink_path"],
            out_base_path=options["pruned_plink_path"],
            pruned_vars_list_path=options["pruned_vars"],
            variant_filters=VariantFilters(max_missing_rate=options["max_missing_rate"]),
            window_size=options["windows_size"],
            step_size=options["step_size"],
            bad_ld=options["bad_ld"]
        )