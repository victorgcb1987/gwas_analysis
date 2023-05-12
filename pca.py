from pathlib import Path
import argparse

from src.plink import VariantFilters, do_pca 



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
    help_pca = "(Required) PCA output path"
    parser.add_argument("--PCA", "-c",
                        type=str, help=help_pca,
                        required=True)
    help_max_missing_rate = "(Optional) maximun variant missing rate"
    parser.add_argument("--max_missing_rate", "-m",
                        type=float, help=help_max_missing_rate,
                        default=0.1)
    return parser


def get_options():
    parser = parse_arguments()
    options = parser.parse_args()
    base_plink_path = Path(options.plink)
    pruned_vars = Path(options.pruned_vars)
    pruned_plink_path = Path(options.pruned_plink)
    max_missing_rate = options.max_missing_rate
    pca_base_path = Path(options.PCA)
    return {'base_plink_path': base_plink_path,
            'pruned_vars': pruned_vars,
            "pruned_plink_path": pruned_plink_path,
            "max_missing_rate": max_missing_rate,
            "pca_base_path": pca_base_path}

if __name__ == '__main__':
    options = get_options()
    pca_variant_filters = VariantFilters(
            max_missing_rate=options["max_missing_rate"],
            lists_of_vars_to_keep_paths=options["pruned_vars"],
        )
    pca_res = do_pca(
            options["base_plink_path"],
            out_base_path=options["pca_base_path"],
            variant_filters=pca_variant_filters
        )