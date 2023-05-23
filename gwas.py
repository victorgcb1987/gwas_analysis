import argparse
from pathlib import Path

import pandas as pd

from src.gwas import do_gwas_analysis
from src import normalization 

def parse_arguments():
    desc = "Create PLINK format file for GWAS"
    parser = argparse.ArgumentParser(description=desc)
    help_plink_input = "(Required) PLINK basename path"
    parser.add_argument("--plink", 
                        "-i", type=str,
                        help=help_plink_input,
                        required=True)
    help_traits = "(Required) traits path"
    parser.add_argument("--traits", 
                        "-t", type=str,
                        help=help_traits,
                        required=True)
    help_normalization = "(Optional) traits normalization method. Available methods are yeo-johnson, quantile"
    parser.add_argument("--normalization", 
                        "-n", type=str,
                        help=help_normalization,
                        default="None")
    help_pca = "(Optional) Add PCA poblation structure to GWAs analysis"
    parser.add_argument("--pca",
                        "-p", type=str, help=help_pca,
                        default="False")
    help_output = "(Required) GWAs output path"
    parser.add_argument("--out", "-o",
                        type=str, help=help_output,
                        required=True)
    help_base_output = "(Required) GWAs output base name"
    parser.add_argument("--out", "-o",
                        type=str, help=help_base_output,
                        required=True)
    help_faidx = "(Required) genome faidx file"
    parser.add_argument("--faidx", "-f",
                        type=str, help=help_faidx,
                        required=True)
    help_traits_to_analyze = "(Optional) traits to analyze"
    parser.add_argument("--traits", "-t",
                        help=help_traits_to_analyze,
                        nargs="*",required=False)
    help_qualitative = "(Optional) set traits to analyze as qualitative"
    parser.add_argument("--qualitative", "-q",
                        help=help_qualitative,
                        action="store_true")
    

def get_options():
    parser = parse_arguments()
    options = parser.parse_args()
    base_plink_path = Path(options.plink)
    traits_path = Path(options.traits)
    normalization_method = options.normalization
    pca_structure = options.pca
    if pca_structure == "False":
        pca_structure = False
    else:
        pca_structure = Path(pca_structure)
    if normalization_method not in ["yeo-johnson", "quantile", "None"]:
        raise ValueError ("Normalization method not available: {}".format(normalization_method))
    
    output_path = Path(options.out)
    output_base = options.base
    faidx = Path(options.faidx)
    traits = options["traits"]
    is_qualitative = options["is_qualitative"]
    return {'base_plink_path': base_plink_path,
            'traits_path': traits_path,
            'normalization_method': normalization_method,
            "output_path": output_path,
            "output_base": output_base,
            "faidx": faidx,
            "traits": traits,
            "is_qualitative": is_qualitative
            }


if __name__ == "__main__":
    options = get_options()
    phenotype_dframe = pd.read_csv(options["traits_path"], sep="\t")
    normalization_method = options["normalization_method"]
    output_dir = options["output_path"]
    if not output_dir.exists():
         output_dir.mkdir()
    if normalization_method:
        if normalization_method == "yeo_johnson":
            phenotype_dframe = normalization.normalize_dframe(
            phenotype_dframe, method=normalization_method)
        elif normalization_method == "quantile":
            phenotype_dframe = normalization.quantile_transform_dframe(
                phenotype_dframe, n_quantiles=10
                )
    gwas_kwargs = {
            "bfiles_base_path": options["base_plink_path"],
            "phenotype_dframe": phenotype_dframe,
            "out_base_name": options["out_base"],
            "out_dir": options["out_path"],
            "traits": options["traits"],
            "genome_fai_path": options["faidx"],
            "qualitative": options["is_qualitative"],
        }
    if options["pca"]:
            gwas_kwargs["covars_path"] = options["pca"]
    do_gwas_analysis(**gwas_kwargs)
