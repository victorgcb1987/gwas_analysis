from pathlib import Path

import numpy

import qmplot

from src import normalization
from src import plink
from src import plot

from src.genome_coord_transform import (
    GenomeCoordinateConverter,
    get_genome_sizes,
)


def create_phenotype_from_df(df, trait):
    print(df)
    names = [name for key, name in df["SAMPLE_NAME"].items()]
    values = [value for key, value in df[trait].items()]
    return {accesion[0]: accesion[1] for accesion in zip(names, values)}



def _do_gwas_analysis(
    bfiles_base_path,
    phenotype_dframe,
    out_dir,
    out_base_name,
    qualitative,
    genome_fai_path,
    covars_path=None,
    traits=None,
    desired_accs=None,
):
    out_dir.mkdir(exist_ok=True, parents=True)

    if traits is None:
        traits = list(phenotype_dframe.keys())

    for trait in traits:
        print(f"Doing GWAS for trait: {trait}")
        trait_out_dir = out_dir / trait
        trait_out_dir.mkdir(exist_ok=True, parents=True)

        phenotypes = create_phenotype_from_df(phenotype_dframe, trait)

        phenotypes_path = trait_out_dir / "phenotypes.pheno"
        with phenotypes_path.open("wt") as fhand:
            plink.write_phenotype_file(phenotypes, fhand, quantitative=qualitative)

        if qualitative:
            test_type = "logistic"
        else:
            test_type = "linear"

        out_base_path = trait_out_dir / f"{bfiles_base_path.name}.{trait}"
        kwargs = {
            "bfiles_base_path": bfiles_base_path,
            "phenotypes_path": phenotypes_path,
            "test_type": test_type,
            "out_base_path": out_base_path,
        }
        if covars_path:
            kwargs["covars_path"] = covars_path
        else:
            kwargs["allow_no_covars"] = True
        res = plink.do_gwas(**kwargs)

        pvals = [
            ("pval", "non_adjusted_pval"),
            ("sidak_step_down_pval", "sidak_step_down_pval"),
            ("benjamini_yekutieli_pval", "benjamini_yekutieli_pval"),
        ]

        pvalues_dframe = res.get("adjusted_pvalues")
        if pvalues_dframe is None:
            print(f"{trait}: Zero valid tests")
            continue

        for pval, pval_tag in pvals:
            pvalues = pvalues_dframe[pval]

            chroms, poss = zip(*[snp_id.split(":") for snp_id in pvalues.index])
            chroms = numpy.array(chroms)
            poss = numpy.array(poss, dtype=int)

            index_to_sort = chroms.argsort()
            chroms = chroms[index_to_sort]
            poss = poss[index_to_sort]
            pvalues = pvalues[index_to_sort]

            if pval == "pval":
                fig_qq = plot.SimpleFigure()
                qmplot.qqplot(data=pvalues, title=trait, ax=fig_qq.axes)
                fig_qq.save_fig(
                    trait_out_dir / f"{out_base_name}_{trait}_{pval_tag}_qq_plot.png"
                )

            coord_converter = GenomeCoordinateConverter(
                chrom_lens=get_genome_sizes(genome_fai_path=genome_fai_path)
            )
            fig = plot.SimpleFigure(fig_size=(10, 6))
            log_pvalues = -numpy.log10(pvalues)
            plot.plot_values_along_genome(
                log_pvalues,
                chroms,
                poss,
                coord_converter=coord_converter,
                axes=fig.axes,
                chroms_are_sorted=True,
            )
            fig.axes.set_title(trait)
            fig.axes.set_ylabel(f"-log10({pval_tag})")
            fig.save_fig(
                trait_out_dir / f"{out_base_name}_{trait}_{pval_tag}_along_genome.png"
            )

        csv_path = trait_out_dir / f"{out_base_name}_{trait}_pvalues_along_genome.csv"
        pvalues = pvalues_dframe
        last_row_idx = pvalues.iloc[99:100, 0].index[0]
        some_pvalues = pvalues_dframe.loc[:last_row_idx, [pval for pval, _ in pvals]]
        some_pvalues.to_csv(csv_path)


def do_gwas_analysis(
    bfiles_base_path,
    phenotype_dframe,
    out_base_name,
    out_dir,
    genome_fai_path,
    qualitative,
    covars_path=None,
    traits=None,
    desired_accs=None,
):
    return _do_gwas_analysis(
        bfiles_base_path,
        phenotype_dframe,
        out_base_name=out_base_name,
        out_dir=out_dir,
        traits=traits,
        desired_accs=desired_accs,
        covars_path=covars_path,
        genome_fai_path=genome_fai_path,
        qualitative=qualitative,
    )


# def _filter_traits_data(traits_data, selected_accs):
#     if selected_accs is None:
#         return traits_data

#     selected_accs = set(selected_accs)

#     accs = selected_accs.intersection(traits_data.index)
#     traits_data = traits_data.loc[accs, :]
#     return traits_data
# from pandas import DataFrame