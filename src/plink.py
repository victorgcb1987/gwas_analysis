import os
import subprocess
from pathlib import Path

from src.config import EXECUTABLES_REQUIREMENTS as exec_recs
from src.dependencies import get_executables


class VariantFilters:
    def __init__(
        self,
        max_major_freq=0.01,
        max_missing_rate=0.1,
        lists_of_vars_to_keep_paths=None,
    ):
        self.max_major_freq = max_major_freq
        self.max_missing_rate = max_missing_rate

        if lists_of_vars_to_keep_paths is None:
            lists_of_vars_to_keep_paths = []
        self.lists_of_vars_to_keep_paths = lists_of_vars_to_keep_paths

    def create_cmd_arg_list(self):
        args = []
        if self.max_major_freq is not None:
            args.extend(["--maf", str(self.max_major_freq)])
        if self.max_missing_rate is not None:
            args.extend(["--geno", "0.1"])
        if self.lists_of_vars_to_keep_paths:
            args.append("--extract")
            args.extend(map(str, self.lists_of_vars_to_keep_paths))
        return args


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

def run_cmd(cmd, stdout_path, stderr_path):
    stdout_fhand = stdout_path.open("wt")
    stderr_fhand = stderr_path.open("wt")
    try:
        subprocess.run(
            cmd,
            stdout=stdout_fhand,
            stderr=stderr_fhand,
            check=True,
        )
    except subprocess.CalledProcessError:
        print("stdout")
        stdout_fhand = stdout_path.open("rt")
        print(stdout_fhand.read())
        stdout_fhand.close()
        print("stderr")
        stderr_fhand = stderr_path.open("rt")
        print(stderr_fhand.read())
        stderr_fhand.close()
        raise
    stderr_fhand.close()
    stdout_fhand.close()


def create_ld_indep_variants_file(
    bfiles_base_path,
    out_base_path,
    pruned_vars_list_path,
    variant_filters: VariantFilters,
    window_size=50,
    step_size=5,
    sizes_are_in_number_of_vars=True,
    r2_threshold=0.5,
    bad_ld=False
):

    pruned_vars_list_path = Path(pruned_vars_list_path)

    if not variant_filters.max_major_freq:
        raise ValueError("It is really important to filter out the low freq variants")

    current_working_dir = os.getcwd()
    working_dir = out_base_path if out_base_path.is_dir() else out_base_path.parent
    os.chdir(working_dir)

    cmd = [get_executables(exec_recs["plink2"])]
    cmd.extend(["--bfile", str(bfiles_base_path)])

    cmd.append("--allow-extra-chr")
    # Note: --allow-no-sex no longer has any effect.  (Missing-sex samples are
    # automatically excluded from association analysis when sex is a covariate, and
    # treated normally otherwise.)
    # cmd.append("--allow-no-sex")

    if variant_filters is not None:
        cmd.extend(variant_filters.create_cmd_arg_list())

    cmd.append("--indep-pairwise")
    if bad_ld:
        cmd.extend("--bad-ld")
    if sizes_are_in_number_of_vars:
        cmd.extend([str(window_size), str(step_size), str(r2_threshold)])
    else:
        cmd.extend([f"{window_size}kb", str(step_size), str(r2_threshold)])

    out_base_path = Path(str(out_base_path) + ".ld")

    stderr_path = Path(str(out_base_path) + ".stderr")
    stdout_path = Path(str(out_base_path) + ".stdout")

    run_cmd(cmd, stdout_path, stderr_path)

    plink_out_path = working_dir / "plink2.prune.in"
    if not pruned_vars_list_path.exists() or not os.path.samefile(
        plink_out_path, pruned_vars_list_path
    ):
        os.rename(plink_out_path, pruned_vars_list_path)

    os.chdir(current_working_dir)