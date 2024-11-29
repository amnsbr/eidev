import time
import numpy as np
from cubnm import optimize, _setup_opts
import os
import json
import argparse

INPUT_DIR = os.environ["INPUT_DIR"]
OUTPUT_DIR = os.environ["OUTPUT_DIR"]

# ensure cubnm is installed with the CUBNM_NOISE_WHOLE env variable set to 1
assert not _setup_opts.noise_segment_flag, (
    "cubnm must be installed with the CUBNM_NOISE_WHOLE env variable set to 1"
)

def get_paths(fc_sub, sc_sub, fc_ses, sc_ses, args):
    if sc_sub == "fc_sub":
        sc_sub = fc_sub
    if sc_sub == "micamics":
        sc_path = f"{INPUT_DIR}/micamics/SC/group-all/ctx_parc-{args.parc}_approach-median_{args.sc_config}_desc-strength.txt"
        if args.dataset == "pnc":
            out_path = f"{OUTPUT_DIR}/{args.dataset}/group-micamics/ctx_parc-{args.parc}_{args.sc_config}/{fc_sub}"
        else:
            raise NotImplementedError
    else:
        if args.dataset == "pnc":
            out_path = f"{OUTPUT_DIR}/{args.dataset}/{sc_sub}/ctx_parc-{args.parc}_{args.sc_config}"
            sc_path = f"{INPUT_DIR}/{args.dataset}/SC/{sc_sub}/ctx_parc-{args.parc}_{args.sc_config}_thresh-1_desc-strength.txt"
        else:
            out_path = f"{OUTPUT_DIR}/{args.dataset}/{fc_ses}/{sc_sub}/sc-{sc_ses}_ctx_parc-{args.parc}_{args.sc_config}"
            sc_path = f"{INPUT_DIR}/{args.dataset}/{sc_ses}/SC/{sc_sub}/ctx_parc-{args.parc}_{args.sc_config}_thresh-1_desc-strength.txt"
    out_path = f"{out_path}/{args.maps_name}"
    # specify fc and fcd path
    fc_prefix = f"ctx_parc-{args.parc}_hemi-LR_highpass-013_lowpass-none"
    if args.exc_inter:
        fc_prefix = f"{fc_prefix}_exc-inter"
    if args.dataset == "pnc":
        fc_path = f"{INPUT_DIR}/{args.dataset}/FC/{fc_sub}/{fc_prefix}_desc-FCtril.txt"
        fcd_path = (
            f"{INPUT_DIR}/{args.dataset}/FC/{fc_sub}/{fc_prefix}_desc-FCDtril.txt"
        )
    elif args.dataset == "imagen":
        fc_path = f"{INPUT_DIR}/{args.dataset}/{fc_ses}/FC/{fc_sub}/{fc_prefix}_desc-FCtril.txt"
        fcd_path = f"{INPUT_DIR}/{args.dataset}/{fc_ses}/FC/{fc_sub}/{fc_prefix}_desc-FCDtril.txt"
    return fc_path, fcd_path, sc_path, out_path


def get_subs(args):
    # get list of all subjects/groups
    if args.dataset == "pnc":
        if args.fc_level == "sub":
            if args.subsample:
                all_subs = np.loadtxt(f"{INPUT_DIR}/pnc_subsample_200.txt", dtype=str)
            else:
                all_subs = np.loadtxt(f"{INPUT_DIR}/pnc_subs.txt", dtype=str)
        elif args.fc_level == "age_group":
            all_subs = np.loadtxt(f"{INPUT_DIR}/pnc_age_groups.txt", dtype=str)
    elif args.dataset == "imagen":
        if args.sc_level == "ses":
            all_subs = np.loadtxt(f"{INPUT_DIR}/imagen_subs_BLnFU2.txt", dtype=str)
            args.sc_ses = args.fc_ses
        elif args.sc_level == "sub":
            all_subs = np.loadtxt(f"{INPUT_DIR}/imagen_subs_FU2.txt", dtype=str)
            args.sc_ses = "FU2"
    # loop through subjects/groups and select the ones
    # that have the input data, are not excluded
    # and CMAES for them has not been completed (all runs)
    selected_subs = []
    for fc_sub in all_subs:
        if os.path.exists(f"{INPUT_DIR}/{args.dataset}/excluded/{fc_sub}.txt"):
            continue
        fc_path, fcd_path, sc_path, out_path = get_paths(
            fc_sub, args.sc_sub, args.fc_ses, args.sc_ses, args
        )
        in_exist = (
            os.path.exists(fc_path)
            and os.path.exists(fcd_path)
            and os.path.exists(sc_path)
        )
        # check if n runs of CMAES with the same options are done
        # for current subject/group
        seeds_exist = args.n_runs * [False]
        if os.path.exists(out_path):
            for seed in range(1, args.n_runs + 1):
                for subdir in os.listdir(out_path):
                    if os.path.exists(f"{out_path}/{subdir}/optimizer.json"):
                        with open(f"{out_path}/{subdir}/optimizer.json", "r") as f:
                            optimizer_conf = json.load(f)
                        if (
                            (optimizer_conf.get("seed", None) == seed)
                            and (optimizer_conf.get("popsize", None) == args.popsize)
                            and (optimizer_conf.get("n_iter", None) == args.n_iter)
                        ):
                            seeds_exist[seed - 1] = True
        out_exist = all(seeds_exist)
        if in_exist and not out_exist:
            selected_subs.append(fc_sub)
        if len(selected_subs) == args.n_subs:
            break
    print(*selected_subs)


def run(args):
    # specify in and out paths
    if args.dataset == "imagen":
        if args.sc_level == "ses":
            args.sc_ses = args.fc_ses
        elif args.sc_level == "sub":
            args.sc_ses = "FU2"
    # specify heterogeneity
    maps_path = None
    node_grouping = None
    het_params = ["wEE", "wEI"]
    if args.maps_name == "homo":
        het_params = []
    elif args.maps_name == "node":
        node_grouping = "node"
    else:
        maps_path = os.path.join(INPUT_DIR, f"{args.maps_name}_{args.parc}_zscore.txt")
    # set up problem objects for all subjects
    problems = []
    for fc_sub in args.fc_subs:
        fc_path, fcd_path, sc_path, out_path = get_paths(
            fc_sub, args.sc_sub, args.fc_ses, args.sc_ses, args
        )
        # skip this CMA-ES run if it's already done
        skip = False
        if os.path.exists(out_path):
            for subdir in os.listdir(out_path):
                if os.path.exists(f"{out_path}/{subdir}/optimizer.json"):
                    with open(f"{out_path}/{subdir}/optimizer.json", "r") as f:
                        optimizer_conf = json.load(f)
                    if (
                        (optimizer_conf.get("seed", None) == args.opt_seed)
                        and (optimizer_conf.get("popsize", None) == args.popsize)
                        and (optimizer_conf.get("n_iter", None) == args.n_iter)
                    ):
                        print(
                            f"Skipping {fc_sub} seed {args.opt_seed} because it was already run"
                        )
                        skip = True
        if skip:
            continue
        print(
            f"Running {fc_sub} seed {args.opt_seed} with:",
            "\n\tfc_path:",
            fc_path,
            "\n\tfcd_path:",
            fcd_path,
            "\n\tsc_path:",
            sc_path,
            "\n\tout_path:",
            out_path,
        )

        problem = optimize.BNMProblem(
            model='rWW',
            params={
                "G": (0.5, 4.0),
                "wEE": (0.05, 0.75),
                "wEI": (0.05, 0.75),
            },
            het_params=het_params,
            max_fic_trials=args.max_fic_trials,
            duration=450,
            TR=3,
            window_size=10,
            window_step=2,
            sc=sc_path,
            emp_fc_tril=np.loadtxt(fc_path),
            emp_fcd_tril=np.loadtxt(fcd_path),
            gof_terms=['+fc_corr', '-fc_diff', '-fcd_ks'],
            maps=maps_path,
            node_grouping=node_grouping,
            out_dir=out_path,
            bw_params="heinzle2016-3T",
            exc_interhemispheric=args.exc_inter,
            rand_seed=410,
        )
        problems.append(problem)
    cmaes = optimize.CMAESOptimizer(
        popsize=args.popsize,
        n_iter=args.n_iter,
        seed=args.opt_seed,
        algorithm_kws=dict(tolfun=5e-3),
    )
    start = time.time()
    if len(problems) == 1:
        cmaes.setup_problem(problem)
        cmaes.optimize()
        cmaes.save()
    else:
        optimizers = optimize.batch_optimize(cmaes, problems)
    print(
        f"CMAES for {len(problems)} subjects with {args.popsize} particles and {args.n_iter} iterations took a total walltime of {time.time()-start}s"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cmd", type=str, choices=["get_subs", "run"])
    # subject (and session)
    parser.add_argument("dataset", type=str, choices=["pnc", "imagen"])
    parser.add_argument("-fc_subs", type=lambda s: s.split(','), default=None)
    parser.add_argument(
        "-sc_sub", type=str, default="fc_sub", choices=["fc_sub", "micamics"]
    )
    parser.add_argument(
        "-fc_ses", type=str, default=None, help="imagen only", choices=["BL", "FU2"]
    )
    parser.add_argument(
        "-sc_level", type=str, default="sub", help="imagen only", choices=["sub", "ses"]
    )  # sub-> uses FU2 SC, ses-> uses session-specific SC
    parser.add_argument(
        "-sc_ses",
        type=str,
        default="FU2",
        help="imagen only; ignored if sc_level is provided",
        choices=["BL", "FU2"],
    )
    # simulation options
    parser.add_argument("-parc", type=str, default="schaefer-100")
    parser.add_argument("-exc_inter", type=int, default=1)
    parser.add_argument("-sc_config", type=str, default="mean001")
    parser.add_argument("-maps_name", type=str, default="6maps")
    # CMA-ES options
    parser.add_argument("-opt_seed", type=int, default=1)
    parser.add_argument("-popsize", type=int, default=210)
    parser.add_argument("-n_iter", type=int, default=81)
    parser.add_argument("-n_runs", type=int, default=2)
    # get_subs specific options
    parser.add_argument("-n_subs", type=int, default=0)
    parser.add_argument(
        "-fc_level",
        type=str,
        default="sub",
        help="pnc only",
        choices=["sub", "age_group"],
    )
    parser.add_argument("-subsample", action="store_true")
    parser.add_argument("-max_fic_trials", type=int, default=10)

    args = parser.parse_args()
    # post-process args
    args.exc_inter = bool(args.exc_inter)

    if args.cmd == "run":
        run(args)
    elif args.cmd == "get_subs":
        get_subs(args)
