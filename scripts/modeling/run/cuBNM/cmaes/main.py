import time
import numpy as np
from cuBNM import optimize
import os
import json
import argparse

PROJECT_DIR = os.environ["PROJECT_DIR"]
INPUT_DIR = os.environ["INPUT_DIR"]
OUTPUT_DIR = os.environ["OUTPUT_DIR"]


def get_paths(fc_sub, sc_sub, fc_ses, sc_ses, args):
    # specify maps, sc and out path
    maps_path = os.path.join(INPUT_DIR, f"{args.n_maps}maps_{args.parc}_zscore.txt")
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
    return fc_path, fcd_path, sc_path, out_path, maps_path


def get_subs(args):
    # get list of all subjects/groups
    if args.dataset == "pnc":
        if args.fc_level == "sub":
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
        fc_path, fcd_path, sc_path, out_path, maps_path = get_paths(
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
    fc_path, fcd_path, sc_path, out_path, maps_path = get_paths(
        args.fc_sub, args.sc_sub, args.fc_ses, args.sc_ses, args
    )
    # skip this CMA-ES run if it's already done
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
                        f"Skipping {args.fc_sub} seed {args.opt_seed} because it was already run"
                    )
                    return
    print(
        f"Running {args.fc_sub} seed {args.opt_seed} with:",
        "\n\tfc_path:",
        fc_path,
        "\n\tfcd_path:",
        fcd_path,
        "\n\tsc_path:",
        sc_path,
        "\n\tout_path:",
        out_path,
        "\n\tmaps_path:",
        maps_path,
    )

    emp_fc_tril = np.loadtxt(fc_path)
    emp_fcd_tril = np.loadtxt(fcd_path)

    start = time.time()
    problem = optimize.RWWProblem(
        params={
            "G": (0.5, 4.0),
            "wEE": (0.05, 0.75),
            "wEI": (0.05, 0.75),
        },
        het_params=["wEE", "wEI"],
        duration=450,
        TR=3,
        window_size=10,
        window_step=2,
        bw_params="heinzle2016-3T",
        exc_interhemispheric=args.exc_inter,
        sc_path=sc_path,
        emp_fc_tril=emp_fc_tril,
        emp_fcd_tril=emp_fcd_tril,
        maps_path=maps_path,
        out_dir=out_path,
    )
    cmaes = optimize.CMAESOptimizer(
        popsize=args.popsize,
        n_iter=args.n_iter,
        seed=args.opt_seed,
        algorithm_kws=dict(tolfun=5e-3),
    )
    cmaes.setup_problem(problem)
    cmaes.optimize()
    cmaes.save()
    print(
        f"CMAES with {args.popsize} particles and {args.n_iter} iterations took a total walltime of {time.time()-start}s"
    )
    return cmaes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("cmd", type=str, choices=["get_subs", "run"])
    # subject (and session)
    parser.add_argument("dataset", type=str, choices=["pnc", "imagen"])
    parser.add_argument(
        "-fc_sub", type=str, default=None, help="sub-<...> or group-<...>"
    )
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
    parser.add_argument("-n_maps", type=int, default=6)
    # CMA-ES options
    parser.add_argument("-opt_seed", type=int, default=1)
    parser.add_argument("-popsize", type=int, default=210)
    parser.add_argument("-n_iter", type=int, default=80)
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

    args = parser.parse_args()
    # post-process args
    args.exc_inter = bool(args.exc_inter)

    if args.cmd == "run":
        run(args)
    elif args.cmd == "get_subs":
        get_subs(args)
