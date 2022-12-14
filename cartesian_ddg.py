import argparse
from pathlib import Path
from os import chdir
from enum import Enum
from typing import Tuple, List
from abc import ABC, abstractmethod
import csv

import pandas as pd

from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.select.movemap import *


class PreRelaxer:
    def __init__(self, scorefxn, num_repeats:int = 20):
        self.scorefxn = scorefxn
        self.num_repeats = num_repeats

    def __call__(self, pose):
        tf = core.pack.task.TaskFactory()
        tf.push_back(core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(core.pack.task.operation.IncludeCurrent())
        # tf.push_back(core.pack.task.operation.NoRepackDisulfides())
        tf.push_back(core.pack.task.operation.RestrictToRepacking())

        fr = pyrosetta.rosetta.protocols.relax.FastRelax(
            scorefxn_in=self.scorefxn, standard_repeats=self.num_repeats
        )
        fr.cartesian(True)
        # As noted in Frenz 2020
        fr.constrain_relax_to_start_coords(False)
        fr.set_task_factory(tf)
        fr.max_iter(200)
        fr.apply(pose)



# pyrosetta.rosetta.protocols.ddg.CartesianddG.run(pose_from_pdb("./example/pdbs/3nir.relax.pdb"))

class CartesianLocalRelaxer:
    def __init__(self):
        pass

    def __call__(self, pose):
        pyrosetta.rosetta.protocols.ddg.CartesianddG.run(pose)



class InputHandler:
    def read_csv(self, mutation_csv: Path) -> Tuple[str, str, List[str]]:
        with open(mutation_csv, "r") as f:
            reader = csv.reader(f)
            next(reader, None)

            for i, row in enumerate(reader):
                pdb = row[0]
                chain_id = row[1]
                mutations = row[2].split(" ")
                yield (pdb.lower() + ".pdb", chain_id, mutations)


def create_mut_file(name: str, pose, chain_id, mutations):
    output = ""
    total = len(mutations)
    output += f"total {total}\n"

    for mutation in mutations:
        old_aa = mutation[0]
        new_aa = mutation[-1]
        pdb_position = int(mutation[1:-1])



        position = pose.pdb_info().pdb2pose(chain_id, pdb_position)
        
        if pose.residue(position).name1() != old_aa:
            print(f"Skipping mutation {mutation} for {name} due to WT mismatch")
            continue   

        output += f"1\n{old_aa} {position} {new_aa}\n"

    with open(name, "w") as f:
        f.write(output)


def cartesian_mutate(
    pdb_dir: Path,
    mutation_csv: Path,
    output_csv: Path,
    num_iterations: int,
    prerelaxer: PreRelaxer,
    localrelaxer: CartesianLocalRelaxer,
    dump_pdb: bool,
    overwrite: bool = True,
):
    assert isinstance(pdb_dir, Path)
    assert isinstance(mutation_csv, Path)
    assert isinstance(output_csv, Path)

    assert pdb_dir.is_dir()
    assert mutation_csv.is_file() and mutation_csv.suffix == ".csv", mutation_csv.resolve()

    pdb_dir = pdb_dir.resolve()
    mutation_csv = mutation_csv.resolve()
    output_csv = output_csv.resolve()

    working_dir = pdb_dir.parent / "working/"
    ddg_out_dir = pdb_dir.parent / "cartesian_ddg"

    working_dir.mkdir(0o774, parents=True, exist_ok=True)
    ddg_out_dir.mkdir(0o774, parents=True, exist_ok=True)

    input_handler = InputHandler()

    mut_rows = list(input_handler.read_csv(mutation_csv))

    with open("flags.txt", "r") as f:
        base_flags = (
            f"-ddg::dump_pdbs {str(dump_pdb).lower()} -ddg:iterations {num_iterations}"
            + " " + " ".join(f.read().split("\n"))
        )

    if "pymp" not in dir():
        import pymp

    with pymp.Parallel() as p:
        for idx in p.xrange(len(mut_rows)):
            pdb, chain_id, mutations = mut_rows[idx]

            out_dir = ddg_out_dir / pdb.split('.')[0]
            out_dir.mkdir(0o774, parents=True, exist_ok=True)

            chdir(out_dir)

            pdb = (pdb_dir / pdb).with_suffix(".pdb")

            assert pdb.is_file(), pdb.resolve()

            cleaned_pdb = pdb.with_suffix(".clean.pdb")
            relaxed_pdb = pdb.with_suffix(".relax.pdb")
            mut_file = working_dir / f"{pdb.stem}_mut_file.txt"

            ddg_mut_out = (Path().cwd() / f"{mut_file.stem}.ddg").resolve()
            pdb_ddg_csv =  (Path().cwd() / f"{pdb.stem}_mut_ddg.csv").resolve()

            if ddg_mut_out.is_file():
                if ddg_mut_out.stat().st_size > 0:
                    if overwrite:
                        ddg_mut_out.unlink()
                    else:
                        continue

                elif ddg_mut_out.stat().st_size == 0:
                    ddg_mut_out.unlink()

            

            # This is where cartesian_ddG runs
            p.print(f"Thread: {p.thread_num:3d} - Running Rosetta Procols for {pdb.name}...")
            

            
            
            try:
                # Check for cleaned (only ATOM and TER lines) file
                if not cleaned_pdb.is_file():
                    p.print(f"Thread: {p.thread_num:3d} - Cleaning PDB: {relaxed_pdb}")
                    pyrosetta.toolbox.cleaning.cleanATOM(str(pdb))
                # Check for cached pre-relaxed pdb
                if relaxed_pdb.is_file():
                    p.print(f"Thread: {p.thread_num:3d} - Using cached pdb: {relaxed_pdb}")
                    pose = pose_from_pdb(str(relaxed_pdb))

                else:
                    p.print(f"Thread: {p.thread_num:3d} - Pre-relaxing pdb: {cleaned_pdb}")
                    pose = pose_from_pdb(str(cleaned_pdb))
                    prerelaxer(pose)
                    pose.dump_pdb(str(relaxed_pdb))
                
                

                create_mut_file(mut_file, pose, chain_id, mutations)
                all_flags = (
                    f"-in:file:s {str(relaxed_pdb)} -ddg::mut_file {str(mut_file)}"
                    + " "
                    + base_flags
                )

            except Exception as e:
                p.print(f"Thread: {p.thread_num:3d} - Exited with Exception during pose relaxation:\n{e}")
                continue

            try:
                # Run Cartesian ddG
                p.print(f"Thread: {p.thread_num:3d} - Running Cartesian ddG on {pdb.name}...")
                #init(all_flags)
                #localrelaxer(pose)

            except Exception as e:
                p.print(f"Thread: {p.thread_num:3d} - Exited with Exception running cartesian ddg:\n{e}")
                continue

            else:
                pass

    print(f"Wrote out Cartesian ddG calculations to: {output_csv.resolve()}")



def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-p", "--pdb-dir", help="Directory containing all input pdbs", type=Path
    )
    parser.add_argument(
        "-m",
        "--mutation-csv",
        help="A csv containing the point mutations to be made (specified in README.md)",
        type=Path,
    )
    parser.add_argument(
        "-o",
        "--output-csv",
        help="What to name the output csv, which contains pdb ids, dGs and ddGs",
        type=Path,
    )

    parser.add_argument(
        "-n",
        "--repeats",
        help="How many times to repeat the protocol on a mutation site",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--dump-pdb", help="Dump Mutated PDBs to the pdb_dir", action="store_true"
    )

    # Init Rosetta
    print("\nInitializing Rosetta")
    init("-mute all -missing_density_to_jump -fa_max_dis 9.0")
    args = parser.parse_args()

    assert args.pdb_dir.is_dir()
    assert args.mutation_csv.exists()

    scorefxn = create_score_function("ref2015_cart")

    prerelaxer = PreRelaxer(scorefxn)

    localrelaxer = CartesianLocalRelaxer()
    
    cartesian_mutate(
        args.pdb_dir,
        args.mutation_csv,
        args.output_csv,
        args.repeats,
        prerelaxer,
        localrelaxer,
        args.dump_pdb,
    )

if __name__ == "__main__":
    import pymp

    main()