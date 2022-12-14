from pyrosetta import *
init("-mute all")
from alignment import *
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.select.movemap import *

class PreRelaxer:
    def __init__(self, scorefxn, num_repeats:int = 5):
        self.scorefxn = scorefxn
        self.num_repeats = num_repeats

    def __call__(self, pose):
        tf = core.pack.task.TaskFactory()
        tf.push_back(core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(core.pack.task.operation.IncludeCurrent())
        tf.push_back(core.pack.task.operation.NoRepackDisulfides())
        tf.push_back(core.pack.task.operation.RestrictToRepacking())

        fr = pyrosetta.rosetta.protocols.relax.FastRelax(
            scorefxn_in=self.scorefxn, standard_repeats=self.num_repeats
        )
        fr.cartesian(True)
        # As noted in Frenz 2020
        fr.constrain_relax_to_start_coords(True)
        fr.set_task_factory(tf)
        fr.max_iter(300)
        fr.apply(pose)

enzs = ["iaaM", "McbB"]
focuses = [513, 46]
trps = ["5b", "5c", "6b", "6c", "7b", "7c"]
jump_nums = [2, 3]


scorefxn = create_score_function("ref2015_cart")
relaxer = PreRelaxer(scorefxn)



for enz, focus in zip(enzs, focuses):
    pose = pose_from_pdb(f"./enzs/{enz}.pdb")
    
    for trp in trps:
        lig = Pose()
        print(f"./trps/{trp}/{trp}.params")
        # uise generate_nonstandard_residue_set
        res_set = generate_nonstandard_residue_set(lig, params_list =  [f"./trps/{trp}/{trp}.params"])
        pose_from_file(lig, res_set, f"./trps/{trp}/{trp}_0001.pdb")
        aid = ["CD2", "CZ2", "CZ3"]
        bid = ["C3", "C8", "C9"]

        # align the ligand to the enzyme
        wpose = Pose()
        wpose.assign(pose)
        w = wpose.residue(len(wpose))
        align_ligand_to_residue(lig, w, bid, aid)
        wpose.delete_residue_slow(len(wpose))
        wpose.append_residue_by_jump(lig.residue(1), focus, "", "", True)
        wpose.dump_pdb(f"{enz}_{trp}.pdb")

        # relax the ligand]
        print("Relaxing: " + enz + "_" + trp + "")
        relaxer(wpose)
        wpose.dump_pdb(f"{enz}_{trp}.relax.pdb")











