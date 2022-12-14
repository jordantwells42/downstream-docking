from pyrosetta import *
from pyrosetta.toolbox.cleaning import cleanATOM

from pyrosetta import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.select.movemap import *

init("-mute all")
# create an array of all pdbs in this directory, harcoded not glob

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

names = ["McbB"]
scorefxn = create_score_function("ref2015_cart")
jump_nums = [2, 3]

relaxer = PreRelaxer(scorefxn)
for name in names:
    print("Relaxin: " + name)
    pose = pose_from_pdb(f"./enzs/{name}.pdb")
    relaxer(pose)
    pose.dump_pdb("{}.relax.pdb".format(name))



