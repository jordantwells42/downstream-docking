from pyrosetta import *
init()
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.select.movemap import *
import pandas as pd

enzs = ["iaaM", "McbB"]
focuses = [513, 46]
trps = ["5b", "5c", "6b", "6c", "7b", "7c"]
interface_strs = ["AB_C", "ABC_D"]
jump_nums = [2, 3]
scores = []
header = ["Enzyme", "TRP", "SC", "IFE", "ddG", "dG_dSASA_ratio", "dSASA", "num_buns"]
sf = create_score_function("ref2015_cart")

for enz, jump_num, interface_str in zip(enzs, jump_nums, interface_strs):
    for trp in trps:
        interface_analyzer = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(jump_num)
        interface_analyzer.set_scorefunction(sf)
        #interface_analyzer.set_pack_input(True)
        interface_analyzer.set_pack_separated(True)
        interface_analyzer.set_pack_rounds(100)
        interface_analyzer.set_interface(interface_str)

        pose = Pose()
        res_set = generate_nonstandard_residue_set(pose, params_list =  [f"./trps/{trp}/{trp}.params"])
        pose_from_file(pose, res_set, f"{enz}_{trp}.relax.pdb")

        interface_analyzer.apply(pose)



        data = interface_analyzer.get_all_data()
        SC = data.sc_value
        IFE = data.complexed_interface_score[1] - data.separated_interface_score[1]
        ddG = interface_analyzer.get_interface_dG()
        dG_dSASA_ratio = data.dG_dSASA_ratio 
        dSASA = data.dSASA_sc[1]/data.dSASA[1]
        num_buns = data.delta_unsat_hbonds

        scores.append([enz, trp, SC, IFE, ddG, dG_dSASA_ratio, dSASA, num_buns])

for enz in enzs:
    interface_analyzer = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover(jump_num)
    interface_analyzer.set_scorefunction(sf)
    #interface_analyzer.set_pack_input(True)
    interface_analyzer.set_pack_separated(True)
    interface_analyzer.set_pack_rounds(300)
    interface_analyzer.set_interface(interface_str)

    pose = pose_from_pdb(f"{enz}.relax.pdb")

    interface_analyzer.apply(pose)



    data = interface_analyzer.get_all_data()
    SC = data.sc_value
    IFE = data.complexed_interface_score[1] - data.separated_interface_score[1]
    ddG = interface_analyzer.get_interface_dG()
    dG_dSASA_ratio = data.dG_dSASA_ratio 
    dSASA = data.dSASA_sc[1]/data.dSASA[1]
    num_buns = data.delta_unsat_hbonds

    scores.append([enz, "wt trp", SC, IFE, ddG, dG_dSASA_ratio, dSASA, num_buns])


df = pd.DataFrame(scores, columns=header)
df.to_csv("scores.csv", index=False)
