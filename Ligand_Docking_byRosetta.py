import os,sys,subprocess,shutil,wget,glob
from pybel import *

from docking_md.Preparation import *

if __name__ == "__main__":
    protein_name = sys.argv[1]
    protein_chain = sys.argv[2]
    protein_path = "./data/protein"
    ligand_path = "./data/ligand"
    output_path = "./data/output"

    p = Protein_prep("/home/yklee/yklee/ex_Tools/rosetta_3.13/main/")
    l = Ligand_prep("/home/yklee/yklee/ex_Tools/rosetta_3.13/main/")

    # Protein Preparation
    os.chdir(protein_path)
    p.clean_by_rosetta(protein_name,protein_chain)
    p.extract_origin_ligand(protein_name)

    os.chdir("../../")

    # Ligand Preparation
    os.chdir(ligand_path)
    ligand_files = sorted(glob.glob("*.smi"))
    for ligand in ligand_files:
        ligand_name = ligand.split(".")[0]
        os.mkdir(ligand_name)

        shutil.copy(ligand,ligand_name)
        shutil.copy("../protein/%s_lig.mol2"%protein_name,"./")
        os.chdir(ligand_name)
        with open(ligand,"r") as F:
            for line in F.readlines():
                tline = line.strip().split()
        smi = tline[0]
        ids = tline[1]
        l.Make_3D_by_pybel(smi,ids)
        l.Move_ligand_to_origin("../../../docking_md",ids,protein_name)
        l.Make_Conformers(ids)
        l.Make_Rosetta_Params("DWP","DWP",ids)

        os.chdir("../")

    os.chdir("../../")

    # Compound Preparation
    os.chdir(output_path)
    ligands = next(os.walk("../ligand/"))[1]

    for ligand in ligands:
        complex_line = []
        os.mkdir("%s_%s"%(protein_name,ligand))
        os.chdir("%s_%s"%(protein_name,ligand))

        shutil.copy("../../protein/%s_%s.pdb" % (protein_name, protein_chain), "./")
        shutil.copy("../../ligand/%s/DWP.params"%ligand,"./")
        shutil.copy("../../ligand/%s/DWP.pdb"%ligand,"./")
        shutil.copy("../../ligand/%s/DWP_conformers.pdb"%ligand,"./")

        with open("%s_%s.pdb"%(protein_name,protein_chain),"r") as F:
            for line in F.readlines():
                complex_line.append(line.strip())
        with open("DWP.pdb","r") as F:
            for line in F.readlines():
                complex_line.append(line.strip())
        with open("compound.pdb","w") as W:
            for line in complex_line:
                W.write(line+"\n")
        os.chdir("../")
    os.chdir("../../")

    # Docking Procedure
    os.chdir("./data/output")
    ligands = next(os.walk("."))[1]

    for ligand in ligands:
        os.chdir(ligand)
        with open("options.txt","w") as W:
            W.write("-in\n")
            W.write("\t-file\n")
            W.write("\t\t-s \'%s_%s.pdb DWP.pdb\'\n"%(protein_name,protein_chain))
            W.write("\t\t-extra_res_fa DWP.params\n")

            W.write("-packing\n")
            W.write("\t-ex1\n")
            W.write("\t-ex2\n")
            W.write("\t-no_optH false\n")
            W.write("\t-flip_HNQ true\n")
            W.write("\t-ignore_ligand_chi true\n")

            W.write("-parser\n")
            W.write("\t-protocol ../../../docking_md/dock.xml\n")

            W.write("-overwrite\n")

            W.write("-mistakes\n")
            W.write("\t-restore_pre_talaris_2013_behavior true\n")

        subprocess.call("/home/yklee/yklee/ex_Tools/rosetta_3.13/main/source/bin/rosetta_scripts.linuxgccrelease @ options.txt -nstruct 1000",shell=True)
        os.chdir("../")
