import os,sys,subprocess,shutil,wget,glob
from pybel import *

from docking_md.Preparation import *
# 0317 - modified ligand preparation
# 0519 - modified to double chain docking
def Make_compound(pn,pc,li):
    complex_line = []
    shutil.copy("../../protein/%s_%s.pdb" % (pn, pc), "./")
    shutil.copy("../../ligand/%s/DWP.params" % li, "./")
    shutil.copy("../../ligand/%s/DWP.pdb" % li, "./")
    shutil.copy("../../ligand/%s/DWP_conformers.pdb" % li, "./")

    with open("%s_%s.pdb" % (pn, pc), "r") as F:
        for line in F.readlines():
            complex_line.append(line.strip())
    with open("DWP.pdb", "r") as F:
        for line in F.readlines():
            complex_line.append(line.strip())
    with open("compound.pdb", "w") as W:
        for line in complex_line:
            W.write(line + "\n")

def Make_rosetta_flags(pn,pc):
    with open("options.txt", "w") as W:
        W.write("-in\n")
        W.write("\t-file\n")
        W.write("\t\t-s \'%s_%s.pdb DWP.pdb\'\n" % (pn, pc))
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

def concatering_proteins(pn,pc):
    temp = []
    flist = sorted(glob.glob("./*.rmv.pdb"))
    for f in flist:
        with open(f,"r") as F:
            for line in F.readlines():
                temp.append(line.strip())
    with open("%s_%s.pdb"%(pn,pc),"w") as W:
        for line in temp:
            W.write(line+"\n")
				
def extract_origin_ligand(pn):
    temp = []
    with open("%s.pdb"%pn,"r") as F:
        for line in F.readlines():
            if line.startswith("HETATM") :
                temp.append(line.strip())
            else:
                pass
    with open("%s_lig.pdb"%pn,"w") as W:
        for line in temp:
            W.write(line+"\n")
    subprocess.call("obabel %s_lig.pdb -O %s_lig.mol2"%(pn,pn),shell=True)
if __name__ == "__main__":
    protein_name = sys.argv[1]
    protein_chain = sys.argv[2]
    protein_path = "./data/protein"
    ligand_path = "./data/ligand"
    output_path = "./data/output"

    os.chdir(protein_path)
    if len(protein_chain)>1:
        for chain in protein_chain:
            p = Protein_prep("/home/yklee/yklee/ex_Tools/rosetta_3.13/main/",protein_name,chain)
            #os.chdir(protein_path)
            p.clean_by_rosetta()
            os.system("mv %s_%s.pdb %s_%s.rmv.pdb"%(protein_name,chain,protein_name,chain))
        concatering_proteins(protein_name,protein_chain)
        extract_origin_ligand(protein_name)
    else:
        p = Protein_prep("/home/yklee/yklee/ex_Tools/rosetta_3.13/main/",protein_name,protein_chain)
        p.clean_by_rosetta()
        p.extract_origin_ligand()

    # Protein Preparation
    """os.chdir(protein_path)
    p.clean_by_rosetta()
    p.extract_origin_ligand()"""

    os.chdir("../../")

    # Ligand Preparation
    l = Ligand_prep("/home/yklee/yklee/ex_Tools/rosetta_3.13/main/")
    os.chdir(ligand_path)
    ligand_files = sorted(glob.glob("*.smi"))
    for ligand in ligand_files:
        q = 0
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
        #try:
        l.Move_ligand_to_origin("../../../docking_md",ids,protein_name)
        #except:
        #    os.chdir("../")
        #    shutil.rmtree(ligand_name)
        #    q = 1
        if q == 1:
            pass
        else:
            l.Make_Conformers(ids)
            tmp = l.Make_Rosetta_Params("DWP","DWP",ids)
            if tmp == 1:
                os.chdir("../")
            else:
                os.chdir("../")
                shutil.rmtree(ligand_name)


    os.chdir("../../")

    # Compound Preparation
    os.chdir(output_path)
    ligands = next(os.walk("../ligand/"))[1]
    a = set(["DWP.params","DWP.pdb","DWP_conformers.pdb"])
    for ligand in ligands:
        ligand_path = os.path.join("../ligand",ligand)
        q = 0
        if os.path.exists(ligand_path):
            aa = set(next(os.walk(ligand_path))[2])
        else:
            q = 1
        os.mkdir("%s_%s"%(protein_name,ligand))
        os.chdir("%s_%s"%(protein_name,ligand))
        if q == 0 :
            if a.issubset(aa):
                Make_compound(protein_name,protein_chain,ligand)
                os.chdir("../")
            else:
                os.chdir("../")
                shutil.rmtree("%s_%s"%(protein_name,ligand))
        else:
            os.chdir("../")
            shutil.rmtree("%s_%s"%(protein_name,ligand))
    os.chdir("../../")

    # Docking Procedure
    os.chdir("./data/output")
    ligands = next(os.walk("."))[1]

    for ligand in ligands:
        os.chdir(ligand)
        Make_rosetta_flags(protein_name,protein_chain)
        try:
            subprocess.call("/home/yklee/yklee/ex_Tools/rosetta_3.13/main/source/bin/rosetta_scripts.linuxgccrelease @ options.txt -nstruct 10",shell=True)
        except:
            pass
        os.chdir("../")
        fs = next(os.walk(ligand))[2]
        if len(fs) == 1007:
            pass
        else:
            with open("./NoOut_ID.txt","a") as F:
                F.write(ligand + "\n")
