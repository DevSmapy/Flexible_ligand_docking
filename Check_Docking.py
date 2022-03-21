import shutil,os,sys

if __name__ == "__main__":
    d_list = next(os.walk("./data/output/"))[1]

    for d in d_list:
        p = os.path.join("./data/output/",d)
        fs = len(next(os.walk(p))[2])
        if fs == 1007:
            with open("Ok.list","a") as F:
                F.write(d + "\n")
        else:
            with open("no_out.list","a") as F:
                F.write(d + "\n")

