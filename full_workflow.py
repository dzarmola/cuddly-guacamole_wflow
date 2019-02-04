#!/usr/bin/env python2.7

""" Workflow for aligning protein family profiles.
Uses hhsuite for profile creation and pairwise alignments, cd-hit for representative selection,
mcl for clustering [optional], and clustalo for adding single sequences to the family alignments
[subject to change].

(C) Aleksandra Jarmolinska 2018-2019 a.jarmolinska@mimuw.edu.pl
"""

import sys,os,glob
import argparse,time
import re,subprocess
import shutil,tempfile
import scripts.extract_subfastas
import scripts.make_M_N
import scripts.extract_ev_for_clustering
import scripts.my_little_merger_error_cor
import scripts.get_representatives_from_clusters
import scripts.my_little_replacer



cwd=os.getcwd()
p2ch = "{}/scripts/cdhit/cdhit/".format(cwd)
p2hh = "${}/scripts/hhsuite/hhsuite/bin/".format(cwd)

def main(runname,data_folders,clustering=False,MPARAM="50",EV=1e-3,INF=1.8,NUM_REPR=3):
    """Each data_folder must contain org_fastas directory, optionally org_profiles directory,
    json with ranges to convert, structures with hhms to use for that"""


    dirname = "{}/{}".format(cwd, runname)
    os.makedirs(dirname)

    all_profiles = []
    all_fastas = []
    all_families = []

    for i,folder in enumerate(data_folders): # Each family is first processed separately to extract parts of profiles
        fdirname = "{}/group{}".format(dirname,i+1)
        os.makedirs(fdirname)
        profiles_dir = "{}/org_profiles".format(folder)
        fastas_dir = "{}/org_fastas".format(folder)
        struct_hhms_dir = "{}/struct_hhms".format(folder)
        json_path = "{}/ranges.json".format(folder)

        families = [re.sub("\.fa[sta]*$","",os.path.basename(f)) for f in glob.glob("{}/*.fa*".format(fastas_dir))]

        if not os.path.exists(profiles_dir) or \
                len(glob.glob("{}/*.hhm".format(profiles_dir)))!=len(glob.glob("{}/*.fa*".format(fastas_dir))):
            profiles_dir = "{}/org_profiles".format(fdirname)
            os.makedirs(profiles_dir)

            print "Making profiles, group",i+1

            for fasta in glob.glob("{}/*.fa*".format(fastas_dir)):
                new_name="{}/{}".format(profiles_dir,re.sub("\.fa[sta]*$",".hhm",fasta.split("/")[-1]))
                with open(os.devnull, "wb") as devnull:
                    subprocess.check_call(["hhmake","-M",MPARAM,"-i",fasta,"-o",new_name],stdout=devnull,stderr=subprocess.STDOUT)
            print "Done"
        else:
            old_profiles_dir = profiles_dir
            profiles_dir = "{}/org_profiles".format(fdirname)
            os.symlink(old_profiles_dir,profiles_dir)


        if os.path.exists(struct_hhms_dir) and os.path.exists(json_path):# we should do the profile shortening
            final_profiles_dir = "{}/final_profiles".format(fdirname)
            os.makedirs(final_profiles_dir)

            all_structs = "{}/all_structs.hhm".format(fdirname)
            with open(all_structs, 'wb') as outFile:
                for ifile in glob.glob("{}/*.hhm".format(struct_hhms_dir)):
                    with open(ifile, 'rb') as inFile:
                        shutil.copyfileobj(inFile, outFile)


            print "Looking for cores"
            for fam in families:
                print "Will find core for {}, group {}".format(fam, i + 1)
                with open(os.devnull, "wb") as devnull:
                    subprocess.check_call(["hhsearch","-i", "{}/{}.hhm".format(profiles_dir,fam),"-d",all_structs,
                                            "-o", "{}/{}.hhr".format(final_profiles_dir,fam)],
                                           stdout=devnull,stderr=subprocess.STDOUT)
                s,e = scripts.extract_subfastas.main("{}/{}.hhr".format(final_profiles_dir,fam),json_path)

                ###instead of extract_fastas.sh
                fd, path = tempfile.mkstemp()
                fd2, path2 = tempfile.mkstemp()
                fasta = glob.glob("{}/{}*.fa*".format(fastas_dir,fam))[0]
                with os.fdopen(fd, 'w') as f:
                    f.write(">moja\n{}\n".format(scripts.extract_subfastas.extract_seq_from_hhm("{}/{}.hhm".format(profiles_dir,fam),s,e)))
                with open(os.devnull, "wb") as devnull:
                    subprocess.check_call(["clustalo","--profile2",fasta,"--profile1",path,"-o",path2,"--force","--wrap=10000"],
                            stdout = devnull, stderr = subprocess.STDOUT)
                new_fasta = "{}/{}.fasta".format(final_profiles_dir,fam)
                with open(new_fasta,"w",0) as outFasta:
                    with os.fdopen(fd2,"r") as inFasta:
                        for _,line in enumerate(inFasta):
                            if not _: continue
                            if _==1:
                                s = line.index(line.strip("-"))
                                e = s + len(line.strip("-")) # TODO check
                            else:
                                if line[0]==">":
                                    outFasta.write(">{}{}\n".format(fam+"_" if fam not in line else "",  re.sub("\\[_0-9-]*$","",line[1:].strip()))) #line[1:].strip() ))
                                else:
                                    outFasta.write("{}\n".format(line[s:e+1]))
                os.remove(path)
                os.remove(path2)
                scripts.make_M_N.main(new_fasta,int(MPARAM))
                new_fasta_short = "{}/{}_{}.fasta".format(final_profiles_dir,fam,MPARAM)
                new_hhm_short = "{}/{}_{}.hhm".format(final_profiles_dir, fam, MPARAM)
                with open(os.devnull, "wb") as devnull:
                    subprocess.check_call(["hhmake","-M",MPARAM,"-i",new_fasta_short,"-name",fam,"-o",new_hhm_short],stdout=devnull,stderr=subprocess.STDOUT)
                #####end
        else: # we take profiles as they are
            final_profiles_dir = "{}/final_profiles".format(fdirname)
            os.makedirs(final_profiles_dir)
            for fam in families:
                fasta = glob.glob("{}/{}*.fa*".format(fastas_dir,fam))[0]
                new_fasta = "{}/{}.fasta".format(final_profiles_dir, fam)
                shutil.copyfile(fasta, new_fasta)
                scripts.make_M_N.main(new_fasta, int(MPARAM),family=fam)
                new_fasta_short = "{}/{}_{}.fasta".format(final_profiles_dir, fam, MPARAM)
                new_hhm_short = "{}/{}_{}.hhm".format(final_profiles_dir, fam, MPARAM)
                with open(os.devnull, "wb") as devnull:
                    subprocess.check_call(["hhmake","-M",MPARAM,"-i",new_fasta_short,"-name",fam,"-o",new_hhm_short],stdout=devnull,stderr=subprocess.STDOUT)
        all_profiles.append(final_profiles_dir)
        all_families.append(families)

    cluster_dir = "{}/clustering".format(dirname)
    os.makedirs(cluster_dir)
    all_hhm = "{}/all.hhm".format(cluster_dir)
    with open(all_hhm, "wb", 0) as outFile:
        for dir in all_profiles:
            for file in glob.glob("{}/*_{}.hhm".format(dir,MPARAM)):
                with open(file,"rb") as inFile:
                    shutil.copyfileobj(inFile, outFile)
    _hhrs = []
    for i,dir in enumerate(all_profiles):
        os.makedirs("{}/group{}".format(cluster_dir,i+1))
        for file in glob.glob("{}/*_{}.hhm".format(dir, MPARAM)):
            nfile = file.replace(".hhm", ".hhr").split("/")[-1]
            with open(os.devnull, "wb") as devnull:
                subprocess.check_call(["hhsearch", "-i", file, "-d", all_hhm,
                                       "-o", "{}/group{}/{}".format(cluster_dir,i+1, nfile)], stdout=devnull, stderr=subprocess.STDOUT)
            _hhrs.append("{}/group{}/{}".format(cluster_dir,i+1, nfile))

    cluster_file = ''

    if clustering: # currently just a stub - TODO: rewrite a bit mcl_in.sh
        scripts.extract_ev_for_clustering.parse_results(_hhrs,"{}/all_vs_all.out".format(cluster_dir),EV)
        #bash "$cwd"/scripts/mcl_in.sh all_vs_all_"$EV".out "$INF" clustering_"$part"_"$EV"_"$INF" > /dev/null  ##TODO
        cluster_file = ''

    save_file = "{}/save_{}.txt".format(dirname,EV)
    plot_file = "{}/plot_{}.png".format(dirname,EV)
    scripts.my_little_merger_error_cor.main(_hhrs,cluster_file,save_file,plot_file)

    representatives_dir = "{}/representatives".format(dirname)
    os.makedirs(representatives_dir)
    for i,group in enumerate(all_families):
        for fam in group:
            cdfasta = "{}/{}_{}_cdhit.fasta".format(all_profiles[i], fam, MPARAM)
            normfasta = "{}/{}_{}.fasta".format(all_profiles[i], fam, MPARAM)
            cdout = "{}/{}_{}_cdhit.out".format(all_profiles[i], fam, MPARAM)
            with open(os.devnull, "wb") as devnull:
                subprocess.check_call(["cd-hit", "-i", cdfasta, "-d", '0', "-sc", '1',"-c",'0.7',
                                       "-o", cdout], stdout=devnull, stderr=subprocess.STDOUT)
            scripts.get_representatives_from_clusters.main("{}.clstr".format(cdout), normfasta,
                                                           "{}/{}_representatives.fa".format(representatives_dir,fam),num_repr=NUM_REPR)
    scripts.my_little_replacer.main(save_file,representatives_dir,"{}/representative_alignment.fasta".format(dirname))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Workflow for alignment of family profiles [with extracting just part ofthe profile based on structure].
    Needs hhmake, hhsearch, clustalo, cd-hit and mcl (if clustering) installed and available on system path.""")
    parser.add_argument('directories', type=str, nargs='+', help="""Directories with groups of families. Each directory must contain 
    an 'org_fastas' directory with fasta formatted alignments (.fa or .fasta) of desired families. Additionally, if a 'ranges.json' files and 'struct_hhms'
    directory are presented all family profiles in this group will be aligned to profiles from struct_hhms, and just columns correspnding to ranges specified
    in the ranges.json will be used in further analysis.""")
    parser.add_argument('--run_name',"-n", type=str, default="", help='Name to be given to this run. Directory will be created in the current working dir.')
    parser.add_argument('--evalue',"-e", type=float, default=1e-3, help='E-value cutoff for significant hhsearch hits')
    parser.add_argument('--num_representatives',"-r", type=float, default=3, help='Number of representative sequences for each family')
    parser.add_argument('--mparam', "-m", type=str, default="50", help='Cutoff parameter for columns in hhmake')
    parser.add_argument('--inflation', "-i", type=float, default=1.8, help='Inflation value for mcl clustering')
    #parser.add_argument('--cluster', "-c", action='store_true', help='Should we add clusters to the plot? Currently ')

    args = parser.parse_args()


    mparam = args.mparam
    ev = args.evalue
    inf = args.inflation
    #clustering = args.cluster
    clustering = False
    num_repr = args.num_representatives

    runname = args.run_name
    if not runname:
        cwd = os.getcwd()
        stamp = time.strftime("+%H%M%S-%d%m%y")
        runname = "{}/run_{}".format(cwd, stamp)

    data_folders = args.directories

    main(runname, data_folders, clustering=clustering, MPARAM=mparam, EV=ev, INF=inf, NUM_REPR=num_repr)


