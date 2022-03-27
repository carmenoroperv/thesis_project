import os
import shutil
import sys

samples = os.listdir("../../../DerivedData/main_data/")
samples_pts = [x for x in samples if x[0] == "C"]
print(len(samples_pts))
samples_pts = samples_pts[2:11]

#1) Create dir for each of the samples in our data folder
#2) Create a dir in our folder for each of the subfolders in the patient folder
#4) get a list of files with bam suffix in the folder
#5) copy the bam files 

path_dir = "../data/"
for pt in samples_pts:
    print(f"Copying sample {pt} data")
    sys.stdout.flush()
    if not os.path.exists(os.path.join(path_dir, pt)):
        os.mkdir(os.path.join(path_dir, pt))
    
    
    subfolds_pt = os.listdir("../../../DerivedData/main_data/" + str(pt) + "/")
    for subfold in subfolds_pt:
        subfold_path = "../data/" + str(pt) + "/"
        print(f"Copying subfolder {subfold} files")
        sys.stdout.flush()
        if not os.path.exists(os.path.join(subfold_path, subfold)):
            os.mkdir(os.path.join(subfold_path, subfold))
        
        files = os.listdir("../../../DerivedData/main_data/" + str(pt) + "/" + str(subfold) + "/")
        for file in files: 
            file_suf = os.path.splitext(file)[1]
            if file_suf == ".bam":
                print(f"Copying file {file}")
                sys.stdout.flush()
                file_to_copy = "../../../DerivedData/main_data/" + str(pt) + "/" + str(subfold) + "/" + str(file)
                to_path = "../data/" + str(pt) + "/" + str(subfold) + "/" + str(file)
                shutil.copy2(file_to_copy, to_path)
        
        
        
        
    

