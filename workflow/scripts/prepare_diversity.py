import os
import shutil
import sys
import zipfile

sys.stderr = open(snakemake.log[0], "w")

zipped_files = []
i = 0
os.mkdir(str(snakemake.output))
while i < len(snakemake.input):
    file = str(snakemake.input[i])
    name_path = os.path.splitext(file)[0]
    name = os.path.split(name_path)[1]
    output_path = str(snakemake.output) + "/" + name + ".zip"
    shutil.copy2(file, output_path)
    zipped_files.append(output_path)
    i += 1

j = 0
while j < len(zipped_files):
    with zipfile.ZipFile(zipped_files[j], "r") as zip_ref:
        name = zipped_files[j].split("/")[-1]
        # print(name)
        new_dir = str(snakemake.output) + "/" + name
        # print(new_dir)
        zip_ref.extractall(os.path.splitext(new_dir)[0] + "/")
        os.remove(new_dir)
    j = j + 1

directory = os.listdir(str(snakemake.output))
# print(directory)

z = 0
while z < len(directory):
    subdir = os.listdir(str(snakemake.output) + "/" + directory[z])
    orig_dir = str(snakemake.output) + "/" + directory[z] + "/" + subdir[0]
    new_dir = str(snakemake.output) + "/" + directory[z]
    for f in os.listdir(orig_dir):
        path = orig_dir + "/" + f
        shutil.move(path, new_dir)
    os.rmdir(orig_dir)
    z = z + 1
