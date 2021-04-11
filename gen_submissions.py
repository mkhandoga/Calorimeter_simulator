import os
from os import path

base_dir = "/afs/cern.ch/user/m/mykhando/Calo_ML/"
sub_scripts_dir = "sub_scripts/"
sub_files_dir = "sub_files/"

base_filename = "calosim_run"

output_dir = "/eos/atlas/atlascerngroupdisk/phys-sm/LowMu2017WZ/histograms_mykola/Calo_output2500k/"

def make_bash(base_dirname, sub_scripts_dir, filename, partition, n_clusters):
	if not os.path.exists(base_dirname+sub_scripts_dir):
		print ("Creating directory:", base_dirname+sub_scripts_dir)
		os.mkdir(base_dirname+sub_scripts_dir)
	print ("dirname", base_dirname, "filename", filename, "n_clusters", n_clusters)
	with open(base_dirname+sub_scripts_dir+filename,'w') as f:
		f.write(r'#!/bin/bash'+"\n")
		f.write(r'echo "Starting calorimeter simulation"'+"\n")
		f.write(f'dir={base_dirname}'+"\n")
		f.write(f'cd $dir'+"\n")
		f.write(r'echo "Home directory:" ${dir}'+"\n")
		f.write(r'echo "current directory:" `pwd`'+"\n")
		f.write(f'nclusters={n_clusters}'+"\n")
		f.write(f'job_partition=part{partition}'+"\n")
		f.write(f'outputdir={output_dir}'+"\n")
		f.write(r'echo "N clusters:" $nclusters'+'\n')
		f.write(r'echo "Job partition:" $job_partition'+'\n')
		f.write(r'echo "Output directory:" $outputdir'+'\n')
		complicated_str = """command=${dir}start_sim_calo.C"("$nclusters,'"'$job_partition'"','"'$outputdir'"'')'"""
		f.write(complicated_str+'\n')
		f.write(r'root -l -b -q "${command}"')
		f.close()



def make_subfile(base_dirname, sub_scripts_dir, sub_files_dir, filename):
	if not os.path.exists(base_dirname+sub_files_dir):
		print ("Creating directory:", base_dirname+sub_files_dir)
		os.mkdir(base_dirname+sub_files_dir)
	with open(base_dirname+sub_files_dir+filename+".sub",'w') as f:
		f.write(f'executable		= {base_dirname}{sub_scripts_dir}{filename}.sh'+"\n")
		multiline_str = """arguments       = $(ClusterId)$(ProcId)
universe = vanilla
+JobFlavour = "workday"
output          = output/calo_run.$(ClusterId).$(ProcId).out
error           = error/calo_run.$(ClusterId).$(ProcId).err
log             = log/calo_run.$(ClusterId).$(ProcId).log
requirements 	= (OpSysAndVer =?= "CentOS7")
max_retries 	= 5
queue"""
		f.write(multiline_str)
		f.close()

n_clusters = 2500000
clusters_per_file = 5000

n_submissions = n_clusters//clusters_per_file
tail = n_clusters%clusters_per_file


script_filenames = []
submission_filenames = []

print ("n_submissions",n_submissions)
print ("tail",tail)

#Generating steering scripts and submittion siles 
for s in range (n_submissions):
	print ("Creating bash scripts and sub files ,",s,"out of", n_submissions)
	filename = base_filename + "_part" + str(s) 
	script_filenames.append(filename)
	make_bash(base_dir, sub_scripts_dir, filename+ ".sh", s, clusters_per_file)
	submission_filenames.append(filename+".sub")
	make_subfile(base_dir,sub_scripts_dir,sub_files_dir,filename)

#taking care of the remainder
if (tail):
	print ("Creating bash script and sub file for the tail")
	filename = base_filename + "_part" + str(n_submissions+1) 
	script_filenames.append(filename+".sh")
	make_bash(base_dir, sub_scripts_dir, filename+ ".sh", n_submissions+1, tail)
	submission_filenames.append(filename+".sub")
	make_subfile(base_dir,sub_scripts_dir,sub_files_dir,filename)

#taking care that the auxilliary folders exist
aux_directories = ["output", "error","log"]

for d in aux_directories:
	if not os.path.exists(base_dir+d):
		os.mkdir(base_dir+d)

if not os.path.exists(output_dir):
	os.mkdir(output_dir)


#submitting jobs
for sub in submission_filenames:
	command = "condor_submit " + base_dir + sub_files_dir + sub
	print (command)
	os.system(command)



