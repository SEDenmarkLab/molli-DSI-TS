#!/bin/bash/env python

"""
This class should enable a user to input a series of tasks
that are automatically multi-tasked.

workflow:
1. user submits location for a list of XYZ files (e.g., for xTB/crest/gaussian/orca)
2. user submits instructions for job (input file structure/command line arguments)
3. user initiates process
4. script automatically distributes tasks across user specified # of nodes

inner details:
 
a. submitted XYZ files need to be located - simple error handling here would be useful
b. a default temp folder (user specified if desired) should be written for each input file
c. for each submitted job (determined by the number of procs accessed), a new temp directory should be created
d. xyz files are moved to the new temp directory
e. input files are created for each xyz file. These are determined by the program to be called
f. output file types are loaded
g. tasks are distributed
h. await tasks


###### EXAMPLE USAGE:

import generic_signal as gs
import get_pid as gp

gp.pid("PID")
test1 = gs.MultiTask(function="xtb", endswith=".mol2.xyz", numthreads = 12, cmd = "--gfnff")

test1.executor()


#############################

"""

import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import signal
import sys
import shutil
import uuid
import threading
import time
import glob

def signal_handler(signum, frame):
    stop_event.set()

class MultiTask:
    def __init__(self,  function: str= "", base_name: str= "", endswith: str="", numthreads: int=16, find_xyzs: bool = True, **kwargs):
        self.function = function #the program to be run
        self.base_name = base_name #the base (common) name for the xyz files
        self.endswith = endswith #the file extension/suffix (e.g., "_out.xyz")
        self.cmd = False
        self.numthreads = numthreads
        if "cwd" in kwargs: #user passed a custom working directory
            self.cwd = kwargs["cwd"]
        else:
            self.cwd = os.getcwd()
        if "cmd" in kwargs:
            self.cmd = kwargs["cmd"] #this used to be a boolean
        else:
            self.cmd = False
        
        #now find the XYZ files automatically
        if(find_xyzs):
            self.files = self.load_files()

        self.original_dir = os.path.dirname(self.files[0])
        self.procs = []
        self.stop_event = threading.Event()
        self.lock = threading.Lock()
        self.stop_processing = False
        


        #once the xyz files are initialized, the user should run a classmethod
    #init function
    def load_files(self):
        files = [os.path.join(self.cwd, f) for f in os.listdir('.') if f.endswith(f"{self.endswith}")]
        return files

    def check_files(self, output_files, temp_dir, temp_file_path):
        for src, dest_basename in output_files.items():
            print("src, dest_basename", src, dest_basename)
            src_path = os.path.join(temp_dir, src)
            print("src_path", src_path)
            
            if os.path.exists(src_path):
                print("\n\n\n")
                print(f"\n {self.base_name}")
                
                print("moving file", dest_basename)
                print("moving", self.original_dir)
                print("\n\n")
                shutil.move(src_path, os.path.join(self.original_dir, dest_basename)) ###########for some reason, dest_basename duplicates the name. Swap this with "src"
        print("moving self.files[0]", temp_file_path, self.files[0])
        
        #try to move back the original .xyz file (if it is still there)
        try:
            shutil.move(temp_file_path, self.files[0])
        except:
            print("warning: original input file is gone, blame xtb!")
        
        shutil.rmtree(temp_dir)
        
    ### ksp 5/28/2024
    #this function is meant to be run after the job is completed
    #it examines the temp directory for any produced files and transfers them all over
    #this avoids having to hard code the names of input files
    def get_output_files(self, temp_dir, file): #file is the basename (obtained from the .xyz file)
        filenames = [] #the file in the temp directory
        destname = [] # what we should call the file..
        print(f"./{temp_dir}/*")
        for outputfile in glob.glob(f"{temp_dir}/*"):
            filenames.append(os.path.basename(outputfile)) # what the files are actually called (output from the program)
            
            #sometimes a software outputs a filename that is not labeled with the input file, but instead a generic outputfile name that is the same across all jobs
            # e.g., input.inp produces: nbofile.47 rather than input_nbofile.47, which causes issues with conflicts
            # this situation needs to be dealt with
            # lets check the original input file and see if 
            # need to check to see if the original input filename (file) is in os.path.basename(file)
            print("here look here", os.path.basename(file))
            print("and here", os.path.basename(outputfile))
            
            if os.path.basename(file)[:-(len(self.endswith)+1)] not in os.path.basename(outputfile): # the filename is contained within the output file
                destname.append(f"{os.path.basename(file)}_{os.path.basename(outputfile)}")
                print(f"C1: the filename is called {os.path.basename(file)}_{os.path.basename(outputfile)}")
            elif os.path.basename(file) == os.path.basename(outputfile):
                destname.append(f"{os.path.basename(file)}")
                print(f"C2: the filename is called {os.path.basename(file)}")
            else:
                destname.append(f"{os.path.basename(outputfile)[:-(len(self.endswith)+1)]}")
                print("C3: the filename is called", os.path.basename(outputfile)[:-(len(self.endswith)+1)])
            
        print(dict(zip(filenames,destname)))
        return dict(zip(filenames, destname))
            



    #this function should be called iteratively...
    def process_file_xtb(self, file):
        
        #self.base_name = os.path.basename(temp_file_path)

        # Create a temporary directory for processing this file
        temp_dir = self.create_temp_dir()
        temp_file_path = shutil.move(file, temp_dir)  # Move the file to the temp directory
        
        # Construct the command - runs the whole directory? maybe we should call a specific xyz file
        if self.cmd:
            command = f"xtb {os.path.basename(temp_file_path)} {self.cmd}"
            print(f".. found custom command, sending '{command}' to xtb")
        else:
            command = f"xtb {os.path.basename(temp_file_path)} thermo --temp 298 --hess"
            print(f".. using default command, sending '{command}' to xtb")        
        # Execute the command in the temporary directory
        results = subprocess.run(command, shell=True, cwd=temp_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        

        # Define the files to move back and rename
        """output_files = {
            "xtb" : f"{os.path.basename(temp_file_path)}_xtb_out",
            "hessian": f"{os.path.basename(temp_file_path)}_hess",
            "g98.out": f"{os.path.basename(temp_file_path)}_g98.out",
            "wbo": f"{os.path.basename(temp_file_path)}_wbo",
            "charges": f"{os.path.basename(temp_file_path)}_charges",
            "vibspectrum": f"{os.path.basename(temp_file_path)}_vibspectrum",
        }"""

        #print("lets look at the submitted files")
        #print("output_files", output_files, temp_dir, temp_file_path)

        #self.get_output_files(temp_dir,file)
        #print("examine above \n\n\n")
        
        #self.check_files(output_files, temp_dir, temp_file_path)
        self.check_files(self.get_output_files(temp_dir, file), temp_dir, temp_file_path)
        #with open(os.path.basename(output_files['xtb']), "w+") as f:
        #    f.writelines(results.stdout)
        with open(f"{file}_xtb_out", "w+") as f:
            f.writelines(results.stdout)
        
        return f"Processed: {file}"
        
    def process_file_orca(self, file):
        
        #self.base_name = os.path.basename(temp_file_path)

        # Create a temporary directory for processing this file
        temp_dir = self.create_temp_dir()
        temp_file_path = shutil.move(file, temp_dir)  # Move the file to the temp directory
        
        # Construct the command - runs the whole directory? maybe we should call a specific xyz file
        if self.cmd:
            command = f"/home/_opt/orca/5.0.3/orca {os.path.basename(temp_file_path)} {self.cmd} > output"
            
            print(f".. found custom command, sending '{command}' to orca")
            print(" .. using /home/_opt/orca/5.0.3/orca for orca directory")

        else:
            command = f"/home/_opt/orca/5.0.3/orca {os.path.basename(temp_file_path)}"
            print(f".. using default command, sending '{command}' to orca") 
            print(" .. using /home/_opt/orca/5.0.3/orca for orca directory")
                   
        # Execute the command in the temporary directory
        results = subprocess.run(command, shell=True, cwd=temp_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        

        # Define the files to move back and rename
        """output_files = {
            "xtb" : f"{os.path.basename(temp_file_path)}_xtb_out",
            "hessian": f"{os.path.basename(temp_file_path)}_hess",
            "g98.out": f"{os.path.basename(temp_file_path)}_g98.out",
            "wbo": f"{os.path.basename(temp_file_path)}_wbo",
            "charges": f"{os.path.basename(temp_file_path)}_charges",
            "vibspectrum": f"{os.path.basename(temp_file_path)}_vibspectrum",
        }"""

        #print("lets look at the submitted files")
        #print("output_files", output_files, temp_dir, temp_file_path)

        #self.get_output_files(temp_dir,file)
        #print("examine above \n\n\n")
        
        #self.check_files(output_files, temp_dir, temp_file_path)
        self.check_files(self.get_output_files(temp_dir, file), temp_dir, temp_file_path)
        #with open(os.path.basename(output_files['xtb']), "w+") as f:
        #    f.writelines(results.stdout)
        with open(f"{file[:-(len(self.endswith)+1)]}_out", "w+") as f:
            f.writelines(results.stdout)
        
        return f"Processed: {file}"

    def signal_handler(self, signum, frame):
        with self.lock:
            self.stop_processing = True
        print("Signal received, stopping threads...")
        #self.stop_event.set()

    """     
    def signal_handler(self, signum, frame):
        print("Signal received, stopping threads...")
        self.stop_event.set()
    """
    
    def create_temp_dir(self):
        """Create a temporary directory in the current working directory."""
        current_working_directory = os.getcwd()
        temp_dir_name = os.path.join(current_working_directory, "temp_" + str(uuid.uuid4()))
        os.makedirs(temp_dir_name)
        return temp_dir_name

    def process_wrapper(self, function, file):
        with self.lock:
            if self.stop_processing:
                return "Processing stopped"
        return function(file)

    def threads(self, function):
        with ThreadPoolExecutor(max_workers=self.numthreads) as exec:
            #submit the function through process_wrapper
            futures = [exec.submit(self.process_wrapper,function, file) for file in self.files]
            for future in as_completed(futures):
                if self.stop_processing:
                    break
                print(future.result())
    """
    def threads(self, function):
        with ThreadPoolExecutor(max_workers=self.numthreads) as exec:
            futures = [exec.submit(function, file) for file in self.files]
        for future in as_completed(futures):
            print(future.result())
    """ 
    def executor(self):
        
        signal.signal(signal.SIGINT, self.signal_handler)
        signal.signal(signal.SIGTERM, self.signal_handler)
        if(self.function == "xtb"):
            self.threads(self.process_file_xtb)#process_file_xtb_quitable)
        if(self.function =="orca"):
            self.threads(self.process_file_orca)


"""
# List all .mol2.xyz files in the current directory
files = [os.path.join(os.getcwd(), f) for f in os.listdir('.') if f.endswith('.mol2.xyz')]

# Using ThreadPoolExecutor to run 8 jobs at a time
with ThreadPoolExecutor(max_workers=8) as executor:
    futures = [executor.submit(process_file, file) for file in files]
    
    for future in as_completed(futures):
        print(future.result())
"""
