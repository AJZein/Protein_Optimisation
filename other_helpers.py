# Additional non specific useful functions
import subprocess
import shutil
import os
import numpy as np

def run_parallel(filename, args='', num_cores=1, out=False):
    """Simple function to run a python script in serial or parallel and provide it with some parameters
    
    Parameters:
    ---
    **filename: *str* **  
    The name of the python file
    
    **args: *str* **  
    The parameters to pass to the python script, seperated by spaces
    
    **num_cores: *int* **  
    Number of cores/threads to use when running the script. 1 runs the script in parallel
    
    **out: *bool* **  
    If True, prints the output of the script
    
    """
    
    
    if num_cores == 1:
        bashCommand = f"python {filename} {args}"
    else:
        bashCommand = f"mpirun -np {num_cores} python {filename} {args}"
    normal = subprocess.run(bashCommand.split(),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True,
        text=True)

    if out:
        print(normal.stdout)
        
def read_atoms(filename):
    """Gets the entire 'Atoms' section from a lammps data file and puts it into a numpy array. 
    
    Parameters:
    ---
    **filename: *str* **  
    The name of the lammps data file
    
    Returns:
    ---
    **all_atoms: *str* **
    A numpy array containing the entire atoms section of the data file
    
    Notes:
    ---
    The 'Atoms' section of a lammps data file is formatted as such: [atom_num, molecule-tag, atom-type, q, x, y, z, nx, ny, nz]. nx, ny, nz optional
    
    """
    with open(filename, 'r') as file:
        data = file.readlines()
    
    for i in range(len(data)):
        line = data[i]
        if "atoms" in line:
            num_atoms = int(line.strip().split()[0])
            break
            
    all_atoms = []
    for i in range(len(data)): # Iterating through file to find where atoms start
        line = data[i]
        if "atoms" == line.strip().lower()[:5]:
            for j in range(2, num_atoms+2): # Iterating through atoms section to get all atoms
                line = data[i+j].split('#')[0].split()
                line = [float(x) for x in line]
                all_atoms.append(line)
            
            break
    
    all_atoms = np.array(all_atoms)
    return all_atoms

def clean_up(original, new, identifier='', storage_path='', mode='delete'):
    """Cleans up current directory by moving generated files out into another directory (or deleting them)
    
    This functions takes two lists of filenames, and filenames in `new` that arent in `original` will be removed/copied. Intended to be used by fetching all filenames before and after running some program, alowing generated files to be dealt with.
    
    Parameters:
    ---
    **original: *list* **  
    A list of strings containing a set of filenames to retain
    
    **new: *list* **  
    A list of strings containing a set of filenames. Any filenames not in original will be removed
    
    **identifier: *str, optional* **  
    When moving/copying files, this string is added as a prefix to their filenames
    
    **storage_path: *str, optional* **  
    Where to move/copy the files to
    
    **mode: *str, optional* **  
    Determines which actions the function takes. In 'delete' mode the files are deleted, in 'copy' or 'move' the files are copied or moved to the storage_path directory
    """
    for filename in new:
        if filename not in original:
            if mode == 'move':
                shutil.move(f"./{filename}", f"{storage_path}/{identifier}{filename}")
            elif mode=='delete':
                os.remove(f"./{filename}")
            elif mode=='copy':
                shutil.copyfile(f"./{filename}", f"{storage_path}/{identifier}{filename}")
                
def time_to_steps(time, step_size=2):
    """Given some step size, computes how many steps are in the given time
    
    Parameters:
    ---
    **time: *str* **  
    The time in the format 'xps' or 'xns' where x is some integer
    
    **step_size: *int* **  
    The step size in femtoseconds.
    
    Returns:
    ---
    **steps: *int* **  
    The number of steps in the provided time
    
    Raises:
    ---
    Format Exception
    If the time provided isn't in picoseconds or nanoseconds
    """
    
    time = time.strip()
    if 'ps' in time:
        steps = int(time.split('ps')[0].split()[0])*1000/step_size
    elif 'ns' in time:
        steps = int(time.split('ns')[0].split()[0])*1000000/step_size
    elif time == '0':
        steps = int(time)
    else:
        raise Exception("Incorrect format, specify an integer immediately followed by either 'ps' or 'ns'")
        
    return steps