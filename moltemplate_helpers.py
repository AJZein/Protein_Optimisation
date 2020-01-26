"""Various functions for building a salt in water system using moltemplate. To be later adapted into a general purpose class to allow usage of moltemplate in python
(similar to plumed_file class)"""

import subprocess
import os
import numpy as np

def unique_enum(lower_bound, upper_bound, length=3):
    """Generator which produces all sets of numbers of a given length between defined bounds where no two sets have identical members.
    
    Parameters:
    ---
    **lower_bound: *int* **  
    The integer to start counting from
        
    **upper_bound: *int* **  
    The integer to count up to
        
    **length: *str, optional* **  
    The size of the sets to consider
        
    Yields:
    ---
    **current_index: *list* **   
    A list of integers of length `length`, whose values range from `lower_bound` to `upper_bound`.
    
    Examples
    ---
    >>> for i in unique_enumerator(1, 2, length=2):
        >>> print(i)
    [1, 1] [2, 1] [2, 2]
    
    >>> for i in unique_enumerator(1, 2, length=3):
        >>> print(i)
    [1, 1, 1] [2, 1, 1] [2, 2, 1] [2, 2, 2]

    """
    
    # Error checking
    if lower_bound > upper_bound:
        raise Exception("Lower bound exceeds upper bound")
    if not isinstance(lower_bound, int) or not isinstance(upper_bound, int):
        raise Exception("Lower and upper bound must be integers")
        
    # Initialising
    current_index = [lower_bound]*length
    upper_shape = [upper_bound]*length
    yield current_index
        
    # Starting to iterate
    while current_index != upper_shape: # End condition
        i = 0
        current_index[i] += 1
        while current_index[i] > upper_bound: # Overflow loop
            i += 1
            current_index[i] += 1
        for j in range(i): # Resetting other indices corrects to produce unique results
            current_index[j] = current_index[i]
           
        yield current_index

        
def cubicity(triplet):
    """Measures how 'cubic' a triplet of numbers is by checking how similar they are to one another by
    computing the sum of their distances from their average
    
    Parameters:
    ---
    **triplet: *list* **  
    A list containing three integers
        
    Returns:
    ---
    **output: *float* **   
    A measure of how 'cubic' the triplet is. Value ranges from 0 to inf where 0 is a perfect cube
    
    Notes:
    ---
    The results of this function dont have much significance alone, instead use it when comparing triplets. If cubicity(triplet1) < cubicity(triplet2) then triplet1 is more cubic
    """
    
    
    output = abs(sum(triplet)/3 - triplet[0]) \
    + abs(sum(triplet)/3 - triplet[2]) + abs(sum(triplet)/3 - triplet[2])
    return output


def decompose(number):
    """Computes a triplet of numbers whose product approximates the number provided and where no number in the triplet is 3 times greater than another.
    
    Parameters:
    ---
    **number: *int* **  
    The number to approximate
        
    Returns:
    ---
    **best_combo: *list* **   
    The first element of this list is a list of three numbers whose product approximates `number`. The second is error of the approximation.
    
    Notes:
    ---
    The problem this function addresses is how to build a cube with an arbritrary integer volume using integer side lengths. This is useful for creating a box of molecules, since the total number of molecules and the number of molecules along each axis must be an integer. However it is impossible to do this exactly unless the desired volume is a cubic number e.g. 9,27,64, etc. In the end, the more the cube is deformed into a cuboid the better the volume can be approximated, up to the extreme limit of using a cuboid with sides V*1*1 to create any volume V.
    
    Currently the function allows the side lengths to vary between 0.5*root to 1.7*root where root is the cube root of the volume. This choice balances approximation error and computational  time whithout distorting the cube too much but can be changed if desired. Allowing the sides to vary too much can drastically increase computational time. Computational time also increases drsatically as the volume increases, currently volumes up to 10 million can be approximated in a few seconds. Also for small volumes (<1000) a different set of side length conditions are used to increase accuracy (root +- 5).
    """
    
    # Setting the limits for the side lengths
    root = round(number**(1/3))
    if root <= 5: 
        lower = 0
    else:
        lower = min(round(root*0.5), root-5)
    upper = max(round(root*1.7), root+5)
    
    # With the lower and upper search bounds set, it is time to generate the triplets and check their errors
    min_error = abs(number-lower**3)
    min_factors = [lower]*3
    for factors in unique_enum(lower, upper, 3):
        error = abs(number-factors[0]*factors[1]*factors[2])
        if min_error < error:
            continue
        elif min_error > error:
            min_error = error
            min_factors = list(factors)
        elif cubicity(factors) < cubicity(min_factors): # Condition for resolving factors with equal error
            min_factors = list(factors)
            
    # Returns a list containing the triplet (in a list) and their error
    return [min_factors, min_error]


def modify_build(num_water, num_salt, grid=0, salt_trans=0, salt_sep=0, error=0, filename='system.lt'):
    """Writes a moltemplate file which can be used to build the system with the desired amount of water and salt
    
    Parameters:
    ---
    **num_water: *int* **  
    The number of water molecules the system will have
    
    **num_salt: *int* **  
    The number of salt molecules the system will have
    
    **grid: *float, optional* **  
    The side length of the simulation box. Function will attempt to use a reasonable value if not specified
    
    **salt_trans: *float, optional* **  
    How much to translate the salt molecules upon creation (so that they dont overlap with the water). Only required if the default values produces problems
    
    **salt_sep: *float, optional* **  
    How far apart to space the salt molecules. If not specified they will have the same spacing as the water molecules
    
    **error: *bool, optional* **  
    Prints out the build error when True (since the built system doesnt always have exactly the amount of water and salt specified)
    
    **filename: *str, optional* **  
    The name of the file to write out to
        
    Returns:
    ---
    **[num_water, num_salt]: *list* **   
    The number of water molecules and salt molecules the completed system will contain.
    
    Notes:
    ---
    The code is incapable of creating an exact amount of water and salt, rather is will attempt to approximate the values provided. For very small numbers (<11) it is exact.
    
    
    """
    
    # This part approximates the amount of water/salt desired with roughly cubic integers
    water_factors, water_error = decompose(num_water)
    salt_factors, salt_error = decompose(num_salt)
    num_water = np.prod(water_factors)
    num_salt = np.prod(salt_factors)
    if error:
        print("Build error of {} molecules".format(water_error + salt_error))
    
    # Setting grid size if not specified using an arbitrary reciprocal density (easier to work with)
    if not grid:
        density = 25 # Cubic Angstroms per Water Molecule
        volume = density*num_water
        grid = round(volume**(1/3), 2)
       
    # Setting the spacings so that the water/salt lattice spans the box 
    water_spacing, salt_spacing = [], []
    for i in range(3):
        water_spacing.append(round(grid/water_factors[i], 2))
        salt_spacing.append(round(grid/salt_factors[i], 2))
    
    if not salt_trans: # Setting salt translation to avoid overlap with water, may require tinkering
        salt_trans = water_spacing[0]/2
        
    if salt_sep: # Setting how far apart Na and Cl will be
        salt_sep = round(salt_sep/np.sqrt(3), 2)
    if not salt_sep:
        salt_sep = water_spacing[0]
    
        
    na_move = round(salt_trans, 2)
    cl_move = round(salt_trans + salt_sep, 2)
    
    # Creating lines for writing into build file
    
    grid_lines = [f"0 {grid} xlo xhi\n"]
    grid_lines.append(f"0 {grid} ylo yhi\n")
    grid_lines.append(f"0 {grid} zlo zhi\n")
    
    wat_lines = [f"[{water_factors[0]}].move(0.00, 0.00, {water_spacing[0]})\n"]
    wat_lines.append(f"[{water_factors[1]}].move(0.00, {water_spacing[1]}, 0.01)\n")
    wat_lines.append(f"[{water_factors[2]}].move({water_spacing[2]}, 0.01, 0.01)\n")
    
    na_lines = [f"[{salt_factors[0]}].move(0.00, 0.00, {salt_spacing[0]})\n"]
    na_lines.append(f"[{salt_factors[1]}].move(0.00, {salt_spacing[1]}, 0.01)\n")
    na_lines.append(f"[{salt_factors[2]}].move({salt_spacing[2]}, 0.01, 0.01)\n")
    
    cl_lines = [f"[{salt_factors[0]}].move(0.00, 0.00, {salt_spacing[0]})\n"]
    cl_lines.append(f"[{salt_factors[1]}].move(0.00, {salt_spacing[1]}, 0.01)\n")
    cl_lines.append(f"[{salt_factors[2]}].move({salt_spacing[2]}, 0.01, 0.01)\n")
    
    # Writing into build file
    
    with open(filename, "r") as file:
        data = file.readlines()
        
    for i in range(len(data)):
        line = data[i]
        if "write_once" in line:
            for j in range(1,4):
                data[i+j] = grid_lines[j-1]

        if "new SPCE" in line:
            for j in range(1, 4):
                data[i+j] = wat_lines[j-1]

        if "new NaIon" in line:
            for j in range(1, 4):
                data[i+j] = na_lines[j-1]

        if "new ClIon" in line:
            for j in range(1, 4):
                data[i+j] = cl_lines[j-1]

        if "na[*]" in line:
            data[i] = f"na[*][*][*].move({na_move},{na_move},{na_move})\n"
            data[i+1] = f"cl[*][*][*].move({cl_move},{cl_move},{cl_move})\n"
                
    with open(filename, "w") as file:
        file.writelines(data)
        
    return [num_water, num_salt]
        

    
    
def build_system(filename="system.lt"):
    """Simple function to run a moltemplate script and delete all outputs except the lammps data file
    
    Parameters:
    ---
    **filename: *str, optional* **  
    The name of the moltemplate file
    
    """
    
    
    root = filename[:-3]
    bashCommand = f"moltemplate.sh {filename}"
    normal = subprocess.run(bashCommand.split(),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=True,
        text=True)

    print(normal.stdout)

    os.remove(f"{root}.in.init")
    os.remove(f"{root}.in.settings")
    os.remove(f"{root}.in")
