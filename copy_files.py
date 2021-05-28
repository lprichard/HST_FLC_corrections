# Created this file 11 June 2020
# Last update Aug 28, 2020
# Written by Laura Prichard

import os
import glob
import shutil
import numpy as np

def copy_files_check(src_dir, dst_dir, files='*', rename=False, src_str='', dst_str='', print_cmd=False, move=False):
    """Given a source directory (src_dir) containing files with a specified format (files;
    default is '*' i.e. all files in the src_dir) and a destination directory to copy to (dst_dir), 
    this code will make the dst_dir if it doesn't exists and copy any specified files to the 
    dst_dir if they don't already exist. A destination file can be renamed by setting rename=True, 
    and specifying a string in the source file (src_str) to be replaced with a new string for the 
    destination file (dst_str). If it is slow to copy the files within the notebook, e.g. if running 
    over VPN connection, print the copy command to paste into a terminal with print_cmd=True.
    Needs os, glob, shutil. If move=True, files will be moved rather than copied between the input 
    directories and renamed if specified. Needs os, glob, shutil. Written by Laura Prichard April 2020."""
    
    # Move to source directory to get files
    os.chdir(src_dir)

    # Set mode to copy or move
    if move==False:
        mode = 'copy'
        md = 'cp'
    else:
        mode = '~move~'
        md = 'mv'
    
    # Check for files to copy
    # If a string is provided for the files, make a list with glob, else assume a list if provided 
    if isinstance(files, str):
        files_ext = files
        files = glob.glob(files)
    elif (type(files) is not list) and (type(files) is not np.ndarray):
        print('ERROR: `files` input must be a string (e.g. "*.fits"), list or array!!')
    else:
        files_ext=None
    
    print('=====================================================================')
    print('{} files to {} from {} to {}'.format(len(files), mode, src_dir, dst_dir))
    if rename==True: 
        print('Renaming {} to {} in files'.format(src_str, dst_str))
        rnm_str, rnmd_str=' and renaming ', ' and renamed '
    else: rnm_str=rnmd_str=''
    print('=====================================================================')
    
    # Check for destination directory and make if it doesn't exist
    if os.path.exists(dst_dir):
        print('Destination directory exists:', dst_dir)
    else:
        os.makedirs(dst_dir, 0o774)
        print('Created destination directory:', dst_dir)  
    
    # If printing the copy command, set the starting cmd variable
    if print_cmd==True: cmd=''
        
    # Loop over each file, create a source and destination, check if it exists, and copy if not 
    n=0
    for file in files:
        if os.path.isfile(file):
            # Define source and destination paths for each file
            src = os.path.join(src_dir, file)   #Source filepath of file to be copied

            # If a file is to be renamed with a replacement extension
            if rename==False:
                dst = os.path.join(dst_dir, file)   #Destination filepath of to be copied and renamed
            else: 
                # Replace the specified component with the new component
                dst = os.path.join(dst_dir, file.replace(src_str, dst_str))

            # Check if the file is already in the destination directory, if not it is copied
            if not os.path.exists(dst):
                # Either copy/move the files or add to the cmd variable
                if print_cmd==False:
                    if move==False:
                        print('Copying{} {} to {}'.format(rnm_str, src, dst))
                        shutil.copy(src, dst)
                    else:
                        print('~Moving~{} {} to {}'.format(rnm_str, src, dst))
                        shutil.move(src, dst)
                else:
                    if n>0: cmd+= ' && '
                    # Set command to copy/move the file to print for command line
                    cmd+= '{} {} {}'.format(md, src, dst)
                n+=1
            else:
                if move==False:
                    print('File exists: {}, not copying from {}'.format(dst, src_dir))
                else:
                    print('File exists: {}, not ~moving~ from {}'.format(dst, src_dir))
        else: print('WARNING: {} is not a file, skipping...'.format(file))

    # Print the command to copy the files or a summary of the copied files
    if print_cmd==True:
        # Check if all the specified files with the string extension provided need to be copied/moved
        if isinstance(files_ext, str) and ('*' in files_ext) and (len(files)==n) and (rename==False):
            all_files = os.path.join(src_dir, files_ext)
            print("Command to {} all {} '{}' files:".format(mode, n, files_ext))
            cmd = '{} {} {}'.format(md, all_files, dst_dir)
        else:
            # Else copy/move just the missing files
            print('Command to {} {} files:'.format(mode, n))
        
        print(cmd)
        return cmd
    else:
        if move==False:
            print('Copied{} {} files to {}'.format(rnmd_str, n, dst_dir))
        else: 
            print('~Moved~{} {} files to {}'.format(rnmd_str, n, dst_dir))
        
        
def rename(direc, files='*', src_str='', dst_str=''):
    """Given a directory (direc) containing files with a specified format (files;
    default is '*' i.e. all files in the direc), the code will rename the files
    by specifying a string in the source file (src_str) to be replaced with a new 
    string for the destination file (dst_str). This should run quickly even over VPN.
    Needs os, glob. Written by Laura Prichard May 2020."""
    
    # Move to input directroy
    os.chdir(direc)

    # Check for files to rename
    # If a string is provided for the files, make a list with glob, else assume a list if provided 
    if isinstance(files, str):
        files_ext = files
        files = glob.glob(files)
    elif (type(files) is not list) and (type(files) is not np.ndarray):
        print('ERROR: `files` input must be a string (e.g. "*.fits"), list or array!!')
    else:
        files_ext=None

    # List all specified files and rename 
    for f in files:
        os.rename(f, f.replace(src_str, dst_str))
    print('--------------------------------------------')
    print('Renamed {} files in {} from {} to {}'.format(len(glob.glob('*{}*'.format(dst_str))), direc, src_str, dst_str))
    print('--------------------------------------------')