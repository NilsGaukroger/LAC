# -*- coding: utf-8 -*-
"""Create the htc files and a .bat script to run simulations.

Each htc file that is made is/has:
    - a different wind speed
    - initial rotor speed that is determined based on wind speed and operation.dat file
    - 400 seconds total simulated, but first 200 seconds discarded
    - steady wind (no turbulence or step)
    - no shear
    - with tower shadow

Master htc file:
    - turbulent wind
    - placed in folder "htc_master/"

HAWC folder structure:
    - control/
    - data/
    - htc_master/
        - <master file name>.htc
    - (script creates) htc_steady/
        - (htc files to run)
    - (script creates) run_steady.bat

STEPS IN RUNNING THIS SCRIPT:
    1. Place your master htc file in a folder called "htc_master/". This folder should
       be in your hawc model folder, same level as control/ and data/.
    2. Copy this .py file and _run_hawc2_utils.py into the hawc folder, at the same
       level of control/ and data/.
    3. Update the input parameters in the script.
    4. Run this script: `python make_steady.py`. You should now have the htc files you
       want to run in htc_steady/. The .bat file will be in the same location as this
       script.
    5. Open one of the created htc files and skim it to make sure everything looks ok.
    6. Open a command prompt in the hawc folder. Run the .bat file by entering its name
       into the Command Prompt, then hitting Enter.
    7. When the simulations are done, check a few log files to make sure things ran okay.
    8. If you ran on multiple machines, transfer the results files to a single folder on
       a single computer with NumPy.
    9. Post-process the statistics using the post-processing script.
"""
import os
from _run_hawc2_utils import clean_directory, get_rotation_speed


operation_dat = './data/operation_DL.dat'  # path to operation.dat file
master_name = 'redesign_DL.htc'  # name of file in the htc_master folder
# wsps = [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]  # wind speeds you want to simulate
wsps = [5, 6, 7]
hawc2_exe = 'C:/hawc_simulations/hawc2_12.8_1900/HAWC2MB.exe' # path to HAWC2 executable

# =======================================================================================
# you shouldn't need to change anything below this line :)

clean_htc = False  # clean the htc directory? !!! WILL DELETE ALL FILES IN HTC DIR IF TRUE !!!
htc_dir = './htc_steady/'  # name of output directory !!! END WITH SLASH !!!
master_dir = './htc_master/'  # name of folder with htc master file !!! END WITH SLASH !!!
res_dir = './res_steady/'  # folder to save output !!! END WITH SLASH !!!
time_start = 200  # time to start recording output
time_stop = 400  # time to stop recording output

htc_master = master_dir + master_name  # path to htc master file
master_noext = os.path.splitext(master_name)[0]  # name of master file w/o extension

# clear or create the htc directory
clean_directory(htc_dir, clean_htc)

# load the master file's contents into memory
with open(htc_master, 'r') as htc:
    contents = htc.readlines()
    for line in contents:
        if line.lstrip().startswith('wind_ramp_abs'):
            print('WARNING!!! Your master htc has a wind_ramp_abs command. Are you '
                  + 'sure it is a turbulent simulation?')

# loop over the wind speeds
htc_files = []
for iw, wsp in enumerate(wsps):
    # define path names
    wsp_name = master_noext + '_' + ('%.1f' % (wsp)).zfill(4)  # name of htc file w/o extension
    htc_path = htc_dir + wsp_name + '.htc'

    # get initial rotation speed
    Omega0 = get_rotation_speed(wsp, operation_dat)

    # loop over the contents and update them accordingly
    output = 0
    for i, line in enumerate(contents):

        # if we've passed the dll block, toggle tracker to prevent new_htc_structure
        #   filenames from being overwritten
        if line.lstrip().startswith('end dll'):
            output = 1

        # simulation time
        if line.lstrip().startswith('time_stop'): 
            contents[i] = ('time_stop   %.1f ; \n' % time_stop)
        # log file
        elif line.lstrip().startswith('logfile'):
            contents[i] = ('logfile ./log/%s.log ; \n' % wsp_name)
        # initial rotation speed
        elif line.lstrip().startswith('mbdy2_ini_rotvec_d1'):
            contents[i] = ('mbdy2_ini_rotvec_d1 0.0 0.0 -1.0 %.2f ; \n' % Omega0)
        # wind speed
        elif line.lstrip().startswith('wsp'):
            contents[i] = ('wsp                     %.1f   ; \n' % wsp)
        # turbulence intensity
        elif line.lstrip().startswith('tint'):
            contents[i] = ('tint                    0.0   ; \n')
        # shear format
        elif line.lstrip().startswith('shear_format'):
            contents[i] = ('shear_format            1 %.1f ;  \n' % wsp)
        # turbulence format
        elif line.lstrip().startswith('turb_format'):
            contents[i] = ('turb_format             0     ;\n')
        # tower shadow
        elif line.lstrip().startswith('tower_shadow_method'):
            contents[i] = ('tower_shadow_method     3     ;\n')
        # output filename
        elif line.lstrip().startswith('filename') and output:
            contents[i] = ('filename %s%s ; \n' % (res_dir, wsp_name))
        # output time
        elif line.lstrip().startswith('time') and output:
            contents[i] = ('time %.1f %.1f ; \n' % (time_start, time_stop))

    # write the new htc file
    with open(htc_path, 'w') as htc:
        htc.writelines(contents)
    
    # append it's name to the list to be added to the bat file
    htc_files.append(htc_path)

# write the .bat file
with open('run_steady.bat', 'w') as bat:
    for htc_path in htc_files:
        bat.write('"%s" "%s" \n' % (hawc2_exe, htc_path))

