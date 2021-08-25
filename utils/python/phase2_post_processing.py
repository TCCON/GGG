#!/usr/bin/env python3
from argparse import ArgumentParser
from datetime import datetime
from glob import glob
import os
from pathlib import Path
import re
import shutil
from subprocess import call
import sys

now = datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
windows_to_remove = {'co_4233', 'wco2_6500', 'hcl_5577', 'hcl_5597', 'hcl_5683', 'hcl_5706', 'hcl_5719', 
                     'hcl_5749', 'hcl_5754', 'hcl_5763', 'hcl_5767', 'hcl_5779', 'hcl_5790'}


def backup(filename):
    file_path = Path(filename)
    backup_path = file_path.parent / '{}.bak.{}'.format(file_path.name, now)
    shutil.move(file_path, backup_path)
    return backup_path


def find_one_file(pattern):
    files = glob(str(pattern))
    if len(files) == 1:
        return files[0]
    else:
        raise IOError('Expected one file matching {}, found {}'.format(pattern, len(files)))


def get_multiggg_line_window(line):
    m = re.search(r'gfit (\w+)\.', line)
    return m.group(1)


def update_multiggg(multiggg_path, force_update=False):
    if not force_update and not multiggg_needs_updated(multiggg_path):
        return

    now = datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
    backup_path = backup(multiggg_path)

    with open(multiggg_path, 'w') as wh, open(backup_path, 'r') as rh:
        for line in rh:
            window = get_multiggg_line_window(line)
            if window in windows_to_remove:
                wh.write(':{}'.format(line))
            else:
                wh.write(line)


def multiggg_needs_updated(multiggg_path):
    with open(multiggg_path) as rh:
        for line in rh:
            window = get_multiggg_line_window(line)
            if window in windows_to_remove and not line.startswith(':'):
                return True

    return False


def read_windows_file(windows_file):
    sfs = dict()
    with open(windows_file) as rh:
        nhead = int(rh.readline().split()[0])
        for _ in range(nhead-1):
            rh.readline()

        for line in rh:
            if line.startswith(':'):
                continue

            # Get the center the same way gsetup does
            center = re.search(r'^\s*(\d+)', line).group(1)
            # Get the target gas (first one after the colon)
            gases = line.split(':')[-1].split()
            target_gas = gases[0].strip()
            window = '{}_{}'.format(target_gas, center)
            # Get the scale factor
            m = re.search(r'sf=(\d+\.\d+)', line)
            if m is not None:
                sfs[window] = m.group(1)

    return sfs


def get_col_file_window(col_file):
    col_file = Path(col_file).name
    return col_file.split('.')[0]


def advance_col_to_cmd_line(rh, wh=None):
    def get_line():
        line = rh.readline()
        if wh is not None:
            wh.write(line)
        return line

    line = get_line()
    nhead = int(line.split()[0])
    for _ in range(nhead-3):
        # The command line is the second to last in the header and we've already read one line
        # so advance the pointer so that the next line to read is the command line
        get_line()


def col_file_needs_updated(col_file, sfs):
    window = get_col_file_window(col_file)
    if window not in sfs:
        return

    curr_sf = None
    with open(col_file) as rh:
        advance_col_to_cmd_line(rh)
        line = rh.readline()
        sf = re.search(r'sf=(\d+\.\d+)', line).group(1)
        return not sf == sfs[window]

def update_col_file(col_file, sfs, force_update=False):
    if not force_update and not col_file_needs_updated(col_file, sfs):
        return

    window = get_col_file_window(col_file)
    if window not in sfs:
        return 

    backup_file = backup(col_file)
    with open(col_file, 'w') as wh, open(backup_file, 'r') as rh:
        advance_col_to_cmd_line(rh, wh)
        line = rh.readline()
        new_sf_str = 'sf={}'.format(sfs[window])
        line = re.sub(r'sf=\d+\.\d+', new_sf_str, line)
        wh.write(line)

        # Copy the rest
        for line in rh:
            wh.write(line)


def restore_postproc_sh(run_dir, gggpath, force_update=False):
    new_pp_lines = make_new_postproc_sh_lines(run_dir, gggpath)
    pp_path = Path(run_dir) / 'post_processing.sh'
    if not force_update and not does_postproc_sh_need_updated(pp_path, new_pp_lines):
        return

    backup_path = backup(pp_path)
    with open(pp_path, 'w') as wh:
        wh.writelines(new_pp_lines)


def get_runlog_from_col(col_file):
    with open(col_file) as rh:
        for _ in range(5):
            # The runlog should be on the sixth line
            rh.readline()
        rl_line = rh.readline().strip()

    # The first 34 characters are the 32 hex digit hash and two spaces,
    # then the runlog path
    rl_path = rl_line[34:]
    if not rl_path.endswith('.grl'):
        raise IOError('Expected a runlog on the sixth line of {}, but the file name does not end with .grl')
    rl_name = Path(rl_path).stem

    # Cross check with the .mav file
    mav_filename = find_one_file(Path(col_file).parent / '*.mav')
    mav_rl = Path(mav_filename).stem

    if mav_rl != rl_name:
        raise IOError('Runlog name in the O2 .col file is inconsistent with the .mav file name')

    return rl_name


def make_new_postproc_sh_lines(run_dir, gggpath):
    o2_col_file = find_one_file(Path(run_dir) / 'o2_7885.*.col')
    runlog = get_runlog_from_col(o2_col_file)

    template = Path(__file__).parent / '.postproc_template.sh'
    with open(template) as rh:
        new_lines = []
        for line in rh:
            new_lines.append(line.replace('$GGGPATH', str(gggpath)).replace('$RUNLOG', str(runlog)))

    return new_lines
    


def does_postproc_sh_need_updated(ppsh_path, new_pp_lines):
    with open(ppsh_path) as rh:
        old_pp_lines = rh.readlines()

    if len(old_pp_lines) != len(new_pp_lines):
        return True
    else:
        return any(l1 != l2 for l1, l2 in zip(old_pp_lines, new_pp_lines))


def check_gggpath(gggpath):
    gggpath = Path(gggpath)
    mypath = Path(__file__)
    myggg = mypath.resolve().parents[2]
    # I'm in $GGGPATH/utils/python so the third parent should be the gggpath
    if gggpath.stat().st_ino != myggg.stat().st_ino:
        raise EnvironmentError('The environment GGGPATH ({}) is not the same directory as this script is in ({}). Please set your GGGPATH environmental variable properly'.format(gggpath, myggg))


def driver(run_dir, force_update=False, gggpath=None, run_postproc=True):
    run_dir = Path(run_dir)
    gggpath = Path(os.getenv('GGGPATH')) if gggpath is None else Path(gggpath)

    check_gggpath(gggpath)

    multiggg = run_dir / 'multiggg.sh'
    if not multiggg.exists():
        raise IOError('multiggg.sh file not found')

    col_files = sorted(glob(str(run_dir / '*.col')))
    windows_file = gggpath / 'windows' / 'gnd' / 'tccon.gnd'
    
    if not windows_file.exists():
        raise IOError('tccon.gnd windows file not found. Is your GGGPATH set?')

    sfs = read_windows_file(windows_file)
    # We know the CO 4233 window should be removed. Check if that's present to know if they accidentally
    # pointed to an old GGG
    if 'co_4233' in sfs:
        raise IOError('{} looks like an old windows file. Is your GGGPATH set to the Phase 2 GGG2020 repo?'.format(windows_file))

    print('Updating multiggg.sh and .col files as needed...')
    update_multiggg(multiggg, force_update=force_update)
    for col_file in col_files:
        update_col_file(col_file, sfs, force_update=force_update)
    print('Done.')

    restore_postproc_sh(run_dir, gggpath, force_update=force_update)

    if run_postproc:
        print('Running post processing')
        call(['sh', 'post_processing.sh'], cwd=run_dir)


def parse_cl_args():
    p = ArgumentParser(description='Modify a GGG2020 run directory for Phase 2 post-processing, and run post processing')
    p.add_argument('run_dir', default='.', nargs='?', help='The path to the GGG run directory (which contains the .col files and multiggg.sh). The default is the current directory.')
    p.add_argument('-f', '--force-update', action='store_true', help='Force an update of the multiggg.sh and .col files, even if they have already been updated.')
    p.add_argument('-n', '--no-run-postproc', action='store_false', dest='run_postproc', help='Do not automatically run post processing after the update is complete.')
    p.add_argument('--pdb', action='store_true', help='Launch python debugger')

    return vars(p.parse_args())


def main():
    clargs = parse_cl_args()
    if clargs.pop('pdb'):
        import pdb
        pdb.set_trace()
    driver(**clargs)


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        print('ERROR: {}'.format(err), file=sys.stderr)

