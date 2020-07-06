#!/usr/bin/env python3
# Note that we deliberately use older functions (e.g. subprocess.call instead of
# subprocess.run, os.path instead of pathlib) to better support old versions of 
# Python 3. Ideally this script should run on *any* version of Python 3.
from argparse import ArgumentParser
from datetime import datetime
import gzip
import hashlib
import json
import os
import re
from subprocess import call
import sys
import tarfile
from textwrap import fill
from urllib import request

MY_DIR = os.path.dirname(__file__)
# Read ~1 MB at a time to avoid eating too much memory
CHUNK_SIZE = 1000000
LINELIST_INFO_FILE = '.linelist_info.json'
PREV_LINELISTS_FILE = '.prev_linelist_info.json'
ALWAYS_KEEP_FILES = {
    os.path.relpath(__file__),
    os.path.relpath(LINELIST_INFO_FILE),
    os.path.relpath(PREV_LINELISTS_FILE),
}

class QuitOnInput(Exception):
    pass

def main():
    p = ArgumentParser(description='Download GGG linelists')
    p.add_argument('--always-backup', action='store_true',
                   help='Automatically back up files in the linelists directory that are NOT old linelists before downloading. '
                        'If this is not given, then you will be prompted to choose whether to back these files up.')
    p.add_argument('--always-delete', action='store_true',
                   help='Automatically delete all files in this directory other that this script, '
                        'the JSON files with the lists of linelists, and hidden files. Without this '
                        'flag, you will be prompted to delete these files.')
    p.add_argument('--pdb', action='store_true', help='Launch the Python debugger')
    clargs = vars(p.parse_args())
    if clargs.pop('pdb'):
        import pdb
        pdb.set_trace()
    try:
        driver(**clargs)
    except QuitOnInput:
        print('Aborting linelist download. Rerun master.sh when ready.', file=sys.stderr)
        sys.exit(2)

def driver(always_delete=False, always_backup=False):
    os.chdir(MY_DIR)

    with open(LINELIST_INFO_FILE) as f:
        linelist_info = json.load(f)
    with open(PREV_LINELISTS_FILE) as f:
        prev_linelists = json.load(f)

    to_delete, to_del_no_back, already_present = check_files(linelist_info, prev_linelists)
    cleanup(to_delete, to_del_no_back, auto_delete=always_delete, auto_backup=always_backup)

    failed_downloads = []
    for download, present in zip(linelist_info, already_present):
        needs_extract = check_one_download(download, present)
        if len(needs_extract) > 0:
            failed = do_one_download(download, needs_extract)
            failed_downloads.extend(failed)

    if len(failed_downloads) > 0:
        failed_str = ', '.join(failed_downloads)
        raise RuntimeError('Some linelists failed MD5 checksum validation ({}), please try again. If the problem persists, open a GitHub issue at https://github.com/TCCON/GGG/issues'.format(
            failed_str,
        ))
        
def gunzip(file):
    outfile = re.sub(r'\.gz$', '', file)
    with gzip.open(file) as g, open(outfile, 'wb') as o:
        chunk = g.read(CHUNK_SIZE)
        while len(chunk) > 0:
            o.write(chunk)
            chunk = g.read(CHUNK_SIZE)

def check_md5(file, md5):
    if not os.path.isfile(file):
        return False
    new_checksum = compute_md5(file)
    return new_checksum == md5


def compute_md5(file):
    checksum = hashlib.md5()
    with open(file, 'rb') as f:
        chunk = f.read(CHUNK_SIZE)
        while len(chunk) > 0:
            checksum.update(chunk)
            chunk = f.read(CHUNK_SIZE)
    return checksum.hexdigest()

def wget(url, outfile, use_external=False, extra_args=tuple()):
    if use_external or extra_args:
        wget_args = ('wget', '--quiet', '--output-document', outfile) + tuple(extra_args) + (url,)
        print('(Calling "', ' '.join(wget_args), '")') 
        call(wget_args)
    else:
        print('(Using URL retrieve.)')
        request.urlretrieve(url, outfile)
        

def cleanup(to_delete, to_del_no_backup, auto_delete=False, auto_backup=True):
    nbak = len(to_delete)
    ndel = len(to_del_no_backup)

    if ndel > 0:
        print('{} files are old versions of linelists and will be removed'.format(ndel))
        for file in to_del_no_backup:
            print(' - {}'.format(file))
            
        if auto_delete or user_yn('Remove these {} files now?'.format(ndel)):
            for file in to_del_no_backup:
                os.remove(file)
                print('Removed {}'.format(file))
        else:
            print(fill('Keeping these files for now; note that installation will FAIL if a new version of one or more of these files is required.'))
                

    if nbak > 0:
        print('{} files are not present in the list of new files:'.format(nbak))
        for path in to_delete:
            print(' - {}'.format(path))
        print('It is recommended to move or remove these files before downloading new linelists.')
        if auto_backup or user_yn('Backup these files?'):
            backup_name = 'backup.{}.tgz'.format(datetime.now().strftime('%Y%m%dT%H%M%S'))
            with tarfile.open(backup_name, 'x:gz') as tar:
                for file in to_delete:
                    tar.add(file)
        
        if auto_delete or user_yn('Remove these {} files now?'.format(nbak)):
            for file in to_delete:
                os.remove(file)
                print('Removed {}'.format(file))
        else:
            print('Keeping these files for now.')

    # TODO: remove empty directories (probably needs to use os.walk to handle subdirectories)
            
        

def user_yn(prompt, default=None, allow_quit=True):
    add_help = False
    quit_opt = ', quit = q' if allow_quit else ''
    if default is None:
        opts = '[yn{}]'.format(quit_opt)
    elif default:
        opts = '[yn{}, default = y]'.format(quit_opt)
    else:
        opts = '[yn{}, default = n]'.format(quit_opt)

    while True:
        if add_help:
            print('Enter "y" or "n" only!', end=' ')
        response = input(f'{prompt} {opts} ').lower()
        if len(response) == 0 and default is not None:
            return default
        if allow_quit and response == 'q':
            raise QuitOnInput()
        if response in {'y', 'yes'}:
            return True
        if response in {'n', 'no'}:
            return False

        add_help = True        
        

def check_files(all_downloads, prev_linelists):
    # Four possibilities:
    #    1. The file is not listed in the new collection of downloads (remove, backup)
    #    2. The file is not listed in the new collection, but its checksum matches an old version (remove, no backup)
    #    3. The file is listed but the checksum does not match any expected checksums (remove, backup)
    #    4. The file is listed and the checksum matches (keep)
    #
    # Note that we use `os.path.relpath` to normalize paths as relative to the current (i.e. linelists)
    # directory. This should help with how os.walk returns '.' as the current dir first, meaning that 
    # os.path.join(currdir, basename) would be e.g. "./atm.101" instead of just "atm.101"; os.path.relpath()
    # makes both into "atm.101" so we can use those strings as keys in a dictionary.
    listed_files = dict()
    for dl_idx, download in enumerate(all_downloads):
        for file_dict in download['contents']:
            listed_files[os.path.relpath(file_dict['file'])] = {'md5': file_dict['md5'], 'index': dl_idx}

    to_delete = []
    to_del_no_backup = []
    already_present = [[] for _ in all_downloads]

    # assume we are in the linelist subdirectory
    for currdir, subdirs, basenames in os.walk('.'):
        for basename in basenames:
            fullname = os.path.relpath(os.path.join(currdir, basename))
            if basename.startswith(('.', 'backup')):
                # Ignore hidden files and backups we've made at any directory level
                continue
        
            if fullname in ALWAYS_KEEP_FILES:
                # Ignore this file and the list of files to download
                continue

            new_checksum = compute_md5(fullname)
            is_old_file = any(fullname == os.path.relpath(p['file']) and new_checksum == p['md5'] for p in prev_linelists)
            
            if fullname not in listed_files or new_checksum != listed_files[fullname]['md5']:
                if is_old_file:
                    to_del_no_backup.append(fullname)
                else:
                    to_delete.append(fullname)
            else:
                already_present[listed_files[fullname]['index']].append(fullname)
            
    return to_delete, to_del_no_backup, already_present


def check_one_download(download, already_present):
    needs_downloaded = []
    wrong_checksum = []
    for check in download['contents']:
        if not os.path.exists(check['file']):
            needs_downloaded.append(check)
        elif check['file'] in already_present:
            print('{} has the correct checksum, will keep the current file'.format(check['file']))
        elif not check_md5(check['file'], check['md5']):
            wrong_checksum.append(check['file'])

    if len(wrong_checksum) == 0:
        return needs_downloaded
    else:
        raise RuntimeError(
            'The following linelist(s) already exist but have the wrong MD5 checksums: {}. Remove these file(s) from the linelists directory and rerun.'.format(
                ', '.join(wrong_checksum)
            )
        )

def do_one_download(download, needed):
    print('Downloading', download['download_url'])
    wget(
        download['download_url'],
        download['output_name'],
        download.get('call_wget', False),
        download.get('extra_wget_args', tuple())
    )
    
    expected_files = [c['file'] for c in needed]

    print('Extracting', download['output_name'])
    if download['output_name'].endswith('.tgz'):
        try:
            tf = tarfile.open(download['output_name'])
            to_extract = [m for m in tf.getmembers() if m.path in expected_files]
            tf.extractall(path='.', members=to_extract)
        finally:
            tf.close()
    elif download['output_name'].endswith('.gz'):
        #call(['gunzip', download['output_name']])
        gunzip(download['output_name'])
    else:
        raise RuntimeError('Unknown file extension on {}'.format(download['output_name']))
    os.remove(download['output_name'])
    
    checksums_failed = []
    for check in needed:
        is_ok = check_md5(check['file'], check['md5'])
        if not is_ok:
            checksums_failed.append(check['file'])

    return checksums_failed


if __name__ == '__main__':
    try:
        main() 
    except RuntimeError as e:
        print(e, file=sys.stderr)
        sys.exit(1)
