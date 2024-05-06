#!/usr/bin/env python3
import gzip
import hashlib
import json
import os
import re
from subprocess import call
import sys
import tarfile
from urllib import request

MY_DIR = os.path.dirname(__file__)
# Read ~1 MB at a time to avoid eating too much memory
CHUNK_SIZE = 1000000

def main():
    os.chdir(MY_DIR)

    with open('.linelist_info.json') as f:
        linelist_info = json.load(f)

    failed_downloads = []
    for download in linelist_info:
        needs_extract = check_one_download(download)
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
    checksum = hashlib.md5()
    with open(file, 'rb') as f:
        chunk = f.read(CHUNK_SIZE)
        while len(chunk) > 0:
            checksum.update(chunk)
            chunk = f.read(CHUNK_SIZE)
    return checksum.hexdigest() == md5


def wget(url, outfile, use_external=False, extra_args=tuple()):
    if use_external or extra_args:
        wget_args = ('wget', '--quiet', '--output-document', outfile) + tuple(extra_args) + (url,)
        print('(Calling "', ' '.join(wget_args), '")') 
        call(wget_args)
    else:
        print('(Using URL retrieve.)')
        request.urlretrieve(url, outfile)
        

def check_one_download(download):
    needs_downloaded = []
    wrong_checksum = []
    for check in download['contents']:
        if not os.path.exists(check['file']):
            needs_downloaded.append(check)
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
