#!/usr/bin/env python3
from argparse import ArgumentParser
from difflib import unified_diff
import os
from subprocess import run
import sys


class ScriptError(Exception):
    pass


def main():
    p = ArgumentParser(description='Diff the headers of two interferograms or spectra')
    p.add_argument('orig_file', help='Original file')
    p.add_argument('new_file', help='New file')
    p.add_argument('--pdb', action='store_true', help='Launch the Python debugger')
    
    clargs = vars(p.parse_args())
    if clargs.pop('pdb'):
        import pdb
        pdb.set_trace()
    try:
        diff = diff_opus_headers(**clargs)
    except ScriptError:
        print('Error: {e}', file=sys.stderr)
        sys.exit(1)
    else:
        if diff:
            sys.stdout.writelines(diff)
        else:
            print(f'{clargs["orig_file"]} and {clargs["new_file"]} have identical headers.')

def diff_opus_headers(orig_file, new_file):
    opushdr = find_opus_hdr()
    old_header = get_header_str(opushdr, orig_file)
    new_header = get_header_str(opushdr, new_file)
    return list(unified_diff(old_header, new_header, fromfile=orig_file, tofile=new_file))


def find_opus_hdr():
    gggpath = os.getenv('GGGPATH')
    if gggpath is None:
        raise ScriptError('GGGPATH not set')
    opushdr = os.path.join(gggpath, 'utils', 'OpusHdr', 'OpusHdr')
    if not os.path.isfile(opushdr):
        raise ScriptError(f'OpusHdr program not found at the expected path ({opushdr})')
    return opushdr


def get_header_str(opushdr, path):
    output = run([opushdr, path], capture_output=True)
    # Need splitlines to keep the newlines so that the diff is formatted properly.
    return output.stdout.decode().splitlines(True)


if __name__ == '__main__':
    main()