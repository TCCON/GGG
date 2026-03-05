#!/usr/bin/env python3
from argparse import ArgumentParser
from glob import glob
import os
import platform
import re
import shutil
from subprocess import PIPE, STDOUT, Popen, check_call, CalledProcessError
import sys
from textwrap import fill

ME = os.path.basename(__file__)
GGGRS_REPO = 'https://github.com/TCCON/ggg-rs.git'
CONDA_NC_LIB = 'libnetcdf=4.9.1'
CONDA_H5_LIB = 'hdf5=1.12.2'

class GggPathError(Exception):
    pass


class QuitOnInput(Exception):
    pass


class SetupFailed(Exception):
    pass


def main():
    p = ArgumentParser(description='Set up your compilation environment for the Rust components of GGG')
    p.add_argument('--pdb', action='store_true', help='Launch Python debugger')

    subp = p.add_subparsers()
    p_setup = subp.add_parser('install', help='Download, setup, and compile GGG-RS (default)')
    #p_setup.add_argument('--non-interactive', action='store_true', help='Error instead of asking for user input (useful for automating)')
    #p_setup.add_argument('--cargo-bin', default='cargo', help='Path to the "cargo" command. The default assumes is it on your PATH.')
    p_setup.add_argument('--allow-wrong-gggpath', action='store_true', help='Continue even if $GGGPATH does not contain the installation directory')
    p_setup.add_argument('--no-clone', action='store_true', help='Do not try to clone or update GGG-RS')
    p_setup.add_argument('--dry-run-compile', action='store_true', help='Run Make in dry run mode for compilation')
    p_setup.set_defaults(driver=install)

    p_list_nc = subp.add_parser('list-nc-opts', help='List options for providing the netCDF library')
    p_list_nc.set_defaults(driver=show_nc_sources)

    clargs = vars(p.parse_args())
    if clargs.pop('pdb'):
        import pdb
        pdb.set_trace()

    driver = clargs.pop('driver', install)
    try:
        driver(**clargs)
    except GggPathError as e:
        print('GGG-RS setup failed: {}. Fix GGGPATH and rerun {} to complete the necessary setup.'.format(e, ME), file=sys.stderr)
    except QuitOnInput:
        print('Aborting {} for now. Rerun when ready to complete the necessary setup.'.format(ME), file=sys.stderr)
        sys.exit(2)
    except SetupFailed as e:
        print('GGG-RS setup failed: {}. Correct the underlying problem and rerun {} to finish setup.'.format(e, ME), file=sys.stderr)
        sys.exit(1)



def install(cargo_bin=None, no_clone=False, allow_wrong_gggpath=False, nc_source=None, dry_run_compile=False, **kwargs):
    acknowledge()
    env_vars = setup(cargo_bin=cargo_bin, no_clone=no_clone, allow_wrong_gggpath=allow_wrong_gggpath, nc_source=nc_source, **kwargs)
    compile(env_vars, allow_wrong_gggpath, dry_run=dry_run_compile)


def acknowledge():
    msg = ('GGG-RS programs are still experimental, and will remain so until GGG2020.2 at the earliest. '
           'You should only use them to process standard TCCON or EM27/SUN data if you have communicated with '
           'the algorithm team and been instructed to do so. Use for non-standard TCCON or EM27/SUN data is '
           'at your discretion.')
    print(fill(msg))
    if not user_yn('Enter y to acknowledge these caveats and continue with setup', allow_quit=False):
        raise QuitOnInput()


def setup(cargo_bin=None, no_clone=False, allow_wrong_gggpath=False, nc_source=None, **kwargs):
    notes = []
    notes.append(clone_or_update_gggrs(allow_wrong_gggpath=allow_wrong_gggpath, no_clone=no_clone))
    cargo_bin = get_cargo_bin(cargo_bin)
    nc_source_class = choose_nc_source()
    nc_source = nc_source_class(**kwargs)
    nc_source.setup(allow_wrong_gggpath)
    notes.append(nc_source.final_user_notes())

    notes = [n for n in notes if n]
    if notes:
        print('\n====== NOTE ======\n')
        for idx, note in enumerate(notes, start=1):
            print('{:3d}. '.format(idx), end='')
            print(fill(note, subsequent_indent='     '))

    env_vars = nc_source.env_vars(gggpath_mismatch=allow_wrong_gggpath)
    env_vars['CARGOCMD'] = cargo_bin
    return env_vars


def compile(env_vars, allow_wrong_gggpath=False, dry_run=False):
    gp = gggpath(allow_wrong_gggpath)
    log_file = os.path.join(gp, 'install', 'compile_rs_messages.out')
    working_dir = _gggrs_dir(allow_wrong_gggpath)
    working_env = os.environ.copy()
    working_env.update(env_vars)
    if dry_run:
        print('--- Start compilation dry run ---')
        p = Popen(['make', 'check-args'], cwd=working_dir, env=working_env)
        p.wait()
        print('')
        p = Popen(['make', 'install', '--dry-run'], cwd=working_dir, env=working_env)
        p.wait()
        print('--- End compilation dry run ---')
    else:
        print(f'Compiling GGG-RS. Compilation messages will be piped to {log_file}')
        with open(log_file, 'w') as f:
            p = Popen(['make', 'install'], stdout=f, stderr=STDOUT, cwd=working_dir, env=working_env)
            p.wait()
            if p.returncode != 0:
                raise SetupFailed(f'GGG-RS compilation failed. See {log_file} to diagnose why.')
        bin_dir = os.path.join(gp, 'bin')
        print(f'GGG-RS compilation complete. Its programs should have been added to {bin_dir}')


def clone_or_update_gggrs(allow_wrong_gggpath, no_clone=False):
    gggrs_dest = _gggrs_dir(allow_wrong_gggpath)
    gggrs_exists = os.path.isdir(gggrs_dest)
    if no_clone and not gggrs_exists:
        return 'GGG-RS was not found at {}, rerun {} without --no-clone to download'.format(gggrs_dest, ME)
    elif no_clone:
        return 'GGG-RS was not updated.'
    elif gggrs_exists:
        _update_gggrs(gggrs_dest)
    else:
        _clone_gggrs(gggrs_dest)


def _gggrs_dir(allow_wrong_gggpath):
    ggg_root = gggpath(allow_wrong_gggpath)
    return os.path.join(ggg_root, 'src-rs')


def _clone_gggrs(gggrs_dest):
    print('Cloning GGG-RS into {}'.format(gggrs_dest))
    try:
        check_call(['git', 'clone', GGGRS_REPO, gggrs_dest])
    except CalledProcessError as err:
        raise SetupFailed('could not clone GGG-RS (reason was: {})'.format(err))
    else:
        print('Clone successful.')

def _update_gggrs(gggrs_dest):
    print('Checking if GGG-RS needs updated.')

    p = Popen(['git', 'status', '--porcelain'], stdout=PIPE, cwd=gggrs_dest)
    files = [line.split()[1] for line in p.stdout.read().decode().splitlines()]
    if len(files) > 0:
        prompt = '{} untracked or modified files are present in {}.'.format(len(files), gggrs_dest)
        choices = ['Pull anyway', 'Do not pull, proceed with setup']
        opt = user_select(prompt, choices)
        if opt == 1:
            print('Not updating GGG-RS.')
            return

    p = Popen(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], stdout=PIPE, cwd=gggrs_dest)
    branch = p.stdout.read().decode().strip()
    if branch != 'main':
        raise SetupFailed('GGG-RS in {} is not on the "main" branch (current branch = {}). If this is expected, rerun {} with the --no-clone option.'.format(gggrs_dest, branch, ME))

    try:
        check_call(['git', 'pull', '--ff-only'], cwd=gggrs_dest)
    except CalledProcessError as err:
        raise SetupFailed('unable to pull updates for ggg-rs (reason was: {}). If this is expected, rerun {} with the --no-clone option.'.format(err, ME))


def get_cargo_bin(user_input_cargo):
    if user_input_cargo is not None:
        print('Using provided path to cargo, {}'.format(user_input_cargo))
        return user_input_cargo

    if which('cargo') is not None:
        print('"cargo" found on PATH')
        return 'cargo'

    prologue = ('"cargo" was not found on your PATH. If on a managed computing system '
                'like an HPC, you may need to load a module to make a Rust toolchain available. '
                'If a Rust toolchain is not yet installed on your system, one can be obtained '
                'from https://rustup.rs/. If one is already installed but not on your PATH, '
                'enter the path to "cargo" here.')
    print(fill(prologue))
    return user_input_path('Path to cargo')


def choose_nc_source():
    nc_source_preference = NcSource.method_preference()
    available = []
    unavailable = []
    for src in nc_source_preference:
        can_use, why_not = src.can_use()
        if can_use:
            available.append(src)
        else:
            unavailable.append(' - {} ({})'.format(src.label()[0], why_not))

    choices = ['{} - {}'.format(*src.label()) for src in available] + ['Show details']
    epilogue = None
    if len(unavailable) > 0:
        unavail_str = '\n'.join(unavailable)
        epilogue = 'Some options are not available:\n{}'.format(unavail_str)

    while True:
        choice = user_select('Select option for the netCDF programs', choices, epilogue)
        if choice < len(available):
            return available[choice]
        else:
            for src in nc_source_preference:
                is_avail, _ = src.can_use()
                avail_str = '' if is_avail else ' [UNAVAILBLE]'
                descr = src.help_text().splitlines()
                descr = [line.strip() for line in descr]
                descr = fill(' '.join(descr), initial_indent='  ', subsequent_indent='  ')
                name, short = src.label()
                print('{name}{avail_str} - {short}'.format(name=name, avail_str=avail_str, short=short))
                print(descr, end='\n\n')


def show_nc_sources():
    for source in NcSource.method_preference():
        can_use, why_not = source.can_use()
        avail_str = '[ AVAILABLE ]' if can_use else '[UNAVAILABLE]'
        name, short = source.label()
        print('{} {} - {}'.format(avail_str, name, short))
        if not can_use:
            print(' > {}'.format(why_not))

class NcSource:
    def __init__(self, **_):
        pass

    @classmethod
    def method_preference(cls):
        return (NcConda, NcSystem, NcStatic, NcNone)


    @classmethod
    def can_use(cls):
        """Return if this option is available.
        Returns both a boolean and a string; the string explains
        why this is not available if ``False`` and should just be
        empty if ``True``.
        """
        return True, ''

    @classmethod
    def label(cls):
        """Label to use in menus, etc."""
        return "BUG!", "This should not show up!"

    @classmethod
    def help_text(cls):
        """Return a string that describes how this class configures
        ggg-rs to build it's netCDF dependency. Multiline strings will
        have newlines removed and be rewrapped to ~70 characters wide.
        """
        return "Displaying this help text is a bug!"

    def setup(self, gggpath_mismatch):
        """Do any actions needed to set up 
        """
        pass

    def env_vars(self, gggpath_mismatch):
        """Return a dictionary of environmental variable name-value pairs
        that will override existing environmental variables when the
        GGG-RS Makefile is executed.
        """
        return dict()

    def final_user_notes(self):
        """Print any messages to the user that should go at the end of the program
        so that they are easier to see.
        """
        return ''


class NcNone(NcSource):
    def __init__(self, **_):
        pass

    @classmethod
    def help_text(cls):
        return """Compile ggg-rs without netCDF support. This is the simplest
        configuration, but does not provide any of the netCDF-related programs,
        such as write_private_netcdf. Therefore, this option is not viable for
        TCCON, EM27, or other users that must deliver the .private.nc files for
        distribution.
        """

    @classmethod
    def label(cls):
        return "None", "Do not compile netCDF programs"

    def env_vars(self, gggpath_mismatch):
        return {'GGGRS_FEATURES': '', 'GGGRS_NCDIR': 'NONE'}

    def final_user_notes(self):
        return ('If you have previously compiled the netCDF programs, you may still find that their '
                'modification times in $GGGPATH/bin are updated after compiling finished. This seems '
                'to be because cargo will copy compiled programs even if they were not updated by the '
                'current run.') 

class NcSystem(NcSource):
    def __init__(self, netcdf_dir=None, allow_missing_hdf5=False, **_):
        self._netcdf_dir = netcdf_dir
        self._allow_missing_hdf5 = allow_missing_hdf5

    @classmethod
    def label(cls):
        return "System", "Use existing HDF5 and netCDF4 libraries"

    @classmethod
    def help_text(cls):
        return """Compile ggg-rs using an HDF5 and netCDF4 library already available
        on your system. This is probably most useful for HPCs and other managed computing
        environments that have many scientifically-relevant software libraries installed."""

    def setup(self, gggpath_mismatch):
        if self._netcdf_dir is None:
            prologue = ('We need the path for NETCDF_DIR. This will almost always be a directory '
                        'with a "lib" or "lib64" subdirectory containing a file whose name starts with '
                        '"libnetcdf". To support GGG-RS, it must also have the HDF5 library, meaning a '
                        'file whose name starts with "libhdf5" On an HPC, you may need to load a module '
                        'to have this library available.')
            nc_config_path = which('nc-config')
            if nc_config_path is not None:
                p = Popen([nc_config_path, '--prefix'], stdout=PIPE)
                prefix = p.stdout.read().decode()
                prologue += ' (Hint: the nc-config on your PATH thinks that the correct value for NETCDF_DIR is "{}".)'.format(prefix)

            prologue = fill(prologue)
            print(prologue)
            self._netcdf_dir = user_input_path('NETCDF_DIR')

        # Check if libhdf5 is also present. Allow for .a (static archive), .so (shared objects on Linux),
        # and .dylib (dynamic libraries on Mac).
        libhdf5_candidates = (
            list(glob(os.path.join(self._netcdf_dir, 'lib', 'libhdf5*'))) + 
            list(glob(os.path.join(self._netcdf_dir, 'lib64', 'libhdf5*')))
        )
        found_hdf5 = any(os.path.basename(f).endswith(('.a', '.so', '.dylib')) for f in libhdf5_candidates)
        if not found_hdf5:
            msg = (f'There was no "libhdf5*.a", "libhdf5*.so", or "libhdf5*.dylib" file under {self._netcdf_dir}, '
                   'in neither the "lib" nor "lib64" subdirectories')
            if self._allow_missing_hdf5:
                full_msg = f'WARNING: {msg}. GGG-RS may fail to compile.'
                print(fill(full_msg))
            else:
                raise SetupFailed(msg)



    def env_vars(self, gggpath_mismatch):
        # We should only need to communicate to the Makefile that the netCDF libraries
        # exist here; it should handle everything else.
        return {'GGGRS_NCDIR': self._netcdf_dir}

    def final_user_notes(self):
        if platform.system() == 'Darwin':
            dylib_note = ('Since you appear to be running on a Mac, you may need to include {0}/lib in your DYLD_FALLBACK_LIBRARY_PATH '
                          'environmental variable, e.g. "export DYLD_FALLBACK_LIBRARY_PATH={0}/lib" in bash and have this variable set '
                          '*when running GGG programs* (not just when compiling).').format(self._netcdf_dir)
            return dylib_note

class NcConda(NcSource):
    def __init__(self, program=None, clean=None, **_):
        programs = ['mamba', 'micromamba', 'conda']
        self._program_options = [p for p in programs if which(p) is not None]
        self._program = program
        self._clean = clean

    @classmethod
    def can_use(cls):
        # We'll give the user a chance to specify a custom program during
        # setup, so always allow use to try this class
        return True, ''

    @classmethod
    def label(cls):
        return 'Conda', 'Install HDF5 and netCDF4 with conda, mamba, or micromamba'

    @classmethod
    def help_text(cls):
        return """Install the HDF5 and netCDF4 libraries in a Conda environment
        at and use these to build the netCDF ggg-rs programs. This is the
        recommended approach for TCCON and EM27 users."""

    def setup(self, gggpath_mismatch):
        if self._program is None:
            i = user_select('Which program should we use to manage the environment', self._program_options)
            self._program = self._program_options[i]

    def env_vars(self, gggpath_mismatch):
        return {
            'GGGRS_NCDIR': 'AUTO',  # this tells the GGG-RS installation to manage its own environment in its own directory
            'GGG_ENV_TOOL': self._program
        }

    def final_user_notes(self):
        if platform.system() == 'Darwin':
            lib_dir = os.path.join(gggpath(True), 'src-rs', '.condaenv')
            dylib_note = ('Since you appear to be running on a Mac, you may need to include {0}/lib in your DYLD_FALLBACK_LIBRARY_PATH '
                          'environmental variable, e.g. "export DYLD_FALLBACK_LIBRARY_PATH={0}/lib" in bash and have this variable set '
                          '*when running GGG programs* (not just when compiling).').format(lib_dir)
            return dylib_note

class NcStatic(NcSource):
    def __init__(self, cmake_path=None, **_):
        self._cmake_path = cmake_path
        self._include_cmake_envvar = None if cmake_path is None else True

    @classmethod
    def can_use(cls):
        # We'll give the user a chance to provide a cmake path during setup
        return True, ''

    @classmethod
    def label(cls):
        return "Static", "Build HDF5 and netCDF4 from source"

    @classmethod
    def help_text(cls):
        return """Compile netCDF and HDF5 libraries from source.
        This makes use of the netcdf and hdf5 crates' ability to
        compile their respective C libraries from source and link
        to that instance. This requires a somewhat-recent version
        of cmake (>= 3 probably, but this has not been rigorously
        tested) to be installed.  This increases compile times, and
        is only recommended if cmake is easier to use on your system
        than a Python environment with the conda package manager installed.
        """

    def setup(self, gggpath_mismatch):
        if self._cmake_path is None:
            cmake = which('cmake')
            cmake_on_path = True
        else:
            cmake = self._cmake_path
            cmake_on_path = False

        if cmake is None:
            cmake = user_input_path('Enter the path to your cmake program')
            cmake_on_path = False

        self._include_cmake_envvar = not cmake_on_path

        p = Popen([cmake, '--version'], stdout=PIPE)
        v = p.stdout.read()
        m = re.search(rb'(\d+)\.\d+(\.\d+)?', v)
        if m is None:
            print('Note: could not determine cmake version')
        else:
            major = int(m.group(1))
            version = m.group().decode()
            if major < 3:
                print('Warning: cmake version ({vers}) is < 3; compiling HDF5 and netCDF4 from source may fail'.format(vers=version))


    def env_vars(self, gggpath_mismatch):
        env_vars = {'GGGRS_NCDIR': 'STATIC'}
        if self._include_cmake_envvar:
            env_vars['CMAKE'] = self._cmake_path
        return env_vars


def gggpath(allow_mismatch, subdir=None):
    try:
        base = os.environ['GGGPATH']
    except KeyError:
        raise GggPathError('GGGPATH environmental variable must be set')

    if not allow_mismatch:
        # assume this file is in $GGGPATH/install
        mydir = os.path.dirname(__file__)
        mydir_real = os.path.realpath(mydir)
        gggdir_real = os.path.join(os.path.realpath(base), 'install')
        if mydir_real != gggdir_real:
            raise GggPathError('GGGPATH ({}) does not contain current install directory ({}), aborting'.format(base, mydir_real))

    if subdir:
        return os.path.join(base, subdir)
    else:
        return base

def which(prog):
    if hasattr(shutil, 'which'):
        return shutil.which(prog)
    else:
        # backwards compatibility before shutil.which was
        # added, but probably less cross-platform compatible
        p = Popen(['which', prog], stdout=PIPE)
        path = p.stdout.read().decode()
        if len(path) == 0:
            return None
        else:
            return path


def user_input_path(prompt, blank_return_none=False, allow_quit=True):
    while True:
        if allow_quit:
            path = input(f'{prompt} (enter q to quit): ')
        else:
            path = input(f'{prompt}: ')
        if len(path) == 0 and blank_return_none:
            return None
        elif len(path) == 0:
            print('Path cannot be empty!', end=' ')
            continue
        elif allow_quit and path == 'q':
            raise QuitOnInput()
        path = os.path.expanduser(path)
        path = os.path.expandvars(path)
        if os.path.exists(path):
            return path
        elif user_yn(f'{path} does not exist, use anyway?'):
            return path

def user_select(prompt, choices, epilogue=None, allow_quit=True):
    n = len(choices)
    print(prompt)
    for i, c in enumerate(choices, start=1):
        print(' {i}: {c}'.format(i=i, c=c))
    if epilogue:
        print(epilogue)
    while True:
        if allow_quit:
            ans = input('Enter 1-{} or q to quit: '.format(n))
        else:
            ans = input('Enter 1-{}: '.format(n))
        if allow_quit and ans == 'q':
            raise QuitOnInput()
        try:
            idx = int(ans)
        except ValueError:
            print('"{}" is not a number.'.format(ans), end=' ')
            continue

        if idx < 1 or idx > n:
            print('{} is not an available choice'.format(idx), end=' ')
        else:
            return idx - 1


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


if __name__ == '__main__':
    main()
