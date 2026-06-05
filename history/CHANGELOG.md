# GGG Changelog

This file is distinct from the `ggg.history` file in that this file summarizes changes between
releases, whereas the `ggg.history` file details individual changes to specific programs.

## GGG2020.1 changes

Note: these changes are relative to the original GGG2020 release through the Mercurial repo,
some of these changes were already included in the GitHub repo.


### Installation

- Now supports use of `conda`, `micromamba`, or `pip` to manage the Python environment.
    - Note that `conda` or `micromamba` are preferred; `pip` is intended only for cases
      where those are not available.
- The Python environment will now be written to `$GGGPATH/install/.condaenv` or `$GGGPATH/install/.venv`,
  rather than a named environment in the central conda directory. This simplifies managing multiple GGG
  installations.
- `check_python.sh` no longer uses a login bash shell.
- The netCDF writer is now cloned from the TCCON GitHub organization.
- `master.sh` renamed to `runme_install_main.sh`, `pymaster.sh` to `run_pyinstall.sh`, and `i2s_master.sh` to `run_i2s_test.sh` to move away from "master" terminology.
- `.mod` and `.vmr` files used in the benchmark have been updated to include the new information about GEOS versions in the header.
- `download_linelists.py` modified to handle updates to caltechData, which stores the linelists.
    - It will now use `requests` if available and has multiple URLs from which to obtain the linelists.

### Interferogram processing

- I2S now requires an extra option to indicate whether to apply a nonlinearity correction and (if so)
  what coefficients to use. The values from this option are also written into the spectrum headers.
    - The example input files have also been updated.
- I2S now compiled with the `-cpp` flag to ensure pre-processing is used where needed.
- Fixed a bug that was outputting the wrong value in the logging information.
- Added Dave Griffith's fix for an incorrect `P2L` header value for firmware version 2.485.
- Added file path length increases requested by Jacob Hedelius and Aaron Meyer
- Added extra checks in I2S for uncommon EM27 instrument strings.


### L2 setup

- Added an option to `utils/python/list_mod_vmr_links` to make links for the alternate
  priors before 1 Apr 2024.
- Added an option to `gsetup` to create a run directory with EM27 post processing.
- Modified `gsetup` to use a new Rust version of `collate_results` if the `GGG_RS_POSTPROC`
  environmental variable is set to `1`.
- Also modified `gsetup` to allow setting up the `post_processing.sh` script with EM27/SUN
  specific options (primarily for `collate_results`).
- Updated the benchmark `.vmr` files to include the newly required `CO_SOURCE` header entry.
- Updated the example InSb windows files to add an `fco2` window (for channel fringes) and
  remove some mid-IR CO windows.
- Updated `create_sunrun_from_Wgong`.
- Updated `create_sunrun_from_darwin_ifs2` to include pressure correction.


### L2 retrieval

- Improved averaging kernel computation: GFIT now outputs Jacobian files for each window
  which must now be concatenated with the new `concat_jac` program. `avg_ker` was updated
  to take the new concatenated Jacobian files as input.
    - This process was used to generate the pre-computed GGG2020 averaging kernels, but
      these programs were not included in the original GGG2020 release.
- Fixed an issue with `gfit` that caused runs to fail on systems that use process ID with 7 or 8 digits.
- Increased output precision of `zmin`.


### Post processing

- Increased maximum number of auxiliary columns to 26 to allow including the O2 global mean DMF as a new auxiliary
  column in the future.
- `collate_results`:
    - For data with spectra following the TCCON naming convention, will now group adjacent observations
      if the spectrum names indicate the observations should be at the same time, even if the ZPD time
      is different - this attempts to address issues with the second detector's ZPD time not matching the
      first detector's.
    - Now accepts `"em27"` as a second command line argument to use EM27/SUN-specific settings.
- Added a bug fix from Jacob Hedelius and Aaron Meyer in `write_aux` to handle paths with periods other than
  at the file extension
- Added various file path length increases requested by Jacob Hedelius and Aaron Meyer
- Updated `xch4`, `zmin`, `fovi`, and `pout` precision in the default `pa_qc.dat` file.
- netCDF files will now be created in the "2020.C" format by default, see https://tccon-wiki.caltech.edu/Main/GGG2020DataChanges#File_format_GGG2020.C
  for a list of the changes.
    - TCCON sites will upload this 2020.C format file to Caltech, and conversion to 2020.1.A format will happen
      in the central data pipeline. This allows us to ensure that the same time-varying O2 mole fraction is used
      for all data.
    - For non-TCCON users (e.g., EM27/SUNs), the program to update to GGG2020.1.A format is included
      (`update_file_format`)

### Utilities

- Added a new utility to difference the headers of two interferograms or spectra
- Added a utility program to transfer updates from the private Mercurial repo to the public Git repo
