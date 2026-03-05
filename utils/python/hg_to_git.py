#!/usr/bin/env python
from argparse import ArgumentParser
import hashlib
import json
import os
from pathlib import Path
import shutil
from subprocess import run

DIR_PLACEHOLDER_FILE = '.dir_placeholder_for_hg'
MY_GGG_PATH = Path(__file__).resolve().parents[2]

WHY_FAILED_NO_FILE = 'does not exist'
WHY_FAILED_EXISTS = 'already exists'
WHY_FAILED_CANT_PATCH = 'cannot patch'

ACT_DIFF = 'Diff files'
ACT_VIEW_REJ = 'Inspect .rej file'
ACT_VIEW_HG = 'Inspect Mercurial copy'
ACT_VIEW_GIT = 'Inspect Git copy'
ACT_COPY_HG_TO_GIT = 'Copy Mercurial file to Git repo (removes .rej file if present)'
MARK_COPY_HG_TO_GIT = 'copied-hg-to-git'
ACT_IS_OK = 'Leave as is, mark as resolved (removes .rej file if present)'
MARK_IS_OK = 'left-as-is'
ACT_NOTHING = 'Follow up afterwards (leaves .rej file if present)'

STEPS = {
    'precheck': -100,
    'make_patch': 10,
    'modify_patch': 20,
    'check_patch': 30,
    'apply_patch': 40,
    'handle_failures': 50,
    'postcheck': 100,
}

def main():
    clargs = parse_args()
    if clargs.pop('pdb'):
        import pdb
        pdb.set_trace()

    driver(**clargs)


def parse_args():
    p = ArgumentParser(description='Program to help export changes from the Mercurial GGG repo to the GitHub one')
    p.add_argument('--git-ggg', required=True, type=Path, help='Path to the Git GGG repo to export to (required)')
    p.add_argument('--hg-start', required=True, help='Tag or other ID to difference from in the Mercurial repo (required)')
    p.add_argument('--hg-end', default='tip', help='Tag or other ID to end the difference for the patch in the Mercurial repo, default = "%(default)s"')
    p.add_argument('--hg-ggg', default=MY_GGG_PATH, type=Path, help='Path to the Mercurial GGG repo to export from, default = "%(default)s"')
    p.add_argument('--allow-hg-unclean', action='store_true', help='Allow the Mercurial repo to be in a non-clean state (including untracked files)')
    p.add_argument('--allow-git-unclean', action='store_true',
                   help='Allow the Git repo to be in a non-clean state (including untracked files). Use with care, as it will be VERY hard to recover '
                        'the pre-patch status.')
    p.add_argument('--start-step', choices=STEPS, default='precheck', help='First step to perform, earlier steps are skipped')
    p.add_argument('--end-step', choices=STEPS, default='postcheck', help='Last step to perform, later steps are skipped')
    p.add_argument('--pdb', action='store_true', help='Launch the Python debugger')
    p.epilog = (
        'This program uses the EDITOR and DIFFPROG environmental variables to determine how to display files and file differences, '
        'respectively, when needed. If EDITOR is not set, files are just printed to the screen. If DIFFPROG is not set, it defaults '
        'to vimdiff. If you want to use other programs, be sure to set these variables.'
    )
    return vars(p.parse_args())


class DoStep:
    def __init__(self, start_step: str, end_step: str):
        self.start_int = STEPS[start_step]
        self.end_int = STEPS[end_step]

    def do_step(self, step_name):
        return self.start_int <= STEPS[step_name] and self.end_int >= STEPS[step_name]


def driver(git_ggg: Path, hg_ggg: Path, hg_start: str, hg_end: str, allow_hg_unclean: bool = False, allow_git_unclean: bool = False,
           start_step='precheck', end_step='postcheck'):
    git_ggg = git_ggg.resolve()
    hg_ggg = hg_ggg.resolve()
    initial_patch_file = hg_ggg / 'initial.patch'
    final_patch_file = hg_ggg / 'final.patch'
    git_check_file = git_ggg / 'git-apply.rpt'

    step_manager = DoStep(start_step, end_step)

    final_messages = []
    if step_manager.do_step('precheck'):
        precheck(hg_ggg=hg_ggg, git_ggg=git_ggg, skip_hg_check=allow_hg_unclean, skip_git_check=allow_git_unclean)

    if step_manager.do_step('make_patch'):
        create_initial_patch_file(hg_ggg=hg_ggg, hg_start=hg_start, hg_end=hg_end, patch_file=initial_patch_file)

    if step_manager.do_step('modify_patch'):
        modify_patch_file(
            initial_file=initial_patch_file,
            hg_ggg=hg_ggg, modified_file=final_patch_file,
            final_messages=final_messages
        )

    if step_manager.do_step('check_patch'):
        do_patch_check(git_ggg=git_ggg, patch_file=final_patch_file, report_file=git_check_file)

    if step_manager.do_step('apply_patch'):
        do_patch_apply(git_ggg=git_ggg, patch_file=final_patch_file)

    if step_manager.do_step('handle_failures'):
        resolutions_file = git_ggg / 'hg-to-git-resolutions.rpt'
        handle_failed_files(
            report_file=git_check_file,
            hg_ggg=hg_ggg, git_ggg=git_ggg,
            resolution_file=resolutions_file,
            final_messages=final_messages
        )

    if step_manager.do_step('postcheck'):
        postcheck_file = git_ggg / 'hg-to-git.rpt'
        postcheck(hg_ggg=hg_ggg, git_ggg=git_ggg, report_file=postcheck_file, final_messages=final_messages)


    if final_messages:
        print('\nIMPORTANT THINGS TO NOTE:')
        for msg in final_messages:
            print(f' - {msg}')


def precheck(hg_ggg: Path, git_ggg: Path, skip_hg_check=False, skip_git_check=False):
    """Check that the Git and Mercurial repos are in a fit state to do the patch.
    """
    if not skip_git_check and not _git_is_clean(git_ggg):
        # We want the git repo to be clean (no uncommitted changes) because if the patch goes wrong,
        # we can revert to the last commit easily with a `git restore . && git clean -fd`. But if there
        # are changes we haven't saved, those would get lost too.
        raise RuntimeError('Git GGG repo is not in a clean state. Restore it to one, or override this check.')
    if not skip_hg_check and not _hg_is_clean(hg_ggg):
        # It's not as critical that the Mercurial repo be in a clean state. But by checking that it is, we can
        # be more sure that the patch will capture all of the changes we want.
        raise RuntimeError('Mercurial GGG repo is not in a clean state. Restore it to one, or override this check.')


def create_initial_patch_file(hg_ggg: Path, hg_start: str, hg_end: str, patch_file: Path):
    """Create the patch file containing all of the Mercurial changes between the two given changesets.
    """
    with open(patch_file, 'w') as f:
        run(
            ['hg', 'diff', '--git', '-r', hg_start, '-r', hg_end],
            stdout=f, cwd=hg_ggg
        )


def modify_patch_file(initial_file: Path, hg_ggg: Path, modified_file: Path, final_messages: list):
    """Modify the patch file to account for differences between the Mercurial and Git repos.
    """
    with open(initial_file) as i, open(modified_file, 'w') as f:
        for diff_block in _iter_patch_blocks(i):
            old_file, new_file = _get_block_files(diff_block)

            if '.hgignore' in old_file or '.hgignore' in new_file:
                # Obviously, git doesn't respect .hgignore files. It also seems
                # that Mercurial really wants you to use a single root .hgignore file -
                # you can put them in subdirectories, but you have to explicitly declare
                # them in the root file. Worse, Mercurial uses regular expressions and
                # Git uses globs by default, so we're almost always going to have to
                # handle this manually.
                diff_block = _replace_hgignore(diff_block)
                final_messages.append(f'.hgignore file changes (a/{old_file}, b/{new_file}) applied to .gitignore. Check that patterns are compatible with Git (Mercurial uses regex patterns by default).')

            elif Path(old_file).name.startswith('.hg') or Path(new_file).name.startswith('.hg'):
                # Other Mercurial files (like .hgtags) just shouldn't be copied at all.
                diff_block = []
                final_messages.append(f'Skipped changes to a .hg* file (a/{old_file}, b/{new_file})')

            elif any(DIR_PLACEHOLDER_FILE in name for name in [old_file, new_file]):
                # GGG traditionally uses .dir_placeholder_for_hg files to keep empty directories.
                # Since Mercurial doesn't make it easy to put .hgignore files in subdirectories,
                # we'll keep following that tradition for now. In the Git repo, I replace these
                # placeholder files with .gitignore files that tell Git to ignore everything else
                # in that directory, this should modify the diff block to do that, if this is a new file.
                diff_block = _replace_dir_placeholder(diff_block)

            elif _is_linelist_file(old_file) or _is_linelist_file(new_file):
                # Don't check in linelist files to Git, they'll burn through the GitHub repo size limits.
                # Instead, we use CaltechData as a manual substitute for Git LFS.
                _check_linelist_file(new_file, hg_ggg=hg_ggg, final_messages=final_messages)
                diff_block = []
            f.writelines(diff_block)


def do_patch_check(git_ggg: Path, patch_file: Path, report_file: Path) -> dict:
    """Check if the patch will apply cleanly and write any files that don't to a report.

    It's much easier to parse this check output that the full output of `git apply`. We will
    use this output later to manually address those files that fail.
    """
    with open(report_file, 'w') as rpt:
        run(
            ['git', 'apply', '--check', str(patch_file.resolve())],
            cwd=git_ggg, stderr=rpt
        )


def do_patch_apply(git_ggg: Path, patch_file: Path) -> list:
    """Apply the patch to the Git repository, allowing files to fail if the patch can't be applied cleanly.
    """
    msg = '* Applying patch *'
    banner = '*' * len(msg)
    print(f'{banner}\n{msg}\n{banner}\n')
    # Use the --reject flag so that we get .rej files when the patch cannot be applied;
    # we'll use those to help decide manual fixes.
    run(
        ['git', 'apply', '--reject', str(patch_file.resolve())],
        cwd=git_ggg
    )


def handle_failed_files(report_file: Path, hg_ggg: Path, git_ggg: Path, resolution_file: Path, final_messages: list):
    """Review the files that ``git apply --check`` reported would fail.

    This will go through the report produced by ``do_patch_check`` and give the user
    an opportunity to examine each failed file. How each one is handled is written to
    a JSON file so that rerunning this step later doesn't force you to go over those
    again, plus it provides traceability."""
    with open(report_file) as f:
        failure_lists_from_check = _parse_git_apply_check(f.read())

    if Path(resolution_file).exists():
        with open(resolution_file) as f:
            resolutions = json.load(f)
    else:
        resolutions = dict()

    file_to_reason = dict()
    n_files = 0
    for key, file_list in failure_lists_from_check.items():
        file_to_reason.update({f: key for f in file_list})
        n_files += len(file_list)
    if n_files != len(file_to_reason):
        raise NotImplementedError('Some files failed to patch for multiple reasons')

    for file in sorted(file_to_reason.keys()):
        if file in resolutions:
            print(f'{file} previously resolved')
            continue

        reason = file_to_reason[file]
        reject_file = git_ggg / f'{file}.rej'
        if not reject_file.exists():
            reject_file = None

        print('')
        if reason == WHY_FAILED_NO_FILE:
            print(f'File "{file}" did not exist in the Git repository, but the patch tried to modify it.')
        elif reason == WHY_FAILED_EXISTS:
            print(f'File "{file}" already existed in the Git repository, but the patch tried to create it as a new file.')
        elif reason == WHY_FAILED_CANT_PATCH:
            print(f'The patch could not be applied to file "{file}".')
        else:
            print(f'File "{file}" failed for an unexpected reason: {reason}')

        # For now, we always resolve the issues the same way. If we want more flexibility in 
        # the future, we can put this call into the if block.
        mark = _compare_and_maybe_copy_failed_file(
            failed_file=file,
            hg_ggg=hg_ggg,
            git_ggg=git_ggg,
            reject_file=reject_file,
            final_messages=final_messages
        )
        if mark is not None:
            resolutions[file] = mark

    with open(resolution_file, 'w') as f:
        json.dump(resolutions, f, indent=2)


def postcheck(hg_ggg: Path, git_ggg: Path, report_file: Path, final_messages: list):
    """Perform a final comparison of the files in the two repos.

    Creates a report on files that are only present in one repo, or which have different
    checksums between the repos."""
    hg_files = _get_all_ggg_files(hg_ggg)
    git_files = _get_all_ggg_files(git_ggg)

    missing_from_git = sorted(set(hg_files).difference(git_files))
    missing_from_hg = sorted(set(git_files).difference(hg_files))
    common_files = set(hg_files).intersection(git_files)
    different = []

    for relative_file in common_files:
        hg_file = hg_ggg / relative_file
        git_file = git_ggg / relative_file
        if not _compare_file_md5(hg_file, git_file):
            different.append(relative_file)

    with open(report_file, 'w') as rpt:
        if different:
            rpt.write('\n')
            rpt.write(_banner('Files with difference checksums'))
            rpt.write('\n\n')
            for relative_file in different:
                rpt.write(f'{relative_file} has different checksums in the two repos\n')

        if missing_from_git:
            rpt.write('\n')
            rpt.write(_banner('Files missing in git repo'))
            rpt.write('\n\n')
            for relative_file in missing_from_git:
                rpt.write(f'{relative_file} present in the Mercurial repo but missing from the Git repo\n')

        if missing_from_hg:
            rpt.write('\n')
            rpt.write(_banner('Files not from hg repo'))
            rpt.write('\n\n')
            for relative_file in missing_from_hg:
                rpt.write(f'{relative_file} present in the Git repo but missing from the Mercurial repo\n')

    n_missing_from_git = len(missing_from_git)
    n_missing_from_hg = len(missing_from_hg)
    n_different = len(different)

    if n_missing_from_hg > 0:
        final_messages.append(
            f'{n_missing_from_hg} files were only in the Git repo, confirm that this is expected (report is at {report_file})'
        )
    if n_missing_from_git > 0:
        final_messages.append(
            f'{n_missing_from_git} files were missing from the Git repo, check the report ({report_file})'
        )
    if n_different > 0:
        final_messages.append(
            f'{n_different} files were different between the repos, check the report ({report_file})'
        )


def _git_is_clean(git_ggg: Path) -> bool:
    """Return True if a Git repo has no uncommitted changes or untracked files."""
    # git status --porcelain returns no lines if there are no outstanding changes.
    out = run(['git', 'status', '--porcelain'], cwd=git_ggg, capture_output=True)
    return len(out.stdout) == 0


def _hg_is_clean(hg_ggg: Path) -> bool:
    """Return True if a Mercurial repo has no uncommitted changes or untracked files."""
    # hg status should just print nothing if no changes, though this isn't as guaranteed
    # as git status --porcelain.
    out = run(['hg', 'status'], cwd=hg_ggg, capture_output=True)
    return len(out.stdout) == 0


def _iter_patch_blocks(f):
    """Iterate over the blocks for each file in a patch"""
    block_lines = []
    for line in f:
        if line.startswith('diff') and len(block_lines) > 0:
            yield block_lines
            block_lines = [line]
        else:
            block_lines.append(line)
    yield block_lines


def _get_block_files(block_lines):
    """Get the old and new files from a patch block.

    Note that this reads the first line, not the +++ and --- lines, because
    the latter don't exist if the file is empty. The first line does not
    use /dev/null if the file was created/deleted."""
    # This won't work if file names happen to have " b/" in them,
    # but I don't see a way to distinguish that case. We can't use
    # the --- and +++ lines, because those don't get shown if the
    # file is empty.
    #
    # Note that unlike the --- and +++ lines, there will not be a /dev/null
    # entry if the file is created or deleted.
    s = block_lines[0].replace('diff --git a/', '').strip()
    f1, f2 = s.split(' b/', maxsplit=1)
    return f1, f2



def _replace_hgignore(block_lines):
    """Replace .hgignore in the header of a diff block"""
    for i in range(len(block_lines)):
        if block_lines[i].startswith('@'):
            # done with the "header", we don't want to change the actual diff
            # since Mercurial doesn't make it as easy to have .hgignore files
            # in subdirectories, we assume that .hgignore should only be listed
            # in the file paths
            return block_lines

        block_lines[i] = block_lines[i].replace('.hgignore', '.gitignore')

    return block_lines


def _replace_dir_placeholder(block_lines):
    """Modify a block for a directory placeholder file to make it into a .gitignore file"""

    if len(block_lines) == 2 and block_lines[1].startswith('new file'):
        # We created a new, empty dir placeholder file. Change the path to be
        # a git ignore file instead, and add content to the file.
        for i in range(len(block_lines)):
            block_lines[i] = block_lines[i].replace(DIR_PLACEHOLDER_FILE, '.gitignore')
        block_lines.extend([
            '--- /dev/null\n',
            '+++ b/.gitignore\n',
            '@@ -0,0 +1,2 @@\n',
            '+*\n',
            '+!.gitignore\n'
        ])
    elif block_lines[1].startswith('deleted file'):
        # If we did want to delete the file, we'd have to get its contents for the unified diff.
        # An empty dir placeholder file wouldn't have that. But we may not want to delete the .gitignore
        # file, so I'm leaving this until I have an actual example.
        raise NotImplementedError('Deleted directory placeholder')
    else:
        # Assume that this was a rename
        for i in range(len(block_lines)):
            block_lines[i] = block_lines[i].replace(DIR_PLACEHOLDER_FILE, '.gitignore')

    return block_lines


def _is_linelist_file(file):
    """Return True if the relative path ``file`` is likely to be a linelist or related file.

    That is, this file should not be moved to the Git repo; linelists should go into CaltechData."""
    file = Path(file)
    if not file.parts[0] == 'linelist':
        return False

    if file.suffix in {'.py', '.json'}:
        return False

    return True


def _check_linelist_file(file, hg_ggg: Path, final_messages: list):
    """Check if the linelist file ``file`` is in the JSON of linelist files to download and it has the right checksum.

    If not, this is an indication that that JSON file needs updated."""
    file_key = Path(file).relative_to('linelist').as_posix()
    linelist_json = hg_ggg / 'linelist' / '.linelist_info.json'
    with open(linelist_json) as f:
        linelist_info = json.load(f)
    for src_group in linelist_info:
        for src_file in src_group['contents']:
            if src_file['file'] == file_key:
                if not _check_md5(file, src_file['md5']):
                    final_messages.append(f'Linelist file {file} has changed, update the linelist data repository and JSON file')
                else:
                    return

    final_messages.append(f'New linelist file {file} found, update the linelist data repository and JSON file')


def _parse_git_apply_check(stderr: str) -> dict:
    """Parse the result of ``git apply --check`` into a dictionary mapping failure reasons to file paths.
    """
    failures = {WHY_FAILED_NO_FILE: set(), WHY_FAILED_EXISTS: set(), WHY_FAILED_CANT_PATCH: set()}
    for line in stderr.splitlines():
        parts = line.split(':')
        _, file_or_note, file_or_reason = parts[:3]
        if file_or_note.strip() == 'patch failed':
            # This should be a line like "error: patch failed: install/check_python.sh:15",
            # get the file name from the third part, don't include the line number.
            # This means that the patch couldn't be matched up with the existing file, so
            # we'll want to possibly copy over the new version manually.
            file = file_or_reason.split(':')[0].strip()
            failures[WHY_FAILED_CANT_PATCH].add(file)
        elif file_or_reason.strip() == 'patch does not apply':
            # This should be a line like "error: install/check_python.sh: patch does not apply"
            # Usually this comes after the "patch failed" line, but catch both just in case.
            # Same issue as the "patch failed" branch above.
            file = file_or_note.strip()
            failures[WHY_FAILED_CANT_PATCH].add(file)
        elif file_or_reason.strip() == 'already exists in working directory':
            # This usually means that the file was previously added in Git and the patch
            # is trying to add it again. We'll want to compare them, and see if we need
            # to manually copy the Mercurial version.
            file = file_or_note.strip()
            failures[WHY_FAILED_EXISTS].add(file)
        elif file_or_reason.strip() == 'No such file or directory':
            # This means that Mercurial only had edits to the file but Git didn't have it
            # at all. Again, we'll want to look at the Mercurial file to see if we need
            # to copy it.
            file = file_or_note.strip()
            failures[WHY_FAILED_NO_FILE].add(file)
        else:
            raise NotImplementedError(f'git apply --check output line: {line}')

    return failures


def _compare_and_maybe_copy_failed_file(
    failed_file: str, hg_ggg: Path, git_ggg: Path,
    final_messages: list, reject_file=None
):
    """Handle comparing a single file that could not be patched between the two repos.

    Asks the user how to handle it.
    """
    actions = []

    hg_file = hg_ggg / failed_file
    git_file = git_ggg / failed_file
    files_exist = hg_file.exists() and git_file.exists()
    files_are_text = _is_text_file(hg_file) and _is_text_file(git_file)
    if files_exist and files_are_text:
        actions.append(ACT_DIFF)

    if reject_file is not None and reject_file.exists():
        actions.append(ACT_VIEW_REJ)

    if hg_file.exists() and _is_text_file(hg_file):
        actions.append(ACT_VIEW_HG)

    if git_file.exists() and _is_text_file(git_file):
        actions.append(ACT_VIEW_GIT)

    actions.append(ACT_COPY_HG_TO_GIT)
    actions.append(ACT_IS_OK)
    actions.append(ACT_NOTHING)

    while True:
        idx = _get_user_selection('What do you want to do?', actions)
        if actions[idx] == ACT_DIFF:
            _show_diff(hg_file, git_file)
        elif actions[idx] == ACT_VIEW_REJ:
            _show_file(reject_file)
        elif actions[idx] == ACT_VIEW_HG:
            _show_file(hg_file)
        elif actions[idx] == ACT_VIEW_GIT:
            _show_file(git_file)
        elif actions[idx] == ACT_COPY_HG_TO_GIT:
            shutil.copy(str(hg_file), str(git_file))
            print(f'Copied {hg_file} to {git_file}')
            if reject_file is not None:
                reject_file.unlink()
            return MARK_COPY_HG_TO_GIT
        elif actions[idx] == ACT_IS_OK:
            print('Left this file alone')
            if reject_file is not None:
                reject_file.unlink()
            return MARK_IS_OK
        elif actions[idx] == ACT_NOTHING:
            print('Doing nothing')
            final_messages.append(f'A file that could not be patched may still need attention: {failed_file}')
            return None


def _get_all_ggg_files(
    root,
    skip_prefixes=('.hg', '.git', '.conda'),
    skip_suffixes=('.patch', '.rpt'),
    skip_dirs=('current_results', 'bin'),
    skip_files=(DIR_PLACEHOLDER_FILE,)
):
    """Get all files under a GGG directory with some exceptions.

    By default, we skip:

    - subdirectories or files starting with ".hg", ".git", or ".conda" as these are
      VCS directories or Python environments.
    - Files ending in ".patch" or ".rpt" as these are usually written by this program,
      and other .rpt files are runtime artifacts.
    - The directories "current_results" and "bin", since these always contain benchmark
      test output and compiled programs, respectively, and would give false positives.
    - The ".dir_placeholder_for_hg" files, these will only exist on the Mercurial side,
      so would also give false positives when looking for missed files.
    """
    all_files = []
    for path, dirnames, filenames in os.walk(root):
        dirs_to_remove = [d for d in dirnames if d.startswith(skip_prefixes) or d in skip_dirs]
        for d in dirs_to_remove:
            dirnames.remove(d)
        path = Path(path).relative_to(root)
        files = [path / f for f in filenames if not f.startswith(skip_prefixes) and f not in skip_files and Path(f).suffix not in skip_suffixes]
        all_files.extend(files)
    return all_files


def _compare_file_md5(file1, file2):
    """Return ``True`` if two files have the same MD5 checksum"""
    sum1 = _compute_md5(file1)
    sum2 = _compute_md5(file2)
    return sum1 == sum2


def _check_md5(file, md5):
    """Return ``True`` if ``file`` has an MD5 checksum equal to ``md5``"""
    if not Path(file).is_file():
        return False
    new_checksum = _compute_md5(file)
    return new_checksum == md5


def _compute_md5(file):
    """Compute the MD5 checksum of ``file``"""
    CHUNK_SIZE = 1000000
    checksum = hashlib.md5()
    with open(file, 'rb') as f:
        chunk = f.read(CHUNK_SIZE)
        while len(chunk) > 0:
            checksum.update(chunk)
            chunk = f.read(CHUNK_SIZE)
    return checksum.hexdigest()


def _is_text_file(file, nbytes=1024):
    """Guess if ``file`` is text or binary

    Returns ``True`` if exists and is text, ``False`` otherwise.
    By default, a file is considered text if the first 1024 bytes
    are all valid UTF-8."""
    if not Path(file).exists():
        return False

    with open(file) as f:
        try:
            f.read(nbytes)
        except UnicodeDecodeError:
            return False
        else:
            return True


def _get_user_selection(prompt, options):
    """Prompt a user to select from ``options``"""
    print(prompt)
    for iopt, opt in enumerate(options, start=1):
        print(f'{iopt}) {opt}')
    while True:
        input_str = input(f'Enter 1-{len(options)}: ')
        try:
            user_index = int(input_str)
        except ValueError:
            print(f'Sorry, "{input_str}" is not a valid input')
        else:
            if 1 <= user_index <= len(options):
                return user_index - 1
            else:
                print(f'Sorry, "{input_str}" is not one of the allowed options')


def _show_diff(file1, file2):
    """Open the difference between two files in an external program"""
    diffprog = os.getenv('DIFFPROG', 'vimdiff')
    run([diffprog, str(file1), str(file2)])


def _show_file(file):
    """Open a file in an external editor, or print the file to the screen directly."""
    editor = os.getenv('EDITOR', None)
    if editor is None:
        with open(file) as f:
            for line in f:
                print(line.rstrip())
    else:
        run([editor, str(file)])


def _banner(title, c='='):
    """Return ``title`` wrapped in the character ``c`` to act as a banner"""
    n = len(title) + 4
    row = c * n
    return f'{row}\n{c} {title} {c}\n{row}'


if __name__ == '__main__':
    main()
