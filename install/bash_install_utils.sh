do_use_micromamba() {
    [ ! -z "$GGG_USE_MICROMAMBA" ] && [ "$GGG_USE_MICROMAMBA" -gt 0 ] 2>/dev/null
}

do_use_pip_only() {
    [ ! -z "$GGG_USE_PIP" ] && [ "$GGG_USE_PIP" -gt 0 ] 2>/dev/null
}

pytool_cfg_conflict() {
    if do_use_micromamba && do_use_pip_only; then
        echo "Both GGG_USE_MICROMAMBA and GGG_USE_PIP set to non-zero values."
        return 0
    else
        return 1
    fi
}

version_ge() {
    # Returns 0 (true) if $1 >= $2, otherwise 1 (false)
    test "$(printf '%s\n' "$@" | sort -V | head -n 1)" == "$2"
}
