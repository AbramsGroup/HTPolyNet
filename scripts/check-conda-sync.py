#!/usr/bin/env python3
# Compare pyproject.toml's runtime dependencies against the
# conda-forge feedstock's `requirements.run` list and report drift.
#
# Why: the conda-forge regro-cf-autotick-bot can only bump version +
# sha256 on a release.  If we add a new runtime dependency in
# pyproject.toml (or rename / remove one), the bot won't notice; the
# next conda-forge release will install but crash at first use.
# This check catches that before tagging — see the release-time
# integration in scripts/release.sh.
#
# It compares two things: the set of package NAMES, and each one's
# lower version bound.  The floor check was added after grayskull, on
# feedstock PR #21, spotted that the recipe carried a bare
# `matplotlib-base` against pyproject's `matplotlib>=3.5` while this
# script reported "In sync" — because it compared names with the
# version specifiers stripped off, and so could not see a floor at all.
# A recipe floor at or above pyproject's is fine; only a missing or
# lower one is drift, since that is what lets conda solve an
# environment pip would have refused.
#
# Usage:
#   ./scripts/check-conda-sync.py             # informational
#   ./scripts/check-conda-sync.py --strict    # exit 1 on drift
#
# Exit codes:
#   0  in sync (or warn-only mode)
#   1  drift detected and --strict was set
#   2  unexpected error (network / parse failure)
import argparse
import re
import sys
import tomllib
import urllib.request
from pathlib import Path

FEEDSTOCK_RAW = (
    'https://raw.githubusercontent.com/conda-forge/'
    'htpolynet-feedstock/main/recipe/meta.yaml'
)

# PyPI → conda-forge package-name remap.  Lower-cased; left side is
# the pyproject.toml dependency name, right side is the conda-forge
# package name.  Add to this when a dep is repackaged under a
# different name on conda-forge.
PYPI_TO_CONDA = {
    'graphviz': 'python-graphviz',  # python wrapper; system 'dot' is `graphviz`
    'matplotlib': 'matplotlib-base',  # conda-forge convention: -base = no GUI
}

# conda-forge-only run deps that have no PyPI equivalent (or that we
# pull from conda-forge specifically to get the system binary rather
# than a Python wheel).  Listing them here suppresses "extra in
# conda recipe" warnings.
CONDA_ONLY = {
    'python',          # always pinned in conda recipes
    'ambertools',      # native MD binaries; conda-forge only
    'graphviz',        # system 'dot' binary (we list both this and python-graphviz)
}


def parse_pyproject_deps(path: Path) -> dict:
    """Return {PyPI package name: lower bound or None} from
    `[project] dependencies`, lower-cased."""
    with path.open('rb') as f:
        data = tomllib.load(f)
    deps = data['project']['dependencies']
    return {_canon(d): _floor(d) for d in deps}


def parse_recipe_run_deps(text: str) -> dict:
    """Return {conda-forge package name: lower bound or None} from a
    meta.yaml's `requirements.run:` list block."""
    in_requirements = False
    in_run = False
    indent = None
    pkgs = {}
    for line in text.splitlines():
        stripped = line.strip()
        if stripped.startswith('requirements:'):
            in_requirements = True
            continue
        if in_requirements and stripped.startswith('run:'):
            in_run = True
            indent = len(line) - len(line.lstrip())
            continue
        if in_run:
            if not stripped or stripped.startswith('#'):
                continue
            line_indent = len(line) - len(line.lstrip())
            if line_indent <= indent and not stripped.startswith('-'):
                in_run = False
                in_requirements = stripped.endswith(':')
                continue
            if stripped.startswith('-'):
                spec = stripped.lstrip('- ')
                if spec.split():
                    pkgs[_canon(spec)] = _floor(spec)
    return pkgs


def _floor(spec: str) -> str | None:
    """Return the `>=` lower bound in a dependency spec, or None.

    'numpy>=1.24' -> '1.24';  'matplotlib-base' -> None;
    'pandas >=2,<3' -> '2'.  Only the lower bound matters here: an
    upper bound in the recipe that pyproject lacks is a deliberate
    conda-side pin, not drift.
    """
    m = re.search(r'>=\s*([0-9][0-9A-Za-z.\-_]*)', spec)
    return m.group(1) if m else None


def _version_key(v: str) -> tuple:
    """Comparable key for a version string, numeric segment by segment.

    Conda and pip agree that 2024.03 and 2024.3 are the same version;
    a plain string compare does not, which would have made the rdkit
    pin look like drift forever.
    """
    parts = []
    for p in re.split(r'[.\-_]', v):
        parts.append((0, int(p), '') if p.isdigit() else (1, 0, p))
    return tuple(parts)


def _floor_is_weaker(recipe_floor: str | None, required_floor: str) -> bool:
    """True if the recipe's lower bound is absent or below what pyproject needs."""
    if recipe_floor is None:
        return True
    a, b = _version_key(recipe_floor), _version_key(required_floor)
    pad = max(len(a), len(b))
    a += ((0, 0, ''),) * (pad - len(a))
    b += ((0, 0, ''),) * (pad - len(b))
    return a < b


def _canon(spec: str) -> str:
    """Strip version specifiers and lower-case a dep name.
    'numpy>=1.24' -> 'numpy';  'python-graphviz >=0.20' -> 'python-graphviz'."""
    name = re.split(r'[<>=!~ ]', spec, maxsplit=1)[0]
    return name.strip().lower()


def map_to_conda(pypi: dict) -> dict:
    return {PYPI_TO_CONDA.get(n, n): floor for n, floor in pypi.items()}


def compare(pypi: dict, conda: dict) -> tuple[set, set, list]:
    """Return (missing names, unexpected names, weakened floors).

    A weakened floor is a dependency the recipe has, but with a lower
    bound that is absent or below what pyproject.toml requires -- the
    case where conda can solve an environment pip would have refused.
    """
    expected = map_to_conda(pypi)
    missing = set(expected) - set(conda)
    extra = set(conda) - set(expected) - CONDA_ONLY
    weaker = []
    for name, required in sorted(expected.items()):
        if required is None or name not in conda:
            continue
        if _floor_is_weaker(conda[name], required):
            weaker.append((name, conda[name], required))
    return missing, extra, weaker


def fetch_recipe(url: str = FEEDSTOCK_RAW) -> str:
    with urllib.request.urlopen(url, timeout=15) as resp:
        return resp.read().decode('utf-8')


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=("Compare pyproject.toml's runtime dependencies "
                     "against the conda-forge feedstock's run-deps."))
    p.add_argument('--strict', action='store_true',
                   help='exit 1 if drift is found (default: warn only)')
    p.add_argument('--pyproject', type=Path, default=Path('pyproject.toml'),
                   help='path to local pyproject.toml')
    p.add_argument('--recipe-url', default=FEEDSTOCK_RAW,
                   help='URL of the conda-forge recipe to compare against')
    args = p.parse_args(argv)

    try:
        pypi_deps = parse_pyproject_deps(args.pyproject)
    except Exception as e:
        print(f'error: failed to parse {args.pyproject}: {e}', file=sys.stderr)
        return 2

    try:
        recipe_text = fetch_recipe(args.recipe_url)
        conda_deps = parse_recipe_run_deps(recipe_text)
    except Exception as e:
        print(f'error: failed to fetch/parse recipe at {args.recipe_url}: {e}',
              file=sys.stderr)
        return 2

    if not conda_deps:
        print(f'error: no run deps parsed from {args.recipe_url}; recipe '
              f'format may have changed', file=sys.stderr)
        return 2

    missing, extra, weaker = compare(pypi_deps, conda_deps)

    print(f'pyproject.toml deps:    {len(pypi_deps)}')
    print(f'conda-forge run deps:   {len(conda_deps)}')
    print()

    if missing:
        print('MISSING from conda recipe (in pyproject, not in recipe):')
        for n in sorted(missing):
            print(f'  - {n}')
        print()

    if extra:
        print('UNEXPECTED in conda recipe (not in pyproject, not in '
              'known conda-only allowlist):')
        for n in sorted(extra):
            print(f'  - {n}')
        print()

    if weaker:
        print('WEAKER version floor in conda recipe than pyproject requires:')
        for name, have, need in weaker:
            shown = have if have is not None else '(no lower bound)'
            print(f'  - {name}: recipe has {shown}, pyproject needs >={need}')
        print()

    if not missing and not extra and not weaker:
        print('In sync.  (conda-forge feedstock recipe matches pyproject.toml '
              "dependencies and version floors; the autotick bot's next version-only PR can "
              'auto-merge.)')
        return 0

    print('Action: when you next bump htpolynet, the autotick-bot PR on the '
          'conda-forge feedstock will be missing these recipe edits.  '
          'Either prepare the recipe update before pushing the release tag, '
          'or plan to supersede the bot PR with a manual update.')

    return 1 if args.strict else 0


if __name__ == '__main__':
    sys.exit(main())
