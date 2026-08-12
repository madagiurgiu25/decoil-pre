#!/usr/bin/env python3

import importlib.metadata
import json
import os
import platform
import subprocess
import sys

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.version import Version, InvalidVersion


# ============================================================
# Configuration
# ============================================================

DISTRIBUTION_NAME = "decoil"


# Conda dependencies required by DECOIL.
from pathlib import Path

def get_conda_requirements_file():
    project_root = Path(__file__).resolve().parent.parent
    system = platform.system()

    if system == "Darwin":
        return project_root / "conda-requirements-osx.txt"

    if system == "Linux":
        return project_root / "conda-requirements-linux.txt"

    if system == "Windows":
        print("WARNING: Windows support has not been tested.")
        return None

    print(f"WARNING: Unsupported operating system: {system}")
    return None

# ============================================================
# Output
# ============================================================

def print_result(name, required, installed, passed):
    status = "PASSED" if passed else "FAILED"

    print(
        f"{name:<20}"
        f"required: {required:<15}"
        f"installed: {installed:<15}"
        f"{status}"
    )


def print_section(title):
    print()
    print(title)
    print("-" * len(title))


# ============================================================
# Python dependencies
# ============================================================

def get_python_requirements():
    """
    Get the dependencies recorded for the installed DECOIL package.

    These originate from install_requires in setup.py.

    setup.py is NOT imported or executed.
    """

    try:
        requirements = importlib.metadata.requires(
            DISTRIBUTION_NAME
        )

    except importlib.metadata.PackageNotFoundError:
        print(
            f"ERROR: Python distribution "
            f"'{DISTRIBUTION_NAME}' is not installed."
        )
        return None

    return requirements or []

def check_python_dependencies():
    """
    Check Python dependencies from the installed DECOIL metadata.
    """

    print_section("Python dependencies")

    requirements = get_python_requirements()

    if requirements is None:
        return False

    if not requirements:
        print("No Python dependencies declared.")
        return True

    success = True

    for requirement_string in requirements:

        try:
            requirement = Requirement(requirement_string)

        except Exception as exc:
            print_result(
                requirement_string,
                "invalid",
                str(exc),
                False,
            )

            success = False
            continue

        # Example:
        #
        # importlib-metadata>=6;
        # python_version < "3.10"
        #
        # If the marker does not apply to the current Python,
        # the requirement should not be checked.
        if (
            requirement.marker is not None
            and not requirement.marker.evaluate()
        ):
            continue

        name = requirement.name

        required = (
            str(requirement.specifier)
            if requirement.specifier
            else "any"
        )

        try:
            installed_version = importlib.metadata.version(name)

        except importlib.metadata.PackageNotFoundError:
            print_result(
                name,
                required,
                "NOT INSTALLED",
                False,
            )

            success = False
            continue

        if requirement.specifier:
            passed = installed_version in requirement.specifier
        else:
            passed = True

        print_result(
            name,
            required,
            installed_version,
            passed,
        )

        if not passed:
            success = False

    return success


# ============================================================
# Conda dependencies
# ============================================================

def parse_conda_requirement(requirement):
    """
    Parse Conda requirements.

    Examples:

        bedtools=2.31.1
        samtools>=1.20
        sniffles=1.0.7
        bwa

    Returns:

        ("bedtools", "==2.31.1")
        ("samtools", ">=1.20")
        ("bwa", "")
    """

    requirement = requirement.strip()

    operators = [
        ">=",
        "<=",
        "!=",
        "==",
        ">",
        "<",
        "=",
    ]

    for operator in operators:

        if operator in requirement:

            name, version = requirement.split(
                operator,
                1,
            )

            name = name.strip()
            version = version.strip()

            if not name or not version:
                raise ValueError(
                    f"Invalid Conda requirement: "
                    f"{requirement}"
                )

            # Conda usually uses:
            #
            # package=1.2.3
            #
            # packaging expects:
            #
            # package==1.2.3
            if operator == "=":
                operator = "=="

            return name, operator + version

    return requirement, ""


def read_conda_requirements(requirements_file):
    """
    Read Conda requirements from a text file.

    Supported examples:
        bedtools=2.31.1
        samtools>=1.7
        sniffles=1.0.7

    Blank lines and comments starting with # are ignored.
    Inline comments are also supported.
    """
    requirements = []

    if requirements_file is None:
        return requirements

    if not requirements_file.exists():
        raise FileNotFoundError(
            f"Conda requirements file not found: {requirements_file}"
        )

    with requirements_file.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            if " #" in line:
                line = line.split(" #", 1)[0].strip()

            if line:
                requirements.append(line)

    return requirements


def get_conda_packages():
    """
    Return packages installed in the currently active Conda
    environment.

    Example result:

        {
            "bedtools": "2.31.1",
            "samtools": "1.21",
            ...
        }
    """

    try:
        result = subprocess.run(
            [
                "conda",
                "list",
                "--json",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    except FileNotFoundError:
        return None, "conda command not found"

    except subprocess.CalledProcessError as exc:

        message = exc.stderr.strip()

        if not message:
            message = "conda list failed"

        return None, message

    try:
        data = json.loads(result.stdout)

    except json.JSONDecodeError:
        return None, "invalid output from 'conda list --json'"

    packages = {}

    for package in data:

        name = package.get("name")
        version = package.get("version")

        if name and version:
            packages[name.lower()] = version

    return packages, None


def conda_version_matches(installed, specifier):
    """
    Check whether a Conda package version satisfies a version
    requirement.

    Most Conda versions can be compared using packaging.Version.
    """

    if not specifier:
        return True

    try:
        version = Version(installed)
        spec = SpecifierSet(specifier)

        return version in spec

    except InvalidVersion:

        # Fallback for non-PEP440 Conda versions.
        if specifier.startswith("=="):
            expected = specifier[2:]
            return installed == expected

        return False


def check_conda_dependencies():
    """
    Check Conda dependencies from the OS-specific requirements file.
    """

    print_section("Conda dependencies (relevant only for decoil-pipeline)")

    requirements_file = get_conda_requirements_file()

    if requirements_file is None:
        print("No Conda requirements configured for this platform.")
        return True

    try:
        requirements = read_conda_requirements(requirements_file)
    except FileNotFoundError as exc:
        print(f"FAILED: {exc}")
        return False

    print(f"Requirements file: {requirements_file.name}")

    if not requirements:
        print("No Conda dependencies declared.")
        return True

    packages, error = get_conda_packages()

    if packages is None:
        print(f"FAILED: {error}")
        return False

    success = True

    for requirement_string in requirements:
        try:
            name, specifier = parse_conda_requirement(requirement_string)
        except ValueError as exc:
            print(f"FAILED: {exc}")
            success = False
            continue

        installed = packages.get(name.lower())
        required = specifier or "any"

        if installed is None:
            print_result(name, required, "NOT INSTALLED", False)
            success = False
            continue

        passed = conda_version_matches(installed, specifier)
        print_result(name, required, installed, passed)

        if not passed:
            success = False

    return success


# ============================================================
# Environment
# ============================================================

def print_environment():
    print()
    print("DECOIL dependency check")
    print("=======================")

    print()
    print(f"Platform:          {platform.system()}")
    print(f"Architecture:      {platform.machine()}")
    print(f"Python:            {platform.python_version()}")

    conda_environment = os.environ.get(
        "CONDA_DEFAULT_ENV"
    )

    if conda_environment:
        print(
            f"Conda environment: {conda_environment}"
        )
    else:
        print(
            "Conda environment: not detected"
        )


# ============================================================
# Main checker
# ============================================================

def run_dependency_check():
    """
    Run all DECOIL dependency checks.

    Returns:
        0 if everything passed
        1 if something failed
    """

    print_environment()

    python_ok = check_python_dependencies()

    conda_ok = check_conda_dependencies()

    print()
    print("=" * 70)

    if python_ok and conda_ok:
        print("Overall result: PASSED")
        return 0

    print("Overall result: FAILED")
    return 1


# ============================================================
# Standalone execution
# ============================================================

if __name__ == "__main__":
    sys.exit(run_dependency_check())