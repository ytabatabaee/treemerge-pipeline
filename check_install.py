#!/usr/bin/env python3

import importlib
import contextlib
import io
import os
import platform
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
SOFTWARE_DIR = REPO_ROOT / "software"

PYTHON_MODULES = [
    "dendropy",
    "ete3",
    "networkx",
    "numpy",
    "scipy",
]

LINUX_X86_64_TOOLS = {
    "ASTRID": SOFTWARE_DIR / "ASTRID-linux",
    "FastME": SOFTWARE_DIR / "fastme-2.1.5-linux64",
    "ASTRAL4": SOFTWARE_DIR / "astral4-linux",
    "TreeMerge": SOFTWARE_DIR / "treemerge.py",
    "PAUP*": [
        SOFTWARE_DIR / "paup4a169_ubuntu64",
        SOFTWARE_DIR / "paup4a168_centos64",
    ],
}

DARWIN_COMMON_TOOLS = {
    "ASTRID": SOFTWARE_DIR / "ASTRID-osx",
    "FastME": SOFTWARE_DIR / "fastme-2.1.5-osx",
    "TreeMerge": SOFTWARE_DIR / "treemerge.py",
    "PAUP*": SOFTWARE_DIR / "paup4a168_osx",
}


def status(ok, message):
    marker = "OK" if ok else "FAIL"
    print(f"[{marker}] {message}")


def paup_smoke_test(path):
    with tempfile.NamedTemporaryFile(
        "w",
        suffix=".nex",
        prefix="treemerge-paup-check-",
        delete=False,
    ) as handle:
        handle.write("#NEXUS\nbegin paup;\nq;\nend;\n")
        check_file = handle.name

    try:
        result = subprocess.run(
            [str(path), "-n", check_file],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=15,
        )
        return result.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        return False
    finally:
        try:
            os.unlink(check_file)
        except OSError:
            pass


def bundled_tools_for_platform():
    system = platform.system()
    machine = platform.machine().lower()

    if system == "Linux" and machine in {"x86_64", "amd64"}:
        return LINUX_X86_64_TOOLS, []

    if system == "Darwin":
        tools = dict(DARWIN_COMMON_TOOLS)
        warnings = [
            "Bundled macOS binaries have mixed x86_64/arm64 architectures.",
            "Apple Silicon users may need Rosetta or explicit native tool paths.",
        ]
        if machine == "arm64":
            tools["ASTRAL4"] = SOFTWARE_DIR / "astral4-osx"
        else:
            warnings.append(
                "No bundled Intel macOS ASTRAL4 binary is present; pass "
                "--astral4-bin manually or run on Linux/Apple Silicon."
            )
        return tools, warnings

    return {}, [
        f"No bundled binary set is defined for {system or 'unknown'} "
        f"{platform.machine() or 'unknown'}."
    ]


def check_python_modules():
    failures = []
    print("Python modules")

    for module in PYTHON_MODULES:
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                with contextlib.redirect_stdout(io.StringIO()):
                    importlib.import_module(module)
            status(True, module)
        except Exception as exc:
            failures.append(module)
            status(False, f"{module}: {exc}")

    return failures


def check_paup_candidates(paths):
    failures = []
    found = False

    for path in paths:
        if not path.exists():
            status(False, f"PAUP*: missing {path}")
            continue

        if not os.access(path, os.X_OK):
            status(False, f"PAUP*: not executable {path}")
            continue

        if paup_smoke_test(path):
            status(True, f"PAUP*: selected {path.relative_to(REPO_ROOT)}")
            found = True
            break

        status(False, f"PAUP*: smoke test failed {path.relative_to(REPO_ROOT)}")

    if not found:
        failures.append("PAUP*")

    return failures


def check_tools():
    failures = []
    tools, warnings = bundled_tools_for_platform()

    print("\nBundled external tools")
    print(f"Platform: {platform.system()} {platform.machine()}")

    for warning in warnings:
        print(f"[WARN] {warning}")

    if not tools:
        failures.append("bundled tools")
        return failures

    for name, path in tools.items():
        if name == "PAUP*" and isinstance(path, list):
            failures.extend(check_paup_candidates(path))
            continue

        if not path.exists():
            failures.append(name)
            status(False, f"{name}: missing {path}")
            continue

        if name != "TreeMerge" and not os.access(path, os.X_OK):
            failures.append(name)
            status(False, f"{name}: not executable {path}")
            continue

        status(True, f"{name}: {path.relative_to(REPO_ROOT)}")

    return failures


def smoke_test_treemerge_import():
    print("\nTreeMerge import smoke test")
    cmd = [
        sys.executable,
        "-c",
        "import sys; sys.path.insert(0, 'software'); import treemerge",
    ]
    result = subprocess.run(
        cmd,
        cwd=REPO_ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode == 0:
        status(True, "software/treemerge.py imports successfully")
        return []

    message = result.stderr.strip() or result.stdout.strip()
    status(False, f"software/treemerge.py import failed: {message}")
    return ["treemerge import"]


def main():
    failures = []
    failures.extend(check_python_modules())
    failures.extend(check_tools())
    failures.extend(smoke_test_treemerge_import())

    if failures:
        print("\nInstall check failed.")
        print("Missing or broken components:")
        for failure in failures:
            print(f"  - {failure}")
        print("\nRecommended fix:")
        print("  python3 -m venv .venv")
        print("  source .venv/bin/activate")
        print("  python -m pip install -r requirements.txt")
        sys.exit(1)

    print("\nInstall check passed.")


if __name__ == "__main__":
    main()
