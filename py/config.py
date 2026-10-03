"""Paths to the native ngstents development package."""

from pathlib import Path


def _package_dir():
    return Path(__file__).resolve().parent


def get_cmake_dir():
    """Return the directory containing ngstentsConfig.cmake."""
    return str(_package_dir() / "cmake")


def get_include_dir():
    """Return the directory containing the installed public headers."""
    return str(_package_dir() / "include")


def get_library_dir():
    """Return the directory containing the ngstents core shared library."""
    return str(_package_dir())


if __name__ == "__main__":
    print(get_cmake_dir())
