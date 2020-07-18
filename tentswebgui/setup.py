from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="tentswebgui",
    version="0.0.1",
    author="Example Author",
    author_email="author@example.com",
    description="Visualization of spacetime tents",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jayggg/ngstents",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    packages=find_packages("src"),
    # packages=find_packages("tentswebgui"),
    package_dir={"": "src"},
    package_data={
        'tentswebgui': ['labextension/plugin.*', 'nbextension/static/*.*'],
    },
)

