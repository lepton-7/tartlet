from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read()

setup(
    name="bio-tart",
    version="0.1.0",
    author="Sachit Kshatriya",
    author_email="sxk1464@case.edu",
    license="MIT",
    description="Trace riboswitches modulating transcription-termination.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lepton-7/tart",
    py_modules=["targeted_pipeline"],
    packages=find_packages(),
    install_requires=[requirements],
    python_requires=">=3.9,<3.12",
    classifiers=[
        "Programming Language :: Python :: 3.11",
        "Operating System :: Linux",
    ],
    entry_points={
        "console_scripts": [
            "tart-targeted = targeted_pipeline:cli",
        ],
    },
)
