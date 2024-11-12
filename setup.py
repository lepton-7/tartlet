from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read()

setup(
    name="tartlet",
    version="0.7.9",
    author="Sachit Kshatriya",
    author_email="sxk1464@case.edu",
    license="MIT",
    description="Trace riboswitches modulating transcription-termination.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lepton-7/tartlet",
    py_modules=[],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[requirements],
    python_requires=">=3.11",
    classifiers=[
        "Programming Language :: Python :: 3.11",
        "Operating System :: Linux",
    ],
    entry_points={
        "console_scripts": [
            "tart-targeted = tartlet.entry_points:targeted",
            "tart-utils = tartlet.entry_points:utils",
        ],
    },
)
