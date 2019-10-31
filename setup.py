#  Developed by Lukas Sauer at the Heidelberg Institute of Theoretical Studies on 10/30/19, 5:04 PM.
#  Contact: lukas.sauer@h-its.org.
#  Last modified on 10/30/19, 5:04 PM.
#  (C) 2019. All rights reserved.

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="flagcomplex",
    version="0.1.0",
    author="Lukas D. Sauer",
    # author_email="author@example.com",
    description="This is a project to visualize tuples of flags in the real projective plane."
                " It was developed at Heidelberg Institute for Theoretical Studies in 2019.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LukasDSauer/flagcomplex",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "numpy>=1.17.3",
        "drawSvg>=1.0.2"
        "scipy>=1.1.0"
    ]
)
