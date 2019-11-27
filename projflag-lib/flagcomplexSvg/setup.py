import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="flagcomplexSvg",
    version="1.0.1",
    author="Lukas D. Sauer",
    # author_email="author@example.com",
    description="This is a project for visualizing tuples of flags in the real projective plane."
                " It was developed at Heidelberg Institute for Theoretical Studies in 2019.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/LukasDSauer/projflag-lib/projflag-lib/flagcomplexSvg",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        # "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "flagcomplex", # ">=2.00"
        "drawSvg==1.0.2",
    ]
)
