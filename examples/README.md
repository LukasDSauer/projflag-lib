# A short explanation for the example tutorials

All of the tutorials are Jupyter notebooks that can be used for interactively trying
the `flagcomplex` code.

## Installation instructions

In order to run the tutorials, we need to install the packages `flagcomplex` and `flagcomplexSvg`.
Even though step 1 is optional, I recommend installing them inside a virtual environment (and not globally).

1. Install a virtual Python3 environment (check the internet, if you don't know how).
2. Install the `flagcomplex` package. A very low-brow installation would be copying the
    package files to your hard drive, and then execute
    ```
   pip3 install -e path/to/flagcomplex
    ```
   This installation leaves everything in place and simply creates a link from your virtual environment
   to the package files.
   
   The more sophisticated installation is
   ```
   pip3 install path/to/flagcomplex
    ```
   which actually copies the files to your virtual environment or
   
   ```
   pip3 install "git+<url-to-git-repo>#egg=flagcomplex&subdirectory=projflag-lib/flagcomplex"
   ```
   if you want to install directly from a git repository (e.g. from GitHub).
   
   If you don't need it anymore
   ```
   pip3 uninstall flagcomplex
   ```
   will remove the package.

3. Install the `flagcomplexSvg` package the same way as the `flagcomplex` package.
   
4. If you want to visualize movies with the tutorial notebooks, you will
   need to install `moviepy` as well.
5. Start the notebook with
   ```
   jupyter-notebook tutorial-<name>.ipynb
   ```
   from **inside your virtual environment** and enjoy!
   
