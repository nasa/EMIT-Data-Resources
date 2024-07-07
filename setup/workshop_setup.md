# Cloud Workspace Setup

> If you plan to use this repository with the **Openscapes 2i2c JupyterHub Cloud Workspace** there are no additional setup requirements for the Python environment. All packages needed are included unless specified within a notebook, in which case a cell will be dedicated to installing the necessary Python libraries using the appropriate package manager.

After completing the [prerequisites](prerequisites.md) you will have access to the Openscapes 2i2c JupyterHub cloud workspace. [Click here to start JupyterLab](https://workshop.openscapes.2i2c.cloud/). Use your email and the provided password to sign in. This password will be provided in the workshop. If you're interested in using the 2i2c cloud workspace outside of the workshop, please [contact us](#contact-info).

After signing in you will be prompted for some server options:

<p align="left">
  <img src="../img/2i2c_server_selection.png" width=100%/>
</p>

Be sure to select the radio button for **Python** and a size of **14.8 GB RAM and up to 3.75 CPUs**.

At this point you can use the terminal to clone the repository.

## Cloning the EMIT-Data-Resources Repository

If you plan to edit or contribute to the EMIT-Data-Resources repository, we recommend following a fork and pull workflow: first fork the repository, then clone your fork to your local machine, make changes, push changes to your fork, then make a pull request back to the main repository. An example can be found in the [CONTRIBUTING.md](../CONTRIBUTING.md) file.

If you just want to work with the notebooks or modules, you can simply clone the repository.

To clone the repository, navigate to the directory where you want to store the repository on your local machine, then type the following:

```cmd
git clone https://github.com/nasa/EMIT-Data-Resources.git
```

## Troubleshooting

We recommend *Shutting down all kernels* after running each notebook. This will clear the memory used by the previous notebook, and is necessary to run some of the more memory intensive notebooks.

<p align="left">
  <img src="../img/shut_down_kernels.png" width=100%/>
</p>

No single notebook exceeds roughly the limit using the provided data, but **if you choose to use your own data in the notebook, or have 2 notebooks open and do not shut down the kernel, you may get an out of memory error.**  

If you elect to try this on your own data/ROI, you may need to select a larger server size. This will often happen if you are using the last EMIT scene from an orbit. In some cases those can be almost double the size of a normal scene. Please select the smallest possible.

---

## Contact Info  

Email: <LPDAAC@usgs.gov>  
Voice: +1-866-573-3222  
Organization: Land Processes Distributed Active Archive Center (LP DAAC)¹  
Website: <https://lpdaac.usgs.gov/>  
Date last modified: 05-24-2024  

¹Work performed under USGS contract 140G0121D0001 for NASA contract NNG14HH33I.
