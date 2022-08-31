# `SIPG_PMMA_demo_2022`

Simulation and comparison of linear elastic and linear viscoelastic response for PMMA using SIPG FEM used (circa 2022) as described below. This code and the `docker` image (see below) can be used to reproduce the data in the paper

>*A Priori Analysis of a Symmetric Interior Penalty Discontinuous Galerkin Finite Element Method for a Dynamic Linear Viscoelasticity Model*
>
>**By Yongseok Jang and Simon Shaw**
>>Submitted to: *Computational Methods in Applied Mathematics*

The DOI and full reference will appear if it is accepted.


# Git usage

Some useful Git tips for using these codes.

```bash
# first set up the bare RERADME, licence and .gitignore) repo on Git Hub
# clone it and then add files
git clone git@github.com:variationalform/SIPG_PMMA_demo_2022.git
```

General usage

```bash
git status
git diff                         # shows differences before staging
git add ./README.md .gitignore   # stage
git commit -m 'useful message'   # commit
git diff --staged                # differences aftre staging
git push                         # push to remote
```




# docker - runtime

The following sequence of commands can be used to generate a suite of results 
for linear elasticity and viscoelasticity which increasing mesh and time-step
precision.

```bash
# to run it with a share but working in a local directory
docker run -ti --name SIPG_PMMA_demo_2022 -v "$PWD":/home/fenics/shared \
        -w /home/fenics/local quay.io/fenicsproject/stable:2017.2.0.r4
cp -r ../shared/runtime/ .
screen -S <name>
# can re-attach with screen -d -r <name> (-d if not cleanly detached)
# run bash in several sessions, and then in each of two sessions
cd ~/local/runtime/ve ; ./longrun_le.sh | tee ./longrun_le_out.txt
cd ~/local/runtime/ve ; ./longrun_ve.sh | tee ./longrun_ve_out.txt
```

These runs could take up to ten days (or even longer) depending on your hardware.
The `*.sh` scripts can edited to get shorter (but less accurate) runs.

You can detach from the screen session while the runs execute. For `docker` these are useful. If `CNTRL-D` was used to exit the container then re-attach like this:

```bash
# list containers - note the ID on the left (e.g. adf...)
docker ps -a
# to start a stopped container
docker start adf
docker attach adf
```
On the other hand if `CNTRL-p CNTRL-q` was used to leave the container then

```bash
docker attach adf
```


# docker - set up

Fir interest, this is how the image was set up. We used the following sequence of commands. The FENICS image was selected from those available at <https://quay.io/repository/fenicsproject/stable?tab=tags>


```bash
docker pull quay.io/fenicsproject/stable:2017.2.0.r4
# to run it when working in the shared directory (not a good idea in general)
docker run -ti --name SIPG_PMMA_demo_2022 -v "$PWD":/home/fenics/shared \
        -w /home/fenics/shared quay.io/fenicsproject/stable:2017.2.0.r4
# we'll need 'screen'
sudo apt-get update
sudo apt-get install screen
screen
# this failed with the following message:
# Cannot open your terminal '/dev/pts/0' - please check.
# solve it with this:
sudo chmod o+rw /dev/pts/0

# to et the FENICS version:
python
>>> import dolfin as df; print (df.__version__); exit()

```

