# WORK IN PROGRESS

# NOT USE-ABLE YET

# `SIPG_PMMA_demo_2022`

Simulation and comparison of linear elastic and linear viscoelastic response for PMMA using SIPG FEM used (circa 2022) as described below. This code and the `docker` image (see below) can be used to reproduce the data in the paper

>*A Priori Analysis of a Symmetric Interior Penalty Discontinuous Galerkin Finite Element Method for a Dynamic Linear Viscoelasticity Model*
>
>**By Yongseok Jang and Simon Shaw**
>>Submitted to: *Computational Methods in Applied Mathematics*

The DOI and full reference will appear __*here*__ if it is accepted.


# Git usage

Some useful Git tips for using these codes.

```bash
# first set up the bare (README.md, LICENSE and .gitignore) repo on Git Hub
# clone it and then add files
cd ~/to/parent/directory
git clone git@github.com:variationalform/SIPG_PMMA_demo_2022.git
```
This will get you the basic codes needed to run the simulation. If you want to get up and running more quickly then the `docker` version described below will suit your needs better.
 
Git: general usage (some memory joggers)

```bash
git status                       # shows what has changed
git diff                         # shows differences before staging
git add ./README.md .gitignore   # stage
git commit -m 'useful message'   # commit
git diff --staged                # differences after staging
git push                         # push to remote
```


# docker - runtime

The following sequence of commands can be used to generate a suite of results 
for linear elasticity and viscoelasticity using the discrete time SIPG FEM as
described in the paper referenced above. These results increase mesh and time-step
precision with each set taking progressively longer to run. 

```bash
# to run with a share but working in a local directory
docker run -ti --name SIPG_PMMA_demo_2022 -v "$PWD":/home/fenics/shared \
        -w /home/fenics/local quay.io/fenicsproject/stable:2017.2.0.r4
# in FEniCS
cp -r ../shared/runtime/ .
screen -S <name>
# can re-attach with screen -d -r <name> (-d if not cleanly detached)
# run bash in several screens sessions, and then in each of two sessions
cd ~/local/runtime/ve ; ./longrun_le.sh | tee ./longrun_le_out.txt
cd ~/local/runtime/ve ; ./longrun_ve.sh | tee ./longrun_ve_out.txt
```

These runs could take up to ten days (or even longer) depending on your hardware.
The `*.sh` scripts can be edited to get shorter (but less accurate) runs.

You can detach from the screen session while the runs execute. For `docker` the commands below are useful. First, if `CNTRL-D` was used to exit the container then you can re-attach like this:

```bash
# list containers - note the ID on the left (e.g. adf...)
docker ps -a
# to start and then attach to a stopped container
docker start adf
docker attach adf
```
On the other hand if `CNTRL-p CNTRL-q` was used to leave the container then

```bash
docker attach adf
```
will get you back in it.


# docker - set up

For interest, this is how the image was set up. We used the following sequence of commands. The FEniCS image was selected from those available at <https://quay.io/repository/fenicsproject/stable?tab=tags>. It's quite old but we've been using it for years and know it to be reliable.


```bash
docker pull quay.io/fenicsproject/stable:2017.2.0.r4
# to run it when working in the shared directory (not a good idea in general)
docker run -ti --name SIPG_PMMA_demo_2022 -v "$PWD":/home/fenics/shared \
        -w /home/fenics/shared quay.io/fenicsproject/stable:2017.2.0.r4
# we'll need 'screen'
sudo apt-get update
sudo apt-get install screen
screen -S <name>
# this failed with the following message:
# Cannot open your terminal '/dev/pts/0' - please check.
# solve it with this:
sudo chmod o+rw /dev/pts/0

# screen helpers:
screen -ls && screen -r <name>
#toggle previous two C-A C-A; next terminal C-AN; switch to shell n C-A n
#C-a (shift to u/c) A <term_name>
#C-AD to detach and exit
#screen -d -r (detaches and reattaches when necessary)
#CA ? gives a list of commands


# to determine the FENICS version:
python
>>> import dolfin as df; print (df.__version__); exit()

```
This has given us a working version of FEniCS in a docker image. It's runtime can be realized with a docker container. The instructions for doing that are above.


## MUTABLE - this below needs update with this content

The tag was created on <https://hub.docker.com/repository/docker/variationalform/puretime> with the following:

```bash
# get the tag
docker ps --filter "status-exited"
# commit it
docker commit 19d1e9956cbd
# gives
sha256:c8ff6b020b5ca02cb686ba19fbd9613f4ae8cf7e9204a0f4d657fb284707decf
# Then
docker login
docker tag c8ff6b020 variationalform/puretime:fouvol
docker push variationalform/puretime:fouvol
```

Now a downloader can execute and run as explained at the top of this section

Alternatively, once the tar file, made like this,

```bash
docker save c8ff6b020 > fouvol_docker_c8ff.tar
docker load < fouvol_docker_c8ff.tar
```

is available,

```bash
docker run -ti c8ff
```

**_*TO DO:*_** update <https://hub.docker.com/repository/docker/variationalform/puretime> with the paper's DOI

**_*TO DO:*_** make the tarfile available (figshare?)

