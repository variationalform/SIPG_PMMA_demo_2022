# WORK IN PROGRESS

# NOT USE-ABLE YET

## `SIPG_PMMA_demo_2022`

SIPG FEM simulation and comparison of linear elastic and linear viscoelastic response for PMMA used (circa 2022) as described below. This code and the `docker` image (see below) can be used to reproduce the data in the paper:

>*A Priori Analysis of a Symmetric Interior Penalty Discontinuous Galerkin Finite Element Method for a Dynamic Linear Viscoelasticity Model*
>
>**By Yongseok Jang and Simon Shaw**
>>Submitted to: *Computational Methods in Applied Mathematics*

The DOI and full reference will appear __*here*__ if it is accepted.


## Git usage

Some useful Git tips for using these codes, and a trail of how the Git and docker resources were created.

```bash
# first set up the bare (README.md, LICENSE and .gitignore) repo on Git Hub
# clone it and then add files
cd ~/to/parent/directory
git clone git@github.com:variationalform/SIPG_PMMA_demo_2022.git
```
This will get you the basic codes needed to run the simulation. If you want to get up and running more quickly then the `docker` version described below will suit your needs better.

The basic commands are something like this:

```bash
cd ./runtime/le && ./longrun_le.sh | tee ./longrun_le_out.txt
cd ./runtime/ve && ./longrun_ve.sh | tee ./longrun_ve_out.txt
```
These will change to a working directory for linear elasticity or viscoelasticity (`le` or `ve`) and then execute the `bash` scripts. The output is `tee`'d off to `stdout` as well as the `*.txt` text files. This could take many days to run - edit the `bash` scripts to just get the coarser results if you want output more quickly.

 
Git: general usage (some memory joggers)

```bash
git status                       # shows what has changed
git diff                         # shows differences before staging
git add ./README.md .gitignore   # stage
git commit -m 'useful message'   # commit
git diff --staged                # differences after staging
git push                         # push to remote
```


## docker - runtime

The following sequence of commands can be used to generate a suite of results 
for linear elasticity and viscoelasticity using the discrete time SIPG FEM as
described in the paper referenced above. These results increase mesh and time-step
precision with each set taking progressively longer to run. 

```bash
# log in to server, cd to where the shared directories will be, for example
cd fenics_docker/shared
docker pull variationalform/fem:SIPG_PMMA_demo_2022
# you might need this to run screen
sudo chmod o+rw /dev/pts/0
screen -S pmma # at least three: one for le, one for ve and one to monitor progress

# in screen 2
newgrp docker # you may not need this
cd fenics_docker/shared
docker run -ti --name le_pmma -v "$PWD":/home/fenics/shared \
        -w /home/fenics/local/runtime/le variationalform/fem:SIPG_PMMA_demo_2022
# in the container
./longrun_le.sh | tee ./longrun_le_out.txt
# docker rm (e.g.) 3e17 can be used for a fresh start with the container
# periodically, to get interim results, use CNTRL-Z...
cd .. && tar cvf le.tar le && mv le.tar /home/fenics/shared/ && cd le
# and then use fg to resume 

# in screen 3
newgrp docker # you may not need this
cd fenics_docker/shared
docker run -ti --name ve_pmma -v "$PWD":/home/fenics/shared \
        -w /home/fenics/local/runtime/ve variationalform/fem:SIPG_PMMA_demo_2022
# in the container
./longrun_ve.sh | tee ./longrun_ve_out.txt
# docker rm (e.g.) 3e17 can be used for a fresh start with the container
# periodically, to get interim results, use CNTRL-Z...
cd .. && tar cvf ve.tar ve && mv ve.tar /home/fenics/shared/ && cd ve
# and then use fg to resume 

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
will get you back in it. If you `CNTRL-D` out of a container and then try to 
`docker run` again you'll get an error because trhe container name is already in use. 
For example,


```bash
docker run -ti --name le_pmma -v "$PWD":/home/fenics/runtime_le \
    -w /home/fenics/local/runtime/le variationalform/fem:SIPG_PMMA_demo_2022
docker: Error response from daemon: Conflict. The container name "/le_pmma" is
already in use by container
"7aef802a84aca5b063b9e10be9ebc17bc5165a62cc0a111f8e47ccd37381dc36".
You have to remove (or rename) that container to be able to reuse that name.
See 'docker run --help'.
```
In this case, note the ID and just remove the container for a fresh start:

```bash
docker rm 7a
```


## docker - set up

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

# to set up the local files ready for running the solves...
cd ~/local
# we copy the directory structure from a remote through the share
cp -r ~/shared/runtime/ .
# it's from a mac so get rid of these things...
rm ./runtime/.DS_Store ./runtime/ve/.DS_Store ./runtime/le/.DS_Store 
ls -R

# alter RUNPATH in the sh files if necessary and then these are the run commands
fenics@329e2fdd6283:~/local/runtime/le$ ./longrun_le.sh | tee ./longrun_le_out.txt
fenics@329e2fdd6283:~/local/runtime/ve$ ./longrun_ve.sh | tee ./longrun_ve_out.txt

```
This has given us a working version of FEniCS in a docker image. Its runtime can be realized with a docker container. The instructions for doing that are above.

For interest, 

```bash
# to run with a share but working in a local directory
docker run -ti --name SIPG_PMMA_demo_2022 -v "$PWD":/home/fenics/shared \
        -w /home/fenics/local quay.io/fenicsproject/stable:2017.2.0.r4
# in FEniCS
cp -r ../shared/runtime/ .
# you might need this to run screen
sudo chmod o+rw /dev/pts/0
screen -S <name>
# can re-attach with screen -d -r <name> (-d if not cleanly detached)
# run bash in several screens sessions, and then in each of two sessions
cd ~/local/runtime/le ; ./longrun_le.sh | tee ./longrun_le_out.txt
cd ~/local/runtime/ve ; ./longrun_ve.sh | tee ./longrun_ve_out.txt
```

## docker hub

The tag was created on <https://hub.docker.com/repository/docker/variationalform/fem> with the following:

```bash
# get the tag
docker ps --filter "status=exited"
# commit it
docker commit 329
# gives
sha256:afcfbfc0646d40a89f300da5d6e4e453be5099a9875940a54835444497b59fda
# login, create and push
docker login
docker tag afcfb variationalform/fem:SIPG_PMMA_demo_2022
docker push variationalform/fem:SIPG_PMMA_demo_2022
```

Now a downloader can execute and run as explained at the top of this section

Alternatively, once the tar file, made like this,

```bash
docker save afcfb > SIPG_PMMA_demo_2022_docker_afcfb.tar
```

is available, the container can be instyanced with this:

```bash
docker load < SIPG_PMMA_demo_2022_docker_afcfb.tar
docker run -ti afcfb
```

**_*TO DO:*_** update <https://hub.docker.com/repository/docker/variationalform/fem> with the paper's DOI

**_*TO DO:*_** make the tarfile available (figshare?)

