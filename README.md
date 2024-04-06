# MD simulation instructions

## Software requirements:

- `GROMACS (2019.6-sp-gcc8.4-mvapich2.3-cuda10.2)`  

## Gain protein and ligand structure

We must download the protein structure file we will be working with.  We will utilize MTNN  (PDB ID: 3DF9)， NFSB(PDB ID:  3X22)， YDHF(PDB ID:  1OG6) and  THIF(PDB ID:  1ZFN). 

ligand structure：https://atb.uq.edu.au

## Workflow

There are 7 conventional sequential steps in the workflow:

1. Generate Topology

2. Define box and Solvate

3. Add Ions

4. Energy Minimization

5. Equilibration

6. Production MD

7. Analysis

   

   PS: Take MTNN as an example

## Step 1: Generate Topology

There should now be a "gromos54a7.ff" subdirectory in your working directory. Write the topology for the mtnn with pdb2gmx:

```sh
#protein gro and itp file
gmx pdb2gmx -f pykf.pdb -o pykf.gro -water spc -ignh -ter
```

ligand gro file anf itp file

```sh
#ligand gro file
gmx editconf -f DF9.pdb -o DF9.gro
```



### Build the Topology

Including the parameters for the ``` DF9```ligand in the system topology . Just insert a line that says `#include "DF9.itp"` into topol.top after the position restraint file is included. The inclusion of position restraints indicates the end of the "Protein" moleculetype section.

```sh
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include water topology
#include "./gromos54a7.ff/spc.itp"
```

becomes...

```sh
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include ligand topology
#include "DF9.itp"

; Include water topology
#include "./gromos54a7_atb.ff/spc.itp"
```

The last adjustment to be made is in the `[ molecules ]` directive. To account for the fact that there is a new molecule in complex.gro, we have to add it here, like so:

```sh
[ molecules ]
; Compound        #mols
Protein_chain_B     1
BWUB                1
```

The topology and coordinate file are now in agreement with respect to the contents of the system.

## Step 2: Define box and Solvate

At this point, the workflow is just like any other MD simulation. We will define the unit cell and fill it with water.

```sh
gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro
```

## Step 3: Add Ions

Use grompp to assemble a .tpr file, using any .mdp file.  ``` ions.mdp ```

```sh
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
```

We now pass our .tpr file to genion:

```shell
gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -np 10
```

The names of the ions specified with -pname and -nname were force field-specific in previous versions of GROMACS, but were standardized in version 4.5. The specified atom names are always the elemental symbol in all capital letters, along with the `[ moleculetype ]`. Residue names may or may not append the sign of the charge (+/-). Refer to ions.itp for proper nomenclature if you encounter difficulties.

Your `[ molecules ]` directive should now look like:

```sh
[ molecules ]
; Compound        #mols
Protein_chain_B     1
BWUB                1
SOL         123865
NA               10
```

## Step 4: Energy Minimization

Now that the system is assembled, create the binary input using grompp using ``` em.mdp``` input parameter file:

```sh
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr 
```

We are now ready to invoke mdrun to carry out the EM:

```sh
gmx mdrun -v -deffnm em
```

## Step 5: Equilibration

### Restraining the Ligand

To restrain the ligand, we will need to generate a position restraint topology for it. First, create an index group for JZ4 that contains only its non-hydrogen atoms:

```shell
gmx make_ndx -f em.gro -o index_DF9.ndx
...
 > 0 & ! a H*
 > q
```

Then, execute the genrestr module and select this newly created index group (which will be group 3 in the index_jz4.ndx file):

```shell
 gmx genrestr -f em.gro -n index_DF9.ndx  -o posre_DF9.itp -fc 1000 1000 1000
```

Now, we need to include this information in our topology. We can do this in several ways, depending upon the conditions we wish to use. If we simply want to restrain the ligand whenever the protein is also restrained, add the following lines to your topology in the location indicated:

```shell
; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include ligand topology
#include "DF9.itp"

; Ligand position restraints
#ifdef POSRES_LIG
#include "posre_DF9.itp"
#endif

; Include water topology
#include "gromos54a7_atb.ff/spc.itp"
```

### Thermostats

Proper control of temperature coupling is a sensitive issue. Coupling every moleculetype to its own thermostatting group is a bad idea. For instance, if you do the following:

```
tc-grps = Protein Non-Protein
```

Here, I am using em.gro, the output (minimized) structure of our system:

```
gmx make_ndx -f em.gro -o index.ndx
  0 System              : 373784 atoms
  1 Protein             :  2141 atoms
  2 Protein-H           :  1709 atoms
  3 C-alpha             :   232 atoms
  4 Backbone            :   696 atoms
  5 MainChain           :   929 atoms
  6 MainChain+Cb        :  1139 atoms
  7 MainChain+H         :  1157 atoms
  8 SideChain           :   984 atoms
  9 SideChain-H         :   780 atoms
 10 Prot-Masses         :  2141 atoms
 11 non-Protein         : 371643 atoms
 12 Other               : 371643 atoms
 13 BWUB                :    38 atoms
 14 SOL                 : 371595 atoms
 15 NA                  :    10 atoms

```

Merge the "Protein" and "fbp" groups with the following, where ">" indicates the make_ndx prompt:

```
> 1 | 13 （ligand)
> q
```

Proceed with *NVT* equilibration using ``` nvt.mdp```  file.

```sh
gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr -maxwarn 2
gmx_mpi mdrun -v -deffnm nvt
```

Once the *NVT* simulation is complete, proceed to *NPT* with ``` npt.mdp``` file

```sh
gmx_mpi grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr -maxwarn 2
gmx_mpi mdrun -v -deffnm npt
```

## Step 6: Production MD

Upon completion of the two equilibration phases, the system is now well-equilibrated at the desired temperature and pressure. We are now ready to release the position restraints and run production MD for data collection. The process is just like we have seen before, as we will make use of the checkpoint file (which in this case now contains preserve pressure coupling information) to grompp. We will run a 100-ns MD simulation, the script for which can be found ``` md.mdp```

```sh
gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top  -o md_0_10.tpr -maxwarn 2
gmx_mpi mdrun -v -deffnm md_0_10
```

## Step 7: Analysis

**Periodic limit**

```shell
 gmx_mpi trjconv -s md_0_10.tpr -f md_0_10.xtc -o md_0_10_center.xtc -center -pbc nojump -ur compact
 
 #Choose "Protein" for centering and "System" for output.将之前的-pbc mol改为-pbc nojump
 
gmx_mpi trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o md_0_10_fit.xtc -fit rot+trans
 
 #Choose "Backbone" to perform least-squares fitting to the protein backbone, and "System" for output. 
```

**RMSD**

```shell
gmx_mpi rms -s md_0_10.tpr -f md_0_10_center.xtc -o rmsd_pro.xvg -tu ns
```

**Extract structures**

```shell
gmx_mpi trjconv -s md_0_10.tpr -f md_0_10_center.xtc -o 100ns.pdb -dump 100000
 
 #Select group for output  "System"
```

**Distance**

```shell
#3df9_BWUB
gmx_mpi make_ndx -f em.gro -o distance.ndx

  0 System              : 263783 atoms
  1 Protein             :  2141 atoms
  2 Protein-H           :  1709 atoms
  3 C-alpha             :   232 atoms
  4 Backbone            :   696 atoms
  5 MainChain           :   929 atoms
  6 MainChain+Cb        :  1139 atoms
  7 MainChain+H         :  1157 atoms
  8 SideChain           :   984 atoms
  9 SideChain-H         :   780 atoms
 10 Prot-Masses         :  2141 atoms
 11 non-Protein         : 261642 atoms
 12 Other               : 261642 atoms
 13 BWUB                :    38 atoms
 14 SOL                 : 261594 atoms
 15 NA                  :    10 atoms

 nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
 'r': residue              'res' nr         'chain' char
 "name": group             'case': case sensitive           'q': save and quit
 'ri': residue index

> 13 & a O5

Copied index group 13 'BWUB'
Found 0 atoms with name O5

 16 BWUB                :    38 atoms

> 13 & a O1

Copied index group 13 'BWUB'
Found 1 atoms with name O1
Merged two groups with AND: 38 1 -> 0
Group is empty

> 13 & a O2

Copied index group 13 'BWUB'
Found 1 atoms with name O2
Merged two groups with AND: 38 1 -> 0
Group is empty

> 13 & a O3

Copied index group 13 'BWUB'
Found 1 atoms with name O3
Merged two groups with AND: 38 1 -> 1

 17 BWUB_&_O3           :     1 atoms

> 1 & r 47 & a HZ2     

Copied index group 1 'Protein'
Merged two groups with AND: 2141 13 -> 13
Found 11 atoms with name HZ2
Merged two groups with AND: 13 11 -> 1

 18 Protein_&_r_47_&_HZ2:     1 atoms

> 17 | 18  

Copied index group 17 'BWUB_&_O3'
Copied index group 18 'Protein_&_r_47_&_HZ2'
Merged two groups with OR: 1 1 -> 2

 19 BWUB_&_O3_Protein_&_r_47_&_HZ2:     2 atoms

> q

gmx_mpi distance -s md_0_10.tpr -f md_0_10_center.xtc -select 19 -oall -oav -n distance.ndx -tu ns -len 2 -binw 0.01

xmgrace distave.xvg 

#3df9K47_BWUB
gmx_mpi make_ndx -f em.gro -o distance.ndx

  0 System              : 263789 atoms
  1 Protein             :  2146 atoms
  2 Protein-H           :  1715 atoms
  3 C-alpha             :   232 atoms
  4 Backbone            :   696 atoms
  5 MainChain           :   929 atoms
  6 MainChain+Cb        :  1139 atoms
  7 MainChain+H         :  1157 atoms
  8 SideChain           :   989 atoms
  9 SideChain-H         :   786 atoms
 10 Prot-Masses         :  2146 atoms
 11 non-Protein         : 261643 atoms
 12 Other               : 261643 atoms
 13 BWUB                :    38 atoms
 14 SOL                 : 261594 atoms
 15 NA                  :    11 atoms

 nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
 'r': residue              'res' nr         'chain' char
 "name": group             'case': case sensitive           'q': save and quit
 'ri': residue index

> 13 & a O3

Copied index group 13 'BWUB'
Found 1 atoms with name O3
Merged two groups with AND: 38 1 -> 1

 16 BWUB_&_O3           :     1 atoms

> 1 & r 47 & a HI1

Copied index group 1 'Protein'
Merged two groups with AND: 2146 18 -> 18
Found 1 atoms with name HI1
Merged two groups with AND: 18 1 -> 1

 17 Protein_&_r_47_&_HI1:     1 atoms

> 16 | 17

Copied index group 16 'BWUB_&_O3'
Copied index group 17 'Protein_&_r_47_&_HI1'
Merged two groups with OR: 1 1 -> 2

 18 BWUB_&_O3_Protein_&_r_47_&_HI1:     2 atoms

> q

gmx_mpi distance -s md_0_10.tpr -f md_0_10_center.xtc -select 18 -oall -oav -n distance.ndx -tu ns -len 2 -binw 0.01

xmgrace distave.xvg 

```



**Hydrogen bonds**

```shell
gmx_mpi hbond -s md_0_10.tpr -f md_0_10_center.xtc -b 50 -e 100 -num hbond_num.xvg -tu ns
```



**Binding free energies**

```shell
gmx_mpi trjconv -f md_0_10_center.xtc -o md_50_100_trj.xtc -b 50000 -e 100000 -dt 100
gmx_mpi make_ndx -f em.gro -o mmpbsa.ndx
  0 System              : 263783 atoms
  1 Protein             :  2141 atoms
  2 Protein-H           :  1709 atoms
  3 C-alpha             :   232 atoms
  4 Backbone            :   696 atoms
  5 MainChain           :   929 atoms
  6 MainChain+Cb        :  1139 atoms
  7 MainChain+H         :  1157 atoms
  8 SideChain           :   984 atoms
  9 SideChain-H         :   780 atoms
 10 Prot-Masses         :  2141 atoms
 11 non-Protein         : 261642 atoms
 12 Other               : 261642 atoms
 13 BWUB                :    38 atoms
 14 SOL                 : 261594 atoms
 15 NA                  :    10 atoms

 nr : group      '!': not  'name' nr name   'splitch' nr    Enter: list groups
 'a': atom       '&': and  'del' nr         'splitres' nr   'l': list residues
 't': atom type  '|': or   'keep' nr        'splitat' nr    'h': help
 'r': residue              'res' nr         'chain' char
 "name": group             'case': case sensitive           'q': save and quit
 'ri': residue index

> 1 | 13

Copied index group 1 'Protein'
Copied index group 13 'BWUB'
Merged two groups with OR: 2141 38 -> 2179

 16 Protein_BWUB        :  2179 atoms

> q
bash gmx_mmpbsa.bsh  >>gmx_mmpbsa.log
```

