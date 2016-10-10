/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Naveen Michaud-Agrawal (Johns Hopkins U)
     K-space terms added by Stan Moore (BYU)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <string.h>
#include "compute_group_group_atom.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "kspace.h"
#include "error.h"
#include <math.h>
#include "comm.h"
#include "domain.h"
#include "math_const.h"
#include "memory.h"

#include <iostream>
using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.00001

ComputeGroupGroupAtom::ComputeGroupGroupAtom(LAMMPS *lmp, int narg, char **arg) :
    ComputeGroupGroup(lmp, narg, arg),
    carray(NULL),
    nmax(0)
{
    if (narg < 4) error->all(FLERR,"Illegal compute group/group command");

    peratom_flag      = 1;  // Indicating a peratom compute
    size_peratom_cols = 4;  // # of Columns per atom.
    extarray          = 0;  // 0/1 if global array is all intensive/extensive
    scalar_flag       = 0;
    vector_flag       = 0;
}


ComputeGroupGroupAtom::~ComputeGroupGroupAtom()
{
    memory->destroy(carray);
}


void ComputeGroupGroupAtom::compute_peratom()
{
    // grow array if necessary
    if (atom->nmax > nmax) {

        memory->destroy(carray);
        nmax = atom->nmax;
        memory->create(carray, nmax, size_peratom_cols, "group/group/atom:carray");
        array_atom = carray;

    }

    if (pairflag) pair_contribution();
    if (kspaceflag) kspace_contribution(); // This doesn't happen though. See compute_group_group.cpp constructor.
}


void ComputeGroupGroupAtom::pair_contribution()
{
    int i,j,ii,jj,inum,jnum,itype,jtype;
    double xtmp,ytmp,ztmp,delx,dely,delz;
    double rsq,eng,fpair,factor_coul,factor_lj;
    int *ilist,*jlist,*numneigh,**firstneigh;

    double **x = atom->x;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    double *special_coul = force->special_coul;
    double *special_lj = force->special_lj;
    int newton_pair = force->newton_pair;
    double *columns;

    // invoke half neighbor list (will copy or build if necessary)

    neighbor->build_one(list);

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    // loop over neighbors of my atoms
    // skip if I,J are not in 2 groups


    for (ii = 0; ii < inum; ii++) {
        i = ilist[ii];

        // skip if atom I is not in either group
        if (!(mask[i] & groupbit || mask[i] & jgroupbit)) continue;

        xtmp = x[i][0];
        ytmp = x[i][1];
        ztmp = x[i][2];
        itype = type[i];
        jlist = firstneigh[i];
        jnum = numneigh[i];

        for (jj = 0; jj < jnum; jj++) {
            j = jlist[jj];
            factor_lj = special_lj[sbmask(j)];
            factor_coul = special_coul[sbmask(j)];
            j &= NEIGHMASK;

            // skip if atom J is not in either group
            if (!(mask[j] & groupbit || mask[j] & jgroupbit)) continue;

            int ij_flag = 0;
            int ji_flag = 0;
            if (mask[i] & groupbit && mask[j] & jgroupbit) ij_flag = 1;
            if (mask[j] & groupbit && mask[i] & jgroupbit) ji_flag = 1;

            // skip if atoms I,J are only in the same group
            if (!ij_flag && !ji_flag) continue;

            delx = xtmp - x[j][0];
            dely = ytmp - x[j][1];
            delz = ztmp - x[j][2];
            rsq = delx*delx + dely*dely + delz*delz;
            jtype = type[j];

            if (rsq < cutsq[itype][jtype]) {
                eng = pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

                // energy only computed once so tally full amount
                // force tally is jgroup acting on igroup

                if (newton_pair || j < nlocal) {
                    array_atom[i][0] += eng;
                    if (ij_flag) {
                        array_atom[i][1] += delx*fpair;
                        array_atom[i][2] += dely*fpair;
                        array_atom[i][3] += delz*fpair;
                    }
                    if (ji_flag) {
                        array_atom[j][1] -= delx*fpair;
                        array_atom[j][2] -= dely*fpair;
                        array_atom[j][3] -= delz*fpair;
                    }

                    // energy computed twice so tally half amount
                    // only tally force if I own igroup atom
                }
                else {
                    array_atom[i][0] += 0.5*eng;
                    if (ij_flag) {
                        array_atom[i][1] += delx*fpair;
                        array_atom[i][2] += dely*fpair;
                        array_atom[i][3] += delz*fpair;
                    }
                }
            }
        }
    }
}



