#include "matrixsolver.h"

//*********************Topology-considered coarsening*****************************/

void matrixsolver::solveamgaggl_ovst_grid_auto(matrix &matlevel0, matrix &matlevel1, amggrid &gridlevel0, amggrid &gridlevel1, int nlev) {

    unsigned long nrows;
    unsigned long ncols;
    unsigned long nci;
    unsigned long ncit;
    unsigned long ncnb;
    unsigned long ncnbt;
    unsigned int nbnum;
    unsigned int nbnfull;

    unsigned long nbci;
    unsigned long nbcit;
    unsigned long cellcur;
    unsigned long cellnb;
    unsigned long cscur;
    unsigned long csnb;

    unsigned long nbmax;

    double epsilon0 = 1e-12;
    double epsilon1 = 1e-16;

    amggrid *csgrid;
    unsigned long csngroup;

    nrows = matlevel0.nrows;
    ncols = matlevel0.ncols;

    matrix ncflag(nrows);

    csgrid = new amggrid;
    csgrid->cellsnum = 0;
    csgrid->cellsmax = nrows;
    csgrid->ngroup = 4;
    csgrid->amgcells = new amgcell[nrows];
    csgrid->xnci = new unsigned long[nrows];
    csgrid->init();

    csngroup = 0;

    /*---------agglomentration and additve correction step---------*/

    nci = 0;

    for (unsigned long i = 0; i < nrows; i++) {

        cellcur = i;

        if (gridlevel0.amgcells[cellcur].ostatus == 0) {

            csngroup++;
            ncflag(cellcur) = 1;

            continue;

        }

        if (gridlevel0.amgcells[cellcur].ostatus == 2) {

            continue;
        }

        if ((int) ncflag(cellcur) > 0) {

            continue;

        }

        nbnum = 0;
        nbnfull = 0;


        for (unsigned int j = 1; j < ncols / 2; j++) {

            ncnb = matlevel0(cellcur, 2 * j);

            if ((fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon0)&& ((int) ncflag(ncnb) == 1)) {

                nbci = csgrid->xnci[ncnb];

                if (csgrid->amgcells[nbci].cellsnum < csgrid->amgcells[nbci].ngroup) {

                    nbnfull = 1;
                    break;

                }

            }

        }

        if (nbnfull == 1) {

            csgrid->amgcells[nbci].cells[csgrid->amgcells[nbci].cellsnum] = cellcur;
            csgrid->amgcells[nbci].cellsnum++;

            ncflag(cellcur) = 1;
            csgrid->xnci[cellcur] = nbci;
            csgrid->amgcells[nbci].blkid = gridlevel0.amgcells[cellcur].blkid;

            continue;

        }


        csgrid->amgcells[nci].cellsnum = 1;
        csgrid->amgcells[nci].cells[0] = cellcur;
        csgrid->amgcells[nci].ostatus = 1;
        ncflag(cellcur) = 1;
        csgrid->xnci[cellcur] = nci;
        csgrid->cellsnum++;


        for (unsigned int j = 1; j < ncols / 2; j++) {

            ncnb = matlevel0(cellcur, 2 * j);

            if ((fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon0) && ((int) ncflag(ncnb) != 1) && ((gridlevel0.amgcells[cellcur].blkid == gridlevel0.amgcells[ncnb].blkid))) {

                if (csgrid->amgcells[nci].cellsnum < csgrid->amgcells[nci].ngroup) {

                    csgrid->amgcells[nci].cells[csgrid->amgcells[nci].cellsnum] = ncnb;
                    csgrid->amgcells[nci].cellsnum++;

                    ncflag(ncnb) = 1;
                    csgrid->xnci[ncnb] = nci;
                    csgrid->amgcells[nci].blkid = gridlevel0.amgcells[ncnb].blkid;
                    nbnum++;

                    if ((gridlevel0.amgcells[cellcur].ostatus == 2)|| (gridlevel0.amgcells[ncnb].ostatus == 2)) {

                        csgrid->amgcells[nci].ostatus = 2;

                    }

                }

            }


        }

        if (nbnum > 0) {

            nci++;

            continue;

        }

        if (nbnum == 0) {

            for (unsigned int j = 1; j < ncols / 2; j++) {

                ncnbt = matlevel0(cellcur, 2 * j);

                if (fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon0) {

                    nbci = csgrid->xnci[ncnbt]; // initialize the coarse cell selection
                    nbnum++;

                    break;

                }
            }

            if (nbnum > 0) {

                for (unsigned int j = 1; j < ncols / 2; j++) {

                    ncnbt = matlevel0(cellcur, 2 * j);

                    if (fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon0) {

                        nbcit = csgrid->xnci[ncnbt];

                        if (csgrid->amgcells[nbcit].cellsnum < csgrid->amgcells[nbci].cellsnum) {

                            nbci = nbcit;

                        }

                    }
                }

                if (csgrid->amgcells[nbci].cellsnum == csgrid->amgcells[nbci].ngroup) {

                    csgrid->amgcells[nbci].reshape(csgrid->amgcells[nbci].ngroup + 1);

                }

                csgrid->amgcells[nci].cellsnum = 0;
                csgrid->amgcells[nci].cells[0] = 0;
                csgrid->cellsnum--;

                csgrid->amgcells[nbci].cells[csgrid->amgcells[nbci].cellsnum] = cellcur;
                csgrid->amgcells[nbci].cellsnum++;

                ncflag(cellcur) = 1;
                csgrid->xnci[cellcur] = nbci;
                csgrid->amgcells[nbci].blkid = gridlevel0.amgcells[cellcur].blkid;

            }
        }

    }

    for (unsigned long i = 0; i < nrows; i++) {

        cellcur = i;

        for (unsigned int j = 1; j < ncols / 2; j++) {

            ncnb = matlevel0(cellcur, 2 * j);

            if (fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon1) {

                if (ncflag(ncnb) == 1) {

                    continue;

                }

                if (((gridlevel0.amgcells[cellcur].blkid == gridlevel0.amgcells[ncnb].blkid)) && (gridlevel0.amgcells[ncnb].ostatus == 2)) {

                    ncit = csgrid->xnci[cellcur];

                    if (csgrid->amgcells[ncit].cellsnum == csgrid->amgcells[ncit].ngroup) {

                        csgrid->amgcells[ncit].reshape(csgrid->amgcells[ncit].ngroup + 1);

                    }

                    csgrid->amgcells[ncit].cells[csgrid->amgcells[ncit].cellsnum] = ncnb;
                    csgrid->amgcells[ncit].cellsnum++;

                    ncflag(ncnb) = 1;
                    csgrid->xnci[ncnb] = ncit;
                    csgrid->amgcells[ncit].blkid = gridlevel0.amgcells[ncnb].blkid;

                    csgrid->amgcells[ncit].ostatus = 2;

		    }

            }

        }

    }


    if (csngroup != 0) {

        csgrid->amgcells[nci].ngroup = csngroup;
        csgrid->amgcells[nci].cellsnum = 0;
        csgrid->amgcells[nci].ostatus = 0;
        csgrid->amgcells[nci].nbsnum = 0;
        csgrid->amgcells[nci].cells = new unsigned long[csngroup];

        for (unsigned long i = 0; i < nrows; i++) {

            cellcur = i;

            if (gridlevel0.amgcells[cellcur].ostatus == 0) {

                csgrid->amgcells[nci].cells[csgrid->amgcells[nci].cellsnum] = cellcur;
                csgrid->xnci[cellcur] = nci;
                csgrid->amgcells[nci].cellsnum++;
                csgrid->amgcells[nci].blkid = gridlevel0.amgcells[cellcur].blkid;

            }

        }

        nci++;
        csgrid->cellsnum++;

    }

    gridlevel1.cellsmax = csgrid->cellsnum;
    gridlevel1.cellsnum = csgrid->cellsnum;
    gridlevel1.ngroup = 4;
    gridlevel1.amgcells = new amgcell[gridlevel1.cellsnum];
    gridlevel1.xnci = new unsigned long[gridlevel0.cellsnum];

    gridlevel1.gridcopy(csgrid);

    for (unsigned long i = 0; i < nrows; i++) {

        cellcur = i;
        gridlevel1.xnci[cellcur] = csgrid->xnci[cellcur];

    }

    //matlevel1.initmatrix(gridlevel1.cellsnum, 2 * (20 + 1)); // nbmax is wrong because it's not the number of neighbours.

    matrix nbs(gridlevel1.cellsnum);
    matrix csmat(gridlevel1.cellsnum, 60);

    for (unsigned long i = 0; i < matlevel0.nrows; i++) {

        cellcur = i;

        cscur = gridlevel1.xnci[cellcur]; // xnci should be part of gridlevel0

        csmat(cscur, 0) = cscur;
        csmat(cscur, 1) = csmat(cscur, 1) + matlevel0(cellcur, 1);

        for (unsigned long j = 1; j < matlevel0.ncols / 2; j++) {

            if (fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon1) {

                cellnb = matlevel0(cellcur, 2 * j);

                csnb = gridlevel1.xnci[cellnb];


                if (csnb == cscur) {

                    csmat(cscur, 1) = csmat(cscur, 1) + matlevel0(cellcur, 2 * j + 1);

                }

                if (csnb != cscur) {

                    for (unsigned k = 1; k < csmat.ncols / 2; k++) {

                        if (((unsigned long) csmat(cscur, 2 * k) == 0) && (fabs(csmat(cscur, 2 * k + 1)) < epsilon1)) {

                            csmat(cscur, 2 * k) = csnb;

                            csmat(cscur, 2 * k + 1) = matlevel0(cellcur, 2 * j + 1);

                            nbs(cscur) = nbs(cscur) + 1;

                            break;

                        }

                        if (((unsigned long) csmat(cscur, 2 * k) == csnb) && (fabs(csmat(cscur, 2 * k + 1)) > epsilon1)) {

                            csmat(cscur, 2 * k + 1) = csmat(cscur, 2 * k + 1) + matlevel0(cellcur, 2 * j + 1);

                            break;

                        }

                    }

                }

            }
        }

    }

    nbmax = 0;

    for (unsigned long i = 0; i < nbs.nrows; i++) {

        nbnum = nbs(i);

        if (nbnum > nbmax) {

            nbmax = nbnum;

        }

    }

    matlevel1.initmatrix(gridlevel1.cellsnum, 2 * (nbmax + 1));

    for (unsigned long i = 0; i < matlevel1.nrows; i++) {

        for (unsigned long j = 0; j < matlevel1.ncols; j++) {

            matlevel1(i, j) = csmat(i, j);

        }

    }

    csgrid->cellsclear();

}


//*********************Topology-neglected coarsening*****************************/

void matrixsolver::solveamgaggl_grid_auto(matrix &matlevel0, matrix &matlevel1, amggrid &gridlevel0, amggrid &gridlevel1) {

    unsigned long nrows;
    unsigned long ncols;
    unsigned long nci;
    unsigned long ncnb;
    unsigned long ncnbt;
    unsigned int nbnum;
    unsigned int nbnfull;
    unsigned long nbci;
    unsigned long nbcit;
    unsigned long cellcur;
    unsigned long cellnb;
    unsigned long cscur;
    unsigned long csnb;

    //    unsigned long cellp0;

    unsigned long nbmax;

    double epsilon0 = 1e-12;
    double epsilon1 = 1e-16;

    amggrid *csgrid;


    nrows = matlevel0.nrows;
    ncols = matlevel0.ncols;

    matrix ncflag(nrows, 1);

    csgrid = new amggrid;
    csgrid->cellsnum = 0;
    csgrid->cellsmax = nrows;
    csgrid->ngroup = 4;
    csgrid->amgcells = new amgcell[nrows];
    csgrid->xnci = new unsigned long[nrows];
    csgrid->init();

    /*---------agglomentration and additve correction step---------*/

    nci = 0;

    for (unsigned long i = 0; i < nrows; i++) {

        cellcur = i;

        if ((int) ncflag(cellcur) > 0) {

            continue;

        }

        nbnum = 0;
        nbnfull = 0;

        for (unsigned int j = 1; j < ncols / 2; j++) {

            ncnb = matlevel0(cellcur, 2 * j);

            if ((fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon0)&& ((int) ncflag(ncnb) == 1)) {
                nbci = csgrid->xnci[ncnb];

                if (csgrid->amgcells[nbci].cellsnum < csgrid->amgcells[nbci].ngroup) {

                    nbnfull = 1;
                    break;

                }

            }

        }

        if (nbnfull == 1) {

            csgrid->amgcells[nbci].cells[csgrid->amgcells[nbci].cellsnum] = cellcur;
            csgrid->amgcells[nbci].cellsnum++;

            ncflag(cellcur) = 1;
            csgrid->xnci[cellcur] = nbci;

            continue;
        }

        csgrid->amgcells[nci].cellsnum = 1;
        csgrid->amgcells[nci].cells[0] = cellcur;
        ncflag(cellcur) = 1;
        csgrid->xnci[cellcur] = nci;
        csgrid->cellsnum++;

        for (unsigned int j = 1; j < ncols / 2; j++) {

            ncnb = matlevel0(cellcur, 2 * j);

            if ((fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon0) && ((int) ncflag(ncnb) != 1)) {

                if (csgrid->amgcells[nci].cellsnum < csgrid->amgcells[nci].ngroup) {

                    csgrid->amgcells[nci].cells[csgrid->amgcells[nci].cellsnum] = ncnb;
                    csgrid->amgcells[nci].cellsnum++;

                    ncflag(ncnb) = 1;
                    csgrid->xnci[ncnb] = nci;

                    nbnum++;

                }

            }

        }

        if (nbnum > 0) {

            nci++;

            continue;

        }

        if (nbnum == 0) {

            for (unsigned int j = 1; j < ncols / 2; j++) {

                ncnb = matlevel0(cellcur, 2 * j);

                if (fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon0) {

                    nbci = csgrid->xnci[ncnb]; // initialize the coarse cell selection

                    nbnum++;
                    break;

                }
            }

            if (nbnum > 0) {

                for (unsigned int j = 1; j < ncols / 2; j++) {

                    ncnbt = matlevel0(cellcur, 2 * j);

                    if (fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon0) {

                        nbcit = csgrid->xnci[ncnbt];

                        if (csgrid->amgcells[nbcit].cellsnum < csgrid->amgcells[nbci].cellsnum) {

                            nbci = nbcit;

                        }

                    }
                }

                if (csgrid->amgcells[nbci].cellsnum == csgrid->amgcells[nbci].ngroup) {

                    csgrid->amgcells[nbci].reshape(csgrid->amgcells[nbci].ngroup + 1);

                }

                csgrid->amgcells[nci].cellsnum = 0;
                csgrid->amgcells[nci].cells[0] = 0;
                csgrid->cellsnum--;

                csgrid->amgcells[nbci].cells[csgrid->amgcells[nbci].cellsnum] = cellcur;
                csgrid->amgcells[nbci].cellsnum++;

                ncflag(cellcur) = 1;
                csgrid->xnci[cellcur] = nbci;

            }

        }

    }


    gridlevel1.cellsmax = csgrid->cellsmax;
    gridlevel1.cellsnum = csgrid->cellsnum;
    gridlevel1.ngroup = 4;
    gridlevel1.amgcells = new amgcell[gridlevel1.cellsnum];
    gridlevel1.xnci = new unsigned long[gridlevel0.cellsnum];

    gridlevel1.gridcopy(csgrid);

    for (unsigned long i = 0; i < nrows; i++) {

        cellcur = i;
        gridlevel1.xnci[cellcur] = csgrid->xnci[cellcur];

    }

    matrix nbs(gridlevel1.cellsnum);
    matrix csmat(gridlevel1.cellsnum, 60);

    for (unsigned long i = 0; i < matlevel0.nrows; i++) {

        cellcur = i;

        cscur = gridlevel1.xnci[cellcur]; // xnci should be part of gridlevel0

        csmat(cscur, 0) = cscur;
        csmat(cscur, 1) = csmat(cscur, 1) + matlevel0(cellcur, 1);

        for (unsigned long j = 1; j < matlevel0.ncols / 2; j++) {

            if (fabs(matlevel0(cellcur, 2 * j + 1)) > epsilon1) {

                cellnb = matlevel0(cellcur, 2 * j);

                csnb = gridlevel1.xnci[cellnb];


                if (csnb == cscur) {

                    csmat(cscur, 1) = csmat(cscur, 1) + matlevel0(cellcur, 2 * j + 1);

                }

                if (csnb != cscur) {

                    for (unsigned k = 1; k < csmat.ncols / 2; k++) {


                        if (((unsigned long) csmat(cscur, 2 * k) == 0) && (fabs(csmat(cscur, 2 * k + 1)) < epsilon1)) {

                            csmat(cscur, 2 * k) = csnb;

                            csmat(cscur, 2 * k + 1) = matlevel0(cellcur, 2 * j + 1);

                            nbs(cscur) = nbs(cscur) + 1;

                            break;

                        }


                        if (((unsigned long) csmat(cscur, 2 * k) == csnb) && (fabs(csmat(cscur, 2 * k + 1)) > epsilon1)) {

                            csmat(cscur, 2 * k + 1) = csmat(cscur, 2 * k + 1) + matlevel0(cellcur, 2 * j + 1);

                            break;

                        }

                    }

                }

            }
        }

    }

    nbmax = 0;

    for (unsigned long i = 0; i < nbs.nrows; i++) {

        nbnum = nbs(i);

        if (nbnum > nbmax) {

            nbmax = nbnum;
        }

    }

    matlevel1.initmatrix(gridlevel1.cellsnum, 2 * (nbmax + 1));

    for (unsigned long i = 0; i < matlevel1.nrows; i++) {

        for (unsigned long j = 0; j < matlevel1.ncols; j++) {

            matlevel1(i, j) = csmat(i, j);
        }

    }

    csgrid->cellsclear();

}
