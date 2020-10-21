#ifndef WORKSPACE_H
#define WORKSPACE_H

/*
 * This file was autogenerated by OSQP-Python on October 21, 2020 at 13:09:50.
 * 
 * This file contains the prototypes for all the workspace variables needed
 * by OSQP. The actual data is contained inside workspace.c.
 */

#include "types.h"
#include "qdldl_interface.h"

// Data structure prototypes
extern csc Pdata;
extern csc Adata;
extern c_float qdata[45];
extern c_float ldata[48];
extern c_float udata[48];
extern OSQPData data;

// Settings structure prototype
extern OSQPSettings settings;

// Scaling structure prototypes
extern c_float Dscaling[45];
extern c_float Dinvscaling[45];
extern c_float Escaling[48];
extern c_float Einvscaling[48];
extern OSQPScaling scaling;

// Prototypes for linsys_solver structure
extern csc linsys_solver_L;
extern c_float linsys_solver_Dinv[93];
extern c_int linsys_solver_P[93];
extern c_float linsys_solver_bp[93];
extern c_float linsys_solver_sol[93];
extern c_float linsys_solver_rho_inv_vec[48];
extern c_int linsys_solver_Pdiag_idx[45];
extern csc linsys_solver_KKT;
extern c_int linsys_solver_PtoKKT[45];
extern c_int linsys_solver_AtoKKT[120];
extern c_int linsys_solver_rhotoKKT[48];
extern QDLDL_float linsys_solver_D[93];
extern QDLDL_int linsys_solver_etree[93];
extern QDLDL_int linsys_solver_Lnz[93];
extern QDLDL_int   linsys_solver_iwork[279];
extern QDLDL_bool  linsys_solver_bwork[93];
extern QDLDL_float linsys_solver_fwork[93];
extern qdldl_solver linsys_solver;

// Prototypes for solution
extern c_float xsolution[45];
extern c_float ysolution[48];

extern OSQPSolution solution;

// Prototype for info structure
extern OSQPInfo info;

// Prototypes for the workspace
extern c_float work_rho_vec[48];
extern c_float work_rho_inv_vec[48];
extern c_int work_constr_type[48];
extern c_float work_x[45];
extern c_float work_y[48];
extern c_float work_z[48];
extern c_float work_xz_tilde[93];
extern c_float work_x_prev[45];
extern c_float work_z_prev[48];
extern c_float work_Ax[48];
extern c_float work_Px[45];
extern c_float work_Aty[45];
extern c_float work_delta_y[48];
extern c_float work_Atdelta_y[45];
extern c_float work_delta_x[45];
extern c_float work_Pdelta_x[45];
extern c_float work_Adelta_x[48];
extern c_float work_D_temp[45];
extern c_float work_D_temp_A[45];
extern c_float work_E_temp[48];

extern OSQPWorkspace workspace;

#endif // ifndef WORKSPACE_H
