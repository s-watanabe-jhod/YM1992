 31.0    132.6   0      :ALAT0,ALNG0 ICORD
 -                      0 FLOBSPD  (input; obs. position and displacement) [FL3]
 -                      0 FLOBSPOS (input; obs. point for fwd-analysis) [FL31]
 -                      0 DELT
 1                      : NSOURCE(number of faults)
----------------------------------------------
 1                      : LS
./flt/block_rot_LM9.flt : FLT(input; fault surface position)
./cof/NE_gsi3+9.cof     : FCOF(LS) OR "PLANE FAULT" (output; B-splined fault surface) [FL1]
 303.00        0.0      :BPHI(LS),BDLT(LS)(in degree)
 36.21  137.92           :BLAT(LS),BLON(LS)
 150.0 100.0            :BLEN(LS),BWID(LS)(in km)
   5  3  3              : JU0(LS),JV0(LS),NDGF(LS)
 - - - - - - - - - - - - - - - - - - - - - - -
 -                      0 HD_JAC(LS) (jacobi header comments) [CJA]
 -                      0 INVJAC(LS)(output; inv-Jacobian, binary) [FL2]
 -                      0 FWDJAC(LS)(output; fwd-Jacobian, binary) [FL21]
  0  0                  0 BLAF(LS),BLOF(LS)
  0  0  0               0 BLEF(LS),BWIF(LS),BDEP(LS)
  0  0  0  0            0 JS0(LS),JT0(LS),JD2(LS),NDGS(LS)
========================================================
  -                     0 [HD] Non-Negative Analysis(result header comments)
-------------------------------------------------------
  0                     0 LSIN
  -                     0 SCD(LS) (output; source displacement)
  0  0  0  0            0 JUST(1,LS),JUST(2,LS),JUST(3,LS),JUST(4,LS)
  0  0  0  0            0 NNLS(1,LS),NNLS(2,LS),NNLS(3,LS),NNLS(4,LS)
  0  0  0  0            0 RSS(LS),ROP(LS),REX(LS)
  0  0  0  0  0         0 IUS(1,LS),IUS(2,LS),IUS(3,LS),IUS(4,LS),MLG(LS)
  0  0  0  0  0         0 IBOC(1,LS),IBOC(2,LS),IBOC(3,LS),IBOC(4,LS)
--------------------------------------------------------
  -                     0 CAL (output; calculated movement at the obs.)
  -                     0 FMP (output; B-splined scd-file)
  -                     0 FWD (output; calculated movement at given points in fl31)
  0.0 0.0 0.0           0 FACT(1),FACT(2),FACT(3)
  0.0 0.0 0.0           0 FABS(1),FABS(2),FABS(3)
  0.0                   0 VRP
--------------------------------------------------------
 0 0 0 0                0 GLAT(MIN,MAX), GLON(MIN,MAX)1
 -                      0 
 -                      0 
 -                      0 
  0.  0.                0 E_GSI, E_JCG (ERROR GSI,JCG [mm])
 100.                   0 ALPHA (to search, set 100.)
