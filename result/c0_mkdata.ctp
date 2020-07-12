 30.5    131.70  0      :ALAT0,ALNG0 ICORD
 -                      :FLOBSPD  (input; obs. position and displacement) [FL3]
./fwd/nankai_jcg.dat    :FLOBSPOS (input; obs. point for fwd-analysis) [FL31]
 1.0                    :DELT
 1                      :NSOURCE(number of faults)
----------------------------------------------
 1                      :LS
 -                      :FLF(input; fault surface position)
./NE_gsi3.cof           :FCOF(LS) OR "PLANE FAULT" (output; B-splined fault surface) [FL1]
  33.00        0.0      :BPHI(LS),BDLT(LS)(in degree)
 30.50  131.70          :BLAT(LS),BLON(LS)
1620.0 337.5            :BLEN(LS),BWID(LS)(in km)
  20   11   3              : JU0(LS),JV0(LS),NDGF(LS)
 - - - - - - - - - - - - - - - - - - - - - - -
/3.11 TOHOKU EQ        /:HD_JAC(LS) (jacobi header comments) [CJA]
../tmp/grdthk4.invjac   :INVJAC(LS)(output; inv-Jacobian, binary) [FL2]
../tmp/grdthk4.fwdjac   :FWDJAC(LS)(output; fwd-Jacobian, binary) [FL21]
 31.00 132.50           :BLAF(LS),BLOF(LS)
 750.0 150.0  0.0       :BLEF(LS),BWIF(LS),BDEP(LS)
   8  3  60  0           :JS0(LS),JT0(LS),JD2(LS),NDGS(LS)
========================================================
/CH BD TOHOKU          /:[HD] Non-Negative Analysis(result header comments)
-------------------------------------------------------
  1                     :LSIN
./fwd/mkdata.scd        :SCD(LS) (output; source displacement)
  1  1  0  0            :JUST(1,LS),JUST(2,LS),JUST(3,LS),JUST(4,LS)
  0  0  0  0     no     :NNLS(1,LS),NNLS(2,LS),NNLS(3,LS),NNLS(4,LS)
 -90.0 0.0 0.0          :RSS(LS),ROP(LS),REX(LS)
  1  8  1  3  36        :IUS(1,LS),IUS(2,LS),IUS(3,LS),IUS(4,LS),MLG(LS)
  0  0  0  0            :IBOC(1,LS),IBOC(2,LS),IBOC(3,LS),IBOC(4,LS)
--------------------------------------------------------
 -                      :CAL (output; calculated movement at the obs.)
./fwd/mkdata.fmp        :FMP (output; B-splined scd-file)
./fwd/mkdata.fwd        :FWD (output; calculated movement at given points in fl31)
  1.0 1.0 1.0           :FACT(1),FACT(2),FACT(3)
  1.0 1.0 1.0           :FABS(1),FABS(2),FABS(3)
  0.0                   :VRP
--------------------------------------------------------
  30.  38. 127.  141.   :GLAT(MIN,MAX), GLON(MIN,MAX)
../share/gps_nankai.dat :FGPS (input)
./fwd/gps_nankai.dat    :FPOS (output)
../share/grid_dummy.dat :FDUM (input; dummy)
   3.   15.    2.       :E_GSI, E_JCG (ERROR GSI,JCG [mm]), Height_error
  100.                  :ALPHA (to search, set 100.)
