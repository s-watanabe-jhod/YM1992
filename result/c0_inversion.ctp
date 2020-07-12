 30.5    131.7   0      :ALAT0,ALNG0 ICORD
./fwd/vel_data_YN.dat   :FLOBSPD  (input; obs. position and displacement) [FL3]
./fwd/vel_data_YN.dat   :FLOBSPOS (input; obs. point for fwd-analysis) [FL31]
 1.0                    :DELT
 1                      :NSOURCE(number of faults)
----------------------------------------------
 1                      :LS
 -                      :FLF(input; fault surface position)
./NE_gsi3.cof           :FCOF(LS) OR "PLANE FAULT" (output; B-splined fault surface) [FL1]
  33.00        0.0      :BPHI(LS),BDLT(LS)(in degree)
 30.50  131.70          :BLAT(LS),BLON(LS)
1620.0 337.5            :BLEN(LS),BWID(LS)(in km)
  20   11  3            : JU0(LS),JV0(LS),NDGF(LS)
 - - - - - - - - - - - - - - - - - - - - - - -
/ SAMPLE DATA          /:HD_JAC(LS) (jacobi header comments) [CJA]
../tmp/grdsmpl.invjac   :INVJAC(LS)(output; inv-Jacobian, binary) [FL2]
../tmp/grdsmpl.fwdjac   :FWDJAC(LS)(output; fwd-Jacobian, binary) [FL21]
 31.00 132.50           :BLAF(LS),BLOF(LS)
 750.0 150.0  0.0       :BLEF(LS),BWIF(LS),BDEP(LS)
  24   8   6   3        :JS0(LS),JT0(LS),JD2(LS),NDGS(LS)
========================================================
/CH BD SAMPLE          /:[HD] Non-Negative Analysis(result header comments)
-------------------------------------------------------
  1                     :LSIN
./data/res.scd          :SCD(LS) (output; source displacement)
  1  1  0  0            :JUST(1,LS),JUST(2,LS),JUST(3,LS),JUST(4,LS)
  0  0  0  0     no     :NNLS(1,LS),NNLS(2,LS),NNLS(3,LS),NNLS(4,LS)
 -90.0 0.0 0.0          :RSS(LS),ROP(LS),REX(LS)
  1 27  1  11  36        :IUS(1,LS),IUS(2,LS),IUS(3,LS),IUS(4,LS),MLG(LS)
  0  0  0  1            :IBOC(1,LS),IBOC(2,LS),IBOC(3,LS),IBOC(4,LS)
--------------------------------------------------------
./data/res.cal          :CAL (output; calculated movement at the obs.)
./data/res.fmp          :FMP (output; B-splined scd-file)
./data/res.fwd          :FWD (output; calculated movement at given points in fl31)
  1.0 1.0 1.0           :FACT(1),FACT(2),FACT(3)
  1.0 1.0 1.0           :FABS(1),FABS(2),FABS(3)
  0.0                   :VRP
--------------------------------------------------------
 0 0 0 0                :GLAT(MIN,MAX), GLON(MIN,MAX)1
 -                      :
 -                      :
 -                      :
  0.  0.                :E_GSI, E_JCG (ERROR GSI,JCG [mm])
 100.                   :ALPHA (to_search,_set_100.)
