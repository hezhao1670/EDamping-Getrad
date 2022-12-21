# getrad-tracking
e-beam parameters calculation in a ring

# Features
    1. radiation damping: synchrotron radiation
    2. IBS: B-M IBS model with vertical dispersion
    3. BBS: beam-beam scattering effect (from hadron beam) 


# Publications
    Phys. Rev. Accel. Beams 24, 043501
    NAPAC2019 paper https://accelconf.web.cern.ch/napac2019/doi/JACoW-NAPAC2019-TUPLM24.html

# How to use
## directions to compile getrad7 in Linux
    1) move to the directory where the source code is located
    2) run the command
    : ifort -o getrad7 getrad7.f90 getIBS.f90 heating.f90
    3) getrad7 is the output run file


## directions for using getrad7 in Linux
    1) create a subdirectory
    :mkdir s1
    2) move to this directory
    :cd s1
    3) copy input and optics file 
    :cp ../para.in .
    :cp ../optics.dat .
    4) edit the input file
    :gedit para.in &
    5) run code
    :../getrad7 

## In windows, just compile the fortran code and run


# Input file explanantion

## The input file looks like:
    150. 273. 170. 1   3.06E-03 betax0,betay0,slen, ibstype, slipfactor  
    293.2  3e11  2.3e-8 1.9e-8 0.1333 10.e-4 10e-6 4000   gamma,pnum,epsx,epsy,sigs,dponp,dtsim,nsim
    100.  100. 0.  0.  betaxI0,betayI0,dispex,dispix
    9.6e-9  1.5e-9  6.8e-4  0.06 6.9e10 		ion: emittx,emitty,dp/p,sigs,nion
    -1.374    2.03788e+3    8.21626e+3  2.14e+01  -4.9365e+01    1.9912e+2  2.01817e+2   I1-I3,I4x,I4y,I5x,I5y  in mad8 notation
    1. 0. bratio,intmeth
    1. 1. 1.  ibson,heaton,dampon

details:

    betax0,betay0 = beta star at the cooling section of e-ring (unit=m)
    slen = length of cooling section (unit=m)
    ibstype = IBS type (useless in this version)
    slipfactor = slippage factor of the electron ring  
    gamma = central gamma for electrons
    pnum = density of electron beam
    epsx,epsy = initial e-beam rms emittance (unit=m)
    sigs = initial e-beam rms length (unit=m)
    dponp = initial e-beam rms momentum spread 
    dtsim = time step in simulation (unit=s)
    nsim = total number of steps in simulation
    betaxI0,betayI0 = beta function of hadron ring at the cooling section (unit=m)
    dispex = electron dispersion at the cooling section (unit=m)
    dispix = hadron dispersion at the cooling section
    emittx,emitty = rms emittance of hadron beam 
    dp/p,sigs = rms dp/p and bunch length (unit=m) of hadron beam 
    nion = density of hadron beam 
    I1,I2,I3,I4x,I4y,I5x,I5y = synchrotron radiation integrals from MAD
    bratio = scaling for synchrotron radiation (0=no radiation, 1=normal, 2=double the radiation)
    intmeth = calculation mode (1= fast method, 0=step by step tracking)
    ibson = turn on/off IBS (1/0)
    heaton = turn on/off BBS (1/0)
    dampon =turn on/off radiation damping (1/0)


## The input lattice file “optics.dat” must have the format as below:
    The first two parameters (for example: 3011 2) define the length of the lattice data and the first row of the data.
    Data format:  S(m)    BETX(m)    ALFX    BETY(m)    ALFY    DX(m)    DPX    DY(m)    DPY 




# Output file explanantion
## The output file emitt.dat: give the evolution of the damping process
1-nturns  = number of turns
2-epsx  = rms e-beam horizontal emittance (unit=m)
3-epsy = rms e-beam vertical emittance (unit=m)
4-sigs = rms e-beam bunch length (unit=m)
5-dponp = rms momentum spread


## The output files ibs.dat: give the evolution of IBS rates, BBS rates and damping rates
1- nturns  = number of turns
2-4-sumx, sumy, sums = IBS rates (emittance) (unit=1/s)
5-7-dampx, dampy, damps = Radiation damping rates (amplitude) (unit=1/s)
8-10-bbsx, bbsy, bbss   =  BBS rates (emittance) (unit=1/s) 



# Example: calculation in a ring-cooler

This simulation is to calculate the equilibrium e-beam parameters in the ring-cooler. Based on the lattice design and the hadron beam parameters, the final e-beam parameters can be calculated under the effects of IBS, BBS and radiation damping.

## Input file:
    150. 273. 170. 1   3.06E-03 betax0,betay0,slen, ibstype, slipfactor  
    293.2 3e11  2.3e-8 1.9e-8 0.1333 10.e-4 10e-6 40   gamma,pnum,epsx,epsy,sigs,dponp,dtsim,nsim
    100.  100. 0.  0.  betaxI0,betayI0,dispex,dispix
    9.6e-9  1.5e-9  6.8e-4  0.06 6.9e10 	ion: emittx,emitty,dp/p,sigs,nion,intesteps
    -1.374    2.03788e+3    8.21626e+3  2.14e+01  -4.9365e+01    1.9912e+2  2.01817e+2   I1-I3,I4x,I4y,I5x,I5y  in mad8 notation
    1. 1. bratio,intmeth
    1. 1. 1.  ibson,heaton,dampon
    




