# #CLHEP
# export CLHEP_BASE_DIR=/products/GEANT4/10.04.p03/CLHEP/clhep2404
# export PATH=$CLHEP_BASE_DIR/bin:$PATH
# export LD_LIBRARY_PATH=$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH
# #GEANT
# export G4INSTALL=/products/GEANT4/10.04.p03/geant4.10.04.p03-install/share/Geant4-10.4.3/geant4make
# # export G4WORKDIR=/home/user/G4WORK
# source $G4INSTALL/geant4make.sh
# source /products/GEANT4/10.04.p03/geant4.10.04.p03-install/bin/geant4.sh

# source /products/root/root-6.20.00/bin/thisroot.sh



#ROOT 6.20.04
source /products/genie/v3_0_6_sl7/ROOT6/root-6.20.04/bin/thisroot.sh

#CLHEP
export CLHEP_BASE_DIR=/products/GEANT4/sl7/CLHEP/clhep2404
export PATH=$CLHEP_BASE_DIR/bin:$PATH
export LD_LIBRARY_PATH=$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH
#GEANT
export G4INSTALL=/products/GEANT4/sl7/geant4.10.04.p03-install/share/Geant4-10.4.3/geant4make
# export G4WORKDIR=/home/user/G4WORK
source $G4INSTALL/geant4make.sh
source /products/GEANT4/sl7/geant4.10.04.p03-install/bin/geant4.sh

