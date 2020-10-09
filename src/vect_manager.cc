#include "vect_manager.h"

InoStrip_Manager ::InoStrip_Manager() {    APointer = this;  }
InoStrip_Manager:: ~InoStrip_Manager(){}
InoStrip_Manager* InoStrip_Manager::APointer;

InoStripX_Manager ::InoStripX_Manager() {    APointer = this;  }
InoStripX_Manager:: ~InoStripX_Manager(){}
InoStripX_Manager* InoStripX_Manager::APointer;

InoStripY_Manager ::InoStripY_Manager() {    APointer = this;  }
InoStripY_Manager:: ~InoStripY_Manager(){}
InoStripY_Manager* InoStripY_Manager::APointer;

InoHit_Manager ::InoHit_Manager() {    APointer = this;  }
InoHit_Manager:: ~InoHit_Manager(){}
InoHit_Manager *InoHit_Manager::APointer;

InoCluster_Manager ::InoCluster_Manager() {    APointer = this;  }
InoCluster_Manager:: ~InoCluster_Manager(){}
InoCluster_Manager *InoCluster_Manager::APointer;

InoTrack_Manager ::InoTrack_Manager() {    APointer = this;  }
InoTrack_Manager:: ~InoTrack_Manager(){}
InoTrack_Manager *InoTrack_Manager::APointer;

InoFittedTrack_Manager ::InoFittedTrack_Manager() {    APointer = this;  }
InoFittedTrack_Manager:: ~InoFittedTrack_Manager(){}
InoFittedTrack_Manager *InoFittedTrack_Manager::APointer;

InoTrackCand_Manager ::InoTrackCand_Manager() {    APointer = this;  }
InoTrackCand_Manager:: ~InoTrackCand_Manager(){}
InoTrackCand_Manager *InoTrackCand_Manager::APointer;

InoGeometry_Manager ::InoGeometry_Manager() 
{    
  cout <<"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx "<<endl;
  icalGeometry        = TGeoManager::Import("geo_mical_world.gdml");//mical_world_alleql.gdml");//geo.gdml");
  // icalGeometry        = TGeoManager::Import("/products/GEANT4.10/ICALDOCS/geo.gdml");
  //icalGeometry        = TGeoManager::Import("detector_world.gdml");

  APointer = this; }
InoGeometry_Manager:: ~InoGeometry_Manager(){}
InoGeometry_Manager *InoGeometry_Manager::APointer;
/*
InoMuRange_Manager :: InoMuRange_Manager(){APointer = this;}
InoMuRange_Manager ::~InoMuRange_Manager(){}
InoMuRange_Manager *InoMuRange_Manager::APointer;
*/


InoRPCStrip_Manager :: InoRPCStrip_Manager(){APointer = this;}
InoRPCStrip_Manager ::~InoRPCStrip_Manager(){}
InoRPCStrip_Manager *InoRPCStrip_Manager::APointer;
