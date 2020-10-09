#include "micalParameterMessenger.hh"
#include <libconfig.h++>

using namespace libconfig;
using namespace std;

micalParameterMessenger *micalParameterMessenger::AnPointer;

micalParameterMessenger::micalParameterMessenger()
{
  AnPointer = this;
  //ctor
  Config cfg;
  // cout<<"===================================================="<<endl;
  // Read the config file. If there is an error,
  // report it and exit.
  try {
    // cfg.readFile("/products/GEANT4.10/ICALDOCS/detector_config.cfg");
    cfg.readFile("detector_config.cfg");
    // cout <<"In file detector_config.cfg =>"<<endl;
  } catch
    (const FileIOException& fioex) {
    G4cerr << "I/O error while reading file." << G4endl;
    G4cerr << "Setting the parameters to default values ...."
           << G4endl;
    SetParameterLocation("local");
    SetGeometryLocation("local");
    SetStripInfoLocation("local");
    SetnFeThickness(56.0);
    SetnAirGap(40.0);
    SetnLayer(150);
  }

  //Setup reading of the detector parameters
  // Get the parameter location.
  try {
    string location = cfg.lookup("parameter_location");
    cout <<"parameter location = "<< location<<endl;
    SetParameterLocation(location);
  } catch
    (const SettingNotFoundException& nfex) {
    SetParameterLocation("local");
    G4cerr << "No 'parameter_location' setting in configuration file."
           << G4endl;
    G4cerr << "Setting the parameters to default values ...."
           << G4endl;
  }

  // Get the geometry location
  try {
    string location = cfg.lookup("geometry_location");
    //    cout <<"geometry location = "<< location<<endl;
    SetGeometryLocation(location);
    cout<<"Setting GeomoetryLocation = "<<location<<" ... "<<endl;
  } catch
    (const SettingNotFoundException& nfex) {
    SetGeometryLocation("local");
    cout<<"Setting GeometryLocation = local ... "<<endl;
    G4cerr << "No 'geometry_location' setting in configuration file."
           << G4endl;
    G4cerr << "Setting the geometry to default value ...."
           << G4endl;
  }

  // Get the StripInfo location
  try {
    string location = cfg.lookup("stripInfo_location");
    //    cout <<"stripInfo location = "<< location<<endl;
    SetStripInfoLocation(location);
    cout<<"Setting StripInfoLocation = "<<location<<" ... "<<endl;
  } catch
    (const SettingNotFoundException& nfex) {
    cout<<"Setting StripInfoLocation = local ... "<<endl;
    SetStripInfoLocation("local");
    G4cerr << "No 'SetStripInfoLocation' setting in configuration file."
           << G4endl;
    G4cerr << "Setting the geometry to default value ...."
           << G4endl;
  }

  double ironThickness;
  bool cfgiron = cfg.lookupValue("iron_thickness",ironThickness);
  if(cfgiron) {
    cout<<"cfgiron "<<cfgiron<<endl;
    SetnFeThickness(ironThickness);
    cout<<"Setting IronThickness = "<<ironThickness<<" mm ... "<<endl;
  } else {
    SetnFeThickness(56.0);
    cout<<"cfgiron "<<cfgiron<<endl;
    cout << "No 'iron_thickness' setting in configuration file."
           << endl;
    cout << "Setting the iron thickness to default value 56 mm ...."
           << endl;
  }

  double airThickness;
  bool cfgair = cfg.lookupValue("air_thickness",airThickness);
  if(cfgair) {
    cout<<"cfgair "<<cfgair<<endl;
    SetnAirGap(airThickness);
    cout<<"Setting AirThickness = "<<airThickness<<" mm ... "<<endl;
  } else {
    SetnAirGap(40.0);
    cout<<"cfgair "<<cfgair<<endl;
    cout << "No 'air_thickness' setting in configuration file." << endl;
    cout << "Setting the air thickness to default value 40 mm ...." << endl;
  }

  int layern;
  bool cfglayern = cfg.lookupValue("n_layer",layern);
  if(cfglayern) {
    cout<<"cfglayern "<<cfglayern<<endl;
    SetnLayer(layern);
    cout<<"Setting number of layers = "<<layern<<" ... "<<endl;
  } else {
    SetnLayer(150);
    cout << "No 'n_layer' setting in configuration file."
           << endl;
    cout << "Setting the number of layers to default value 150 ...." << endl;

  }

  double strwdXY;
  bool cfgstrwdXY = cfg.lookupValue("strwd_xy",strwdXY);
  if(cfgstrwdXY) {
    cout<<"cfglayern "<<cfglayern<<endl;
    SetXYstrwd(strwdXY);
    cout<<"Setting X-Y stripwidth = "<<strwdXY<<" mm ... "<<endl;
  } else {
    SetXYstrwd(30.0);
    cout << "No 'strwd_xy' setting in configuration file." << endl;
    cout << "Setting the stripwidth to default value 30 mm ...." << endl;
  }

//   // Get the file version
//   try {
//     string location = cfg.lookup("file_version");
//     //    cout <<"file_version = "<< location<<endl;
//     SetFileVersion(location);
//   } catch
//     (const SettingNotFoundException& nfex) {
//     SetFileVersion("v1.00");
//     G4cerr << "No 'file_version' setting in configuration file."
//            << G4endl;
//     G4cerr << "Setting the file version to v1.00 ...."
//            << G4endl;
//   }
}

micalParameterMessenger::~micalParameterMessenger()
{
  //dtor
}
