# mIcal simulation

This code is Geant4 simulation for mIcal with 10 layers of RPCs. Total 10 RPCs are placed in the center of each layer.

The code is tested for `root6.20.04`, `clhep2404` and `geant4.10.04.p03`.

*Warning:*
- The code has memory leakage. It may overflow memory. I do not recommend it for running more than 1M events.
- One might need to comment out `G4DataQuestionaire` portions in `src/micalPhysicsList.cc` for later geant4 versions.

The field file is `B_mical_hist.root`, the file name is hard-coded in `src/micalFieldPropagator.cc`. The field map could be scaled down. Please look for `fieldxin->Scale(` in the same file.

InputFlag : 
```
0 ->
1 ->
2 ->
3 ->
4 -> Corsika 3D histogram   (tested and working)
5 -> 
```

InputOutput:
```
0: GEN  -> RECO
1: GEN  -> DIGI    (tested and working)
2: GEN  -> SIM     (tested and working)
3: SIM  -> RECO
4: SIM  -> DIGI    (tested and working)
5: DIGI -> RECO 
```


The name of the geometry file `geo_mical_world.gdml` is also hard-coded in `src/vect_manager.cc`.

Please create the file if running for the first time. To create the `gdml` file uncomment the following portion in `mICAL.cc`.
```
  // //  Write the GDML files
  // parser.Write("detector_world.gdml", G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume());
  // cout <<"11xxx "<<endl;
  // system("rm geo_mical_world.gdml");
  // cout <<"12xxx "<<endl;
  // system("cp detector_world.gdml geo_mical_world.gdml");
  // cout <<"13xxx "<<endl;
  // system("rm detector_world.gdml");
  // cout <<"14xxx "<<endl;
```
Then execute it just once. You can break the execution moment you see the output `14xxx`. You can comment the above portion. Recreate the `gdlm` file if you made any changes in the geometry. 


Source in sim01: `source env.sh`

Requirement only once:
```
rm src/Hitsdict.cc
rm src/HitPosdict.cc
rm src/HitPosdict_rdict.pcm
rm src/Hitsdict_rdict.pcm
cd include
rootcint ../src/HitPosdict.cc -c HitPos.h
rootcint ../src/Hitsdict.cc -c Hits.h
cd ..
```

Compile: Create or empty director named `build` and do `cd build && cmake3 .. && make -j 4 && cd ..`

Run: `./mICAL run.mac`


Batch submission to htcondor:
```
1. Create or empty a directory named *runFiles*.
2. Run *./batch_create_mac 2500*. It creates 2500 run files in *runFiles*.
3. Add the lisf of run files in *batch_surya_job.jdl*.
4. Run *condor_submit batch_surya_job.jdl* to submit jobs.
```

*Warning:*
- Change the parameters (i.e. `initialdir`) before submitting jobs.
- The collated file, once read and uncompressed, occupy a lot of memory (almost 2GB). Please adjust the number of jobs per node as per the total memory available.
