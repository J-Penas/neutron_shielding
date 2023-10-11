//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file Run.hh
/// \brief Definition of the Run class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>

#include "G4Event.hh"
#include "G4THitsMap.hh"

class DetectorConstruction;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*);
   ~Run();

  public:
    void CountProcesses(const G4VProcess* process);                  
    void ParticleCount(G4String, G4double);
    void SumTrackLength (G4int,G4int,G4double,G4double,G4double,G4double);
    
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);    
    void EndOfRun(); 

    virtual void RecordEvent(const G4Event*);        
    virtual void Merge(const G4Run*);
   
  private:
    struct ParticleData {
     ParticleData()
       : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
     ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
       : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
     G4int     fCount;
     G4double  fEmean;
     G4double  fEmin;
     G4double  fEmax;
    };

    G4int nEvent;
    G4int re0ID, re1ID, re2ID, re3ID, re4ID, re5ID, reTID;
    G4THitsMap<G4double> re0;
    G4THitsMap<G4double> re1;
    G4THitsMap<G4double> re2;
    G4THitsMap<G4double> re3;
    G4THitsMap<G4double> re4;
    G4THitsMap<G4double> re5;
    G4THitsMap<G4double> reT;
    
     
  private:
    DetectorConstruction* fDetector;
    G4ParticleDefinition* fParticle;
    G4double              fEkin;
        
    std::map<G4String,G4int>        fProcCounter;            
    std::map<G4String,ParticleData> fParticleDataMap;
        
    G4int    fNbStep1, fNbStep2;
    G4double fTrackLen1, fTrackLen2;
    G4double fTime1, fTime2;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

