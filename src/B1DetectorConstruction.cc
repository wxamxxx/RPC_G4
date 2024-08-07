#include "B1DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

B1DetectorConstruction::~B1DetectorConstruction()
{ }

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // Making Bakelite, CO2, Ar and Gas Mix.

  G4double z, a, density, fractionmass, co2_percent, thickness;
  G4int ncomponents, natoms;
  G4String name, symbol;

  G4Element* elC = new G4Element(name = "Carbon",symbol="C",z = 6., a=   12.0107*g/mole);
  G4Element* elH = new G4Element(name = "Hydrogen",symbol="H",z = 1., a= 1*g/mole);
  G4Element* elO = new G4Element(name = "Oxygen",symbol="O",z = 8., a= 16*g/mole);

  density=1.4*g/cm3;
  G4Material* Bakelite =  new G4Material(name = "Bakelite", density, ncomponents=3);

  density=1.977*kg/m3;

  G4Material* CO2 =  new G4Material(name = "CarbonDioxide", density, ncomponents=2);

  CO2->AddElement(elC,natoms=1);
  CO2->AddElement(elO,natoms=2);

  Bakelite->AddElement(elC,natoms=1);
  Bakelite->AddElement(elH,natoms=4);
  Bakelite->AddElement(elO,natoms=2);

  G4Material* Ar = nist->FindOrBuildMaterial("G4_Ar");


  co2_percent = 30;
  // percent of CO2 in the mixture.
  density = ((3994.8 + 4.062*co2_percent)/2240)*kg/m3; // calculating the density of the mixture of gases at stp.
  G4Material* Gas = new G4Material(name="ArCO2Mix",density,ncomponents=2);
  Gas->AddMaterial(CO2, fractionmass=co2_percent*perCent);
  Gas->AddMaterial(Ar, fractionmass=(100-co2_percent)*perCent);

  // Define Glass
  G4Element* Si = nist->FindOrBuildElement("Si");
  G4Element* O = nist->FindOrBuildElement("O");
  G4Material* Glass = new G4Material("Glass", 2.5*g/cm3, 2);
  Glass->AddElement(Si, 1);
  Glass->AddElement(O, 2);

  // Define R134a
  G4Element* C = nist->FindOrBuildElement("C");
  G4Element* H = nist->FindOrBuildElement("H");
  G4Element* F = nist->FindOrBuildElement("F");
  G4Material* R134a = new G4Material("R134a", 4.55*kg/m3, 3);
  R134a->AddElement(C, 2);
  R134a->AddElement(H, 2);
  R134a->AddElement(F, 4);
  R134a->GetIonisation()->SetMeanExcitationEnergy(12.6*eV);

  // Define SF6
  G4Element* S = nist->FindOrBuildElement("S");
  G4Material* SF6 = new G4Material("SF6", 6.51*kg/m3, 2);
  SF6->AddElement(S, 1);
  SF6->AddElement(F, 6);
  SF6->GetIonisation()->SetMeanExcitationEnergy(16.48*eV);

  // Define Isobutane (C4H10)
  G4Material* Isobutane = new G4Material("Isobutane", 2.59*kg/m3, 2);
  Isobutane->AddElement(C, 4);
  Isobutane->AddElement(H, 10);
  Isobutane->GetIonisation()->SetMeanExcitationEnergy(10.68*eV);

  // Define the gas mixture
  G4Material* RPCGas = new G4Material("RPCGas", 4.46*kg/m3, 3);
  RPCGas->AddMaterial(R134a, 0.947);
  RPCGas->AddMaterial(SF6, 0.003);
  RPCGas->AddMaterial(Isobutane, 0.05);

  G4cout << "Mean Ionization Energy of R134a: " << R134a->GetIonisation()->GetMeanExcitationEnergy() / eV << " eV" << G4endl;
  G4cout << "Mean Ionization Energy of SF6: " << SF6->GetIonisation()->GetMeanExcitationEnergy() / eV << " eV" << G4endl;
  G4cout << "Mean Ionization Energy of C4H10: " << Isobutane->GetIonisation()->GetMeanExcitationEnergy() / eV << " eV" << G4endl;
  G4cout << "Mean Ionization Energy of Gas: " << RPCGas->GetIonisation()->GetMeanExcitationEnergy() / eV << " eV" << G4endl;

  // //Making the geometry now!

  // //making world volume
  // G4double world_x = 0.2*m;
  // G4double world_y = 0.2*m;
  // G4double world_z = 0.2*m;
  // G4Box* worldBox =  new G4Box("World", world_x/2,world_y/2,world_z/2);

  // G4Material* Air = nist->FindOrBuildMaterial("G4_AIR");
  // //making logical volume
  // G4LogicalVolume* worldLog = new G4LogicalVolume(worldBox, Air,"World");
  // //placing it using physical volume
  // G4double pos_x = 0.0*m;
  // G4double pos_y = 0.0*m;
  // G4double pos_z = 0.0*m;
  // G4VPhysicalVolume* World = new G4PVPlacement(0,G4ThreeVector(pos_x, pos_y,pos_z),worldLog,"RPC World",0,false,0);

  thickness = 1;
  // // in mm
  // // Making,logical and placing of three volumes...
  // G4Box* top_bake =  new G4Box("Top_glass", 7.5*cm,7.5*cm,1*mm);
  // G4LogicalVolume* toplog = new G4LogicalVolume(top_bake, Glass,"Top_glass");
  // G4VPhysicalVolume* topbake = new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,(thickness/2 + 1)*mm),toplog,"Top Electrode",worldLog,false,0);

  // G4Box* cavity =  new G4Box("cavity",7.5*cm,7.5*cm,(thickness/2)*mm);
  // G4LogicalVolume* cavitylog = new G4LogicalVolume(cavity, RPCGas,"cavity");
  // G4VPhysicalVolume* Cavity = new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,0*mm),cavitylog,"Gap",worldLog,false,0);

  // G4Box* bot_bake =  new G4Box("Bot_glass", 7.5*cm,7.5*cm,1*mm);
  // G4LogicalVolume* botlog = new G4LogicalVolume(bot_bake , Glass,"Bot_glass");
  // G4VPhysicalVolume* botbake = new G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,-(thickness/2 + 1)*mm),botlog,"Bot Electrode",worldLog,false,0);

  // fScoringVolume = cavitylog;

  // return World;

  //
  // World
  //
  G4double world_sizeXY = 0.2*m;
  G4double world_sizeZ  = 0.2*m;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  auto solidWorld = new G4Box("World",                           // its name
    0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);  // its size

  auto logicWorld = new G4LogicalVolume(solidWorld,  // its solid
    world_mat,                                       // its material
    "World");                                        // its name

  auto physWorld = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                           // at (0,0,0)
    logicWorld,                                // its logical volume
    "World",                                   // its name
    nullptr,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    checkOverlaps);                            // overlaps checking

  //
  // Gas Gap
  //
  G4Material* gas_mat = R134a;
  G4ThreeVector pos2 = G4ThreeVector(0, 0*cm, 0*cm);

  // Gap shape
  G4double gas_dx = 15*cm;
  G4double gas_dy = 15*cm;
  G4double gas_dz  = thickness*mm;
  auto solidGas = new G4Box("Gas",  // its name
    0.5 * gas_dx, 0.5 * gas_dy,
    0.5 * gas_dz);  // its size

  auto logicGas = new G4LogicalVolume(solidGas,  // its solid
    gas_mat,                                        // its material
    "Gas");                                         // its name

  new G4PVPlacement(nullptr,  // no rotation
    pos2,                     // at position
    logicGas,              // its logical volume
    "Gas",                 // its name
    logicWorld,                 // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    checkOverlaps);           // overlaps checking

  // Set Gas as scoring volume
  //
  fScoringVolume = logicGas;

  //
  //always return the physical World
  //
  return physWorld;
}
