#include "Particles.h"

TDatabasePDG* gPDG     = new TDatabasePDG();
TParticlePDG* gProton  = gPDG->GetParticle(2212);
TParticlePDG* gPiM     = gPDG->GetParticle( 211);
TParticlePDG* gKaonP   = gPDG->GetParticle( 321);
TParticlePDG* gLambda  = gPDG->GetParticle(2114);
TParticlePDG* gCascade = gPDG->GetParticle(3324);
TParticlePDG* gHDLL    = new TParticlePDG("HDLL","HDLL",2.250,kFALSE, 0, 0,"Baryon",90003,-90003,0);
TParticlePDG* gHDPC    = new TParticlePDG("HDPC","HDPC",2.265,kFALSE, 0, 0,"Baryon",90004,-90004,0);
