// -*- C++ -*-
//
// Package:    HPStracks/HPStracksProducer
// Class:      HPStracksProducer
// 
/**\class HPStracksProducer HPStracksProducer.cc HPStracks/HPStracksProducer/plugins/HPStracksProducer.cc

 Description: create a Track collection from pions of HPS that can be processed in V0 producer

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hlushchenko Olena
//         Created:  Tue, 21 Mar 2017 11:15:40 GMT
//
//


// system include files
#include <memory>
#include <vector> 

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"


////
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

// #include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/AssociationMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/RecoTauPiZero.h"
#include "DataFormats/TauReco/interface/RecoTauPiZeroFwd.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// My new map
#include "DataFormats/HPStracks/interface/TracksToHPSPionsMap.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
// #include "PhysicsTools/HepMCCandAlgos/plugins/MCTruthDeltaRMatcherNew.cc"

// #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
// #include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

// #include "TrackingTools/TransientTrack/interface/TransientTrack.h"

// #include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
// #include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"

// #include "TFile.h"
// #include "TH1D.h"
// #include "TH2D.h"
// #include "TLorentzVector.h"
// #include "TString.h"

//
// class declaration
//

class HPStracksProducer : public edm::stream::EDProducer<> 
{
   public:
      explicit HPStracksProducer(const edm::ParameterSet&);
      ~HPStracksProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      edm::EDGetTokenT<reco::RecoTauPiZeroCollection> TauPiZeroCollectionToken_;//hpsPFTauProducer
      edm::EDGetTokenT<reco::PFTauCollection> TauHPSCollectionToken_;
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
// see: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookChapter9
  // typedef edm::AssociationMap< edm::OneToValue <reco::PFCandidateCollection,//reco::TrackCollection, 
  // int//reco::PFCandidateCollection
  // > > TracksToHPSPionsMap; //CVal -> CKey
  //WORKS1  
    typedef edm::AssociationMap< edm::OneToValue <reco::TrackCollection, int > > TracksToHPSPionsMap_local_int; 
      // replace by DataFormats TracksToHPSPionsMap 
    typedef edm::AssociationMap< edm::OneToValue <reco::TrackCollection, reco::PFCandidatePtr > > TracksToHPSPionsMap_local; 
//
// static data member definitions
//

//
// constructors and destructor
//
HPStracksProducer::HPStracksProducer(const edm::ParameterSet& iConfig)
{
  //register your products
    /* Examples
       produces<ExampleData2>();

       //if do put with a label
       produces<ExampleData2>("label");
     
       //if you want to put into the Run
       produces<ExampleData2,InRun>();
    */
  //now do what ever other initialization is needed
  TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));
  //patTauIDTokens_ = edm::vector_transform(tauIDSrcs_, [this](NameTag const & tag){return mayConsume<pat::PATTauDiscriminator>(tag.second);});

  //WAS std::vector<reco::Track>  
    produces<reco::TrackCollection>("HPSTracks"); //vector<reco::Track>                   "generalTracks"             ""                "RECO"  
    //////produces<reco::TracksToHPSPionsMap>("reco::TracksToHPSPionsMap");
}


HPStracksProducer::~HPStracksProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HPStracksProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::unique_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(std::move(pOut));
  */

  /* this is an EventSetup example
     //Read SetupData from the SetupRecord in the EventSetup
     ESHandle<SetupData> pSetup;
     iSetup.get<SetupRecord>().get(pSetup);
  */ 
  //std::cout << "Begin" << std::endl;
  using namespace edm;

  edm::Handle<reco::PFTauCollection> PF_hps_taus;// typedef vector< PFTau >  PFTauCollection - the collection of charged pions will be taken from here
  iEvent.getByToken( TauHPSCollectionToken_, PF_hps_taus); //TauHPSCollectionToken_ = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer","","RECO"));

  //std::vector<reco::Track> * outputTaus = new reco::TrackCollection();
  std::auto_ptr<reco::TrackCollection> outputTaus(new reco::TrackCollection());
  //reco::TrackCollection * outputTaus = new reco::TrackCollection();
  //outputTaus->reserve(PF_hps_taus->size());

  std::auto_ptr<reco::TracksToHPSPionsMap> tracksToPionsMap(new reco::TracksToHPSPionsMap());
  std::cout<<"try to ini"<<std::endl;
  std::auto_ptr<TracksToHPSPionsMap_local> tracksToPionsMap_local(new TracksToHPSPionsMap_local());
  std::cout<<"ini:"<< tracksToPionsMap_local.get() <<std::endl;
  std::auto_ptr<TracksToHPSPionsMap_local_int> tracksToPionsMap_local_int(new TracksToHPSPionsMap_local_int());



  int tau_index = 0;
  for(reco::PFTauCollection::const_iterator inputTau = PF_hps_taus->begin(); inputTau != PF_hps_taus->end(); ++inputTau, tau_index++)
  {
    reco::PFTau outputTau(*inputTau);
    reco::PFTauRef pftauref(PF_hps_taus, tau_index);

    //pi+-  typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
      //vector < reco::PFCandidatePtr  > 
      std::vector < reco::PFCandidatePtr  >  tau_signalPFChargedHadrCands    = pftauref->signalPFChargedHadrCands();//typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
      std::vector < reco::PFCandidatePtr  >  tau_isolationPFChargedHadrCands = pftauref->isolationPFChargedHadrCands(); // //typedef edm::Ptr<PFCandidate> reco::PFCandidatePtr
    // All pi+-'s 
    std::vector < reco::PFCandidatePtr > tau_picharge = pftauref->signalPFChargedHadrCands(); 
                                          tau_picharge.insert(tau_picharge.end(), tau_isolationPFChargedHadrCands.begin(), tau_isolationPFChargedHadrCands.end());

    for(unsigned pi_i = 0 ; pi_i < tau_picharge.size() ; pi_i++) 
    {
      if ((tau_picharge[pi_i]->trackRef()).isNonnull()) 
      {
          reco::Track track_pion = *(tau_picharge[pi_i]->trackRef());
          reco::TrackRef track_pionRef = tau_picharge[pi_i]->trackRef();// typedef edm::Ref<TrackCollection> reco::TrackRef //typedef std::vector<Track> reco::TrackCollection
          const reco::Track * track = track_pionRef.get();
          /*
            // std::cout << "Pushed pion " << std::endl;
            // std::cout << "Pushed pion "<< tau_index <<":" << pi_i <<  std::endl;
            // std::cout << " Test:" << track->dxy()  << std::endl;// not All the members are accessible, some are missing because no TrackExtra is stored in AOD
            // //std::cout << *(tau_picharge[pi_i]->trackRef()) << std::endl;
            // //outputTaus->push_back(*(tau_picharge[pi_i]->trackRef()));
          */
          outputTaus->push_back(*track);
          // tracksToPionsMap_local_int->insert( tau_picharge[pi_i]->trackRef(), int(1));
          // tracksToPionsMap->insert( tau_picharge[pi_i]->trackRef(),  (tau_picharge[pi_i]) );

          std::cout<<"try to push: pi_i"<< pi_i<< " " << &track_pion << " " << (tau_picharge[pi_i]);// << std::endl;
           //can't insert 
          tracksToPionsMap_local->insert( tau_picharge[pi_i]->trackRef(),  (tau_picharge[pi_i]) );
          //edm::Ref< reco::TrackCollection >( (tau_picharge[pi_i]->trackRef()).get()),
          //edm::Ref< reco::PFCandidateCollection >( tau_picharge[pi_i] ) 

            //tau_picharge[pi_i]->trackRef(), tau_picharge, pi_i);
          
      }
      else std::cout <<"TRACK NOT VALID "<< tau_index <<":" << pi_i << std::endl;
    }
      // for (reco::PFCandidateRefVector::const_iterator iter = tau_picharge.begin(); iter != tau_picharge.end(); ++iter) 
      // {
      //    if (iter->get()->trackRef().isNonnull())
      //     outputTaus.push_back(reco::Track(iter->get()->trackRef()));

      // }
  }


  // Get Non-Tau tracks
    // reco::TrackCollection nonTauTracks;
    //   // remove tau tracks and only tracks associated with the vertex
    //   unsigned int idx = 0;
    //   for (reco::TrackCollection::const_iterator iTrk = trackCollection->begin(); iTrk != trackCollection->end(); ++iTrk, idx++) 
    //   {
    //     reco::TrackRef tmpRef(trackCollection, idx);
    //     reco::TrackRef tmpRefForBase = tmpRef;
        
    //     nonTauTracks.push_back(*iTrk);
    //   }
    
    ////////////////////////

  // std::cout <<"General check"<< std::endl;
  // //for(std::vector<reco::Track>::iterator itv = outputTaus->begin(); itv != outputTaus->end(); ++itv) 
  // for (reco::TrackCollection::const_iterator itv = outputTaus->begin(); itv != outputTaus->end(); ++itv) 
  // {
  //    std::cout<< (*itv).dxy() << "+-" << (*itv).dxyError() << std::endl;;
  // }
  iEvent.put(outputTaus, "HPSTracks");
  //////iEvent.put(tracksToPionsMap, "TracksToHPSPionsMap");
  //std::cout << "End" << std::endl;
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
HPStracksProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
HPStracksProducer::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
HPStracksProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
HPStracksProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HPStracksProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HPStracksProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HPStracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HPStracksProducer);
