// framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// data product includes
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h" // lar::providerFrom<>()
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Wire.h"

// root includes
#include "TFile.h"
#include "TTree.h"

// c++ includes
#include <vector>
#include <iterator>
#include <typeinfo>
#include <memory>
#include <string>
#include <algorithm>
#include <cmath>
#include <map>

// local includes
#include "DataFormat/EventImage2D.h"
#include "DataFormat/EventROI.h"
#include "DataFormat/IOManager.h"

namespace larhsn {

class LArCVMaker : public art::EDAnalyzer {
public:
  explicit LArCVMaker(fhicl::ParameterSet const & pset);
  void analyze(art::Event const & evt);
  void beginJob();
  void endJob();

private:
  void ClearData();
  void MissingWiresList(std::vector<int> missingWiresList[]);
  void ResetROI();
  void SetROISize();
  bool IsInsideTPC(double x, double y, double z);
  void AssignOrigin(const simb::MCParticle& part, double xyz_origin[]);
  void OriginToChannels(double xyz_origin[], int channels_origin[]);
  void OriginToTicks(double x, double ticks_origin[]);
  void FitBoxInDetector(int centerTick[], int centerChannel[], int realTopEdgeTick[], int realBottomEdgeTick[], int realRightEdgeChannel[], int realLeftEdgeChannel[], int tickOptWidth, int channelOptWidth);
  int FindBestAPA(std::vector<int> apas);
  int FindROI(int apa, int plane);

  larcv::IOManager fMgr;

  std::string fWireModuleLabel;
  std::string fMcTruthModuleLabel;
  int fMaxTick;
  int fADCCut;
  int fEventType;

  const int fNumberChannels[3] = { 2400, 2400, 3456 };
  const int fFirstChannel[3] = { 0, 2400, 4800 };
  const int fLastChannel[3] = { 2399, 4799, 8255 };
  const double minTpcBound[3] = { 10., -105.53, 10.1 };
  const double maxTpcBound[3] = { 246.35, 107.47, 1026.9 };

  const int fImageSizeX;
  const int fImageSizeY;

  int fEvent;
  int fAPA;
  int fNumberWires;
  int fNumberTicks;

  bool fOriginIsInsideTPC = true;
  const bool fConstrainProportions;


  std::map<int,std::vector<float> > fWireMap;
  std::vector<std::vector<float> > fImage;

  geo::GeometryCore const*             fGeometry;           ///< pointer to the Geometry service
  detinfo::DetectorProperties const* fDetectorProperties; ///< Pointer to the detector properties
}; // EOClass LArCVMaker

LArCVMaker::LArCVMaker(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset),
    fMgr(larcv::IOManager::kWRITE),
    fWireModuleLabel(pset.get<std::string>("WireModuleLabel")),
    fMcTruthModuleLabel(pset.get<std::string>("McTruthModuleLabel")),
    fMaxTick(pset.get<int>("MaxTick")),
    fADCCut(pset.get<int>("ADCCut")),
    fEventType(pset.get<int>("EventType")),
    fImageSizeX(pset.get<int>("ImageSizeX")),
    fImageSizeY(pset.get<int>("ImageSizeY")),
    fConstrainProportions(pset.get<bool>("ConstrainProportions"))
{} // EOFunction LArCVMaker::LArCVMaker

void LArCVMaker::beginJob() {
  fGeometry = lar::providerFrom<geo::Geometry>();
  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();
  std::string filename;
  if (std::getenv("PROCESS") != nullptr) filename = "larcv_" + std::string(std::getenv("PROCESS")) + ".root";
  else filename = "larcv.root";
  fMgr.set_out_file(filename);
  fMgr.initialize();
} // EOFunction LArCVMaker::beginJob

void LArCVMaker::endJob() {
  fMgr.finalize();
} // EOFunction LArCVMaker::endJob

void LArCVMaker::ClearData() {
  fWireMap.clear();
} // EOFunction LArCVMaker::ClearData

void LArCVMaker::MissingWiresList(std::vector<int> missingWiresList[]) {
  for (int it_plane = 0; it_plane < 3; ++it_plane) {
    missingWiresList[it_plane].clear();
    for (int it_channel = fFirstChannel[it_plane]; it_channel < fLastChannel[it_plane]; ++it_channel) {
      if (fWireMap.find(it_channel) != fWireMap.end()){
        std::cout << "it_channel: " << it_channel << " is ok" << std::endl;
      }
      else missingWiresList[it_plane].push_back(it_channel);
    } 
  } 
  return;
}

bool LArCVMaker::IsInsideTPC(double x, double y, double z) {
  // Check whether coordinates are inside TPC 
  bool isInsideX = (x>minTpcBound[0] && x<maxTpcBound[0]);
  bool isInsideY = (y>minTpcBound[1] && y<maxTpcBound[1]);
  bool isInsideZ = (z>minTpcBound[2] && z<maxTpcBound[2]);
  if (isInsideX && isInsideY && isInsideZ) return true;
  else return false;
}

void LArCVMaker::AssignOrigin(const simb::MCParticle& part, double xyz_origin[]) {
  // Assign part coordinates to xyz vector and move to nearest edge
  // if particle is outside TPC
  if (part.Vx() < minTpcBound[0]) xyz_origin[0] = minTpcBound[0];
  else if (part.Vx() > maxTpcBound[0]) xyz_origin[0] = maxTpcBound[0];
  else xyz_origin[0] = part.Vx();

  if (part.Vy() < minTpcBound[1]) xyz_origin[1] = minTpcBound[1];
  else if (part.Vy() > maxTpcBound[1]) xyz_origin[1] = maxTpcBound[1];
  else xyz_origin[1] = part.Vy();

  if (part.Vz() < minTpcBound[2]) xyz_origin[2] = minTpcBound[2];
  else if (part.Vz() > maxTpcBound[2]) xyz_origin[2] = maxTpcBound[2];
  else xyz_origin[2] = part.Vz();
  return;
}

void LArCVMaker::OriginToChannels(double xyz_origin[], int channels_origin[]) {
  // Find nearest wires to provided coordinates
  raw::ChannelID_t channels_info_0 = fGeometry->NearestChannel(xyz_origin,0);
  raw::ChannelID_t channels_info_1 = fGeometry->NearestChannel(xyz_origin,1);
  raw::ChannelID_t channels_info_2 = fGeometry->NearestChannel(xyz_origin,2);
  channels_origin[0] = channels_info_0;
  channels_origin[1] = channels_info_1;
  channels_origin[2] = channels_info_2;
  return;
}

void LArCVMaker::OriginToTicks(double x, double ticks_origin[]) {
  // Find nearest ticks to provided coordinates
  ticks_origin[0] = fDetectorProperties->ConvertXToTicks(x, 0, 0, 0);
  ticks_origin[1] = fDetectorProperties->ConvertXToTicks(x, 1, 0, 0);
  ticks_origin[2] = fDetectorProperties->ConvertXToTicks(x, 2, 0, 0);
  return;
}

void LArCVMaker::FitBoxInDetector(int centerTick[], int centerChannel[], int realTopEdgeTick[], int realBottomEdgeTick[], int realRightEdgeChannel[], int realLeftEdgeChannel[], int tickOptWidth, int channelOptWidth) {
  // Solve awkward problem of box lying somehow outside of detector
  // Takes centerTick and centerChannel and try to build a box around those coordinates with tickOptWidth and channelOptWidth dimensions. If that box ends up intersecting the TPC boundaries than it shifts it until if fits inside the TPC. It then returns the coordinates of the edges of the box.
  for (int it_plane = 0; it_plane < 3; ++it_plane) {

    // Tick resizing
    if ((centerTick[it_plane]-tickOptWidth/2) < 0) realBottomEdgeTick[it_plane] = 0;
    else realBottomEdgeTick[it_plane] = centerTick[it_plane]-tickOptWidth/2;

    if ((centerTick[it_plane]+tickOptWidth/2) > fMaxTick){
      realTopEdgeTick[it_plane] = fMaxTick;
      realBottomEdgeTick[it_plane] = realTopEdgeTick[it_plane] - tickOptWidth;
    }
    else realTopEdgeTick[it_plane] = realBottomEdgeTick[it_plane] + tickOptWidth;

    // Channel resizing
    if ((centerChannel[it_plane]-channelOptWidth/2) < fFirstChannel[it_plane]) realLeftEdgeChannel[it_plane] = fFirstChannel[it_plane];
    else realLeftEdgeChannel[it_plane] = centerChannel[it_plane]-channelOptWidth/2;

    if ((centerChannel[it_plane]+channelOptWidth/2) > fLastChannel[it_plane]){
     realRightEdgeChannel[it_plane] = fLastChannel[it_plane];
     realLeftEdgeChannel[it_plane] = realRightEdgeChannel[it_plane] - channelOptWidth;
    }
    else realRightEdgeChannel[it_plane] = realLeftEdgeChannel[it_plane] + channelOptWidth;
  }
  return;
}


void LArCVMaker::analyze(art::Event const & evt) {
  // Clear data before starting analyzing event
  ClearData();
  bool isEmptyEvent = false;

  // Get event number
  fEvent = evt.event();
  // Set larcv manager
  // std::cout << "Generating manager for event:" << std::endl;
  // std::cout <<  "EVT: " << evt.id().run() << " " << evt.id().subRun() << " " << evt.id().event() << std::endl;
  fMgr.set_id(evt.id().run(),evt.id().subRun(),evt.id().event());
  // std::cout <<  "FMGR 1: " << fMgr.event_id().run() << " " << fMgr.event_id().subrun() << " " << fMgr.event_id().event() << std::endl;
  // Get objects from event
  art::Handle<std::vector<recob::Wire>> wireHandle;
  evt.getByLabel(fWireModuleLabel,wireHandle);

  // Code is built on assumption that there's only 1 mcTruth for event.
  // Return error if that's not the case
  // if ((*mctruthHandle).size()!=1){
  //   std::ostringstream oss;
  //   oss << "Unexpected number of mcTruth objects (" << (*mctruthHandle).size() << ") in event " << fEvent << std::endl;
  //   std::string ss = oss.str();
  //   throw std::invalid_argument(ss);
  // }

  // Take first MC particle (make sure it's a neutrino) and take coordinates
  art::Handle<std::vector<simb::MCTruth>> mctruthHandle;
  std::vector<art::Ptr<simb::MCTruth>> TruthList;
  if (evt.getByLabel(fMcTruthModuleLabel,mctruthHandle)) art::fill_ptr_vector(TruthList,mctruthHandle);
  art::Ptr<simb::MCTruth> mctruth = TruthList[0];
  const simb::MCParticle& part = (*mctruth).GetParticle(0);

  // Report error if not neutrino
  // if (part.PdgCode()!=14){
  //   std::ostringstream oss;
  //   oss << "First McParticle not a neutrino (PDG: " << part.PdgCode() << ") in event " << fEvent << std::endl;
  //   std::string ss = oss.str();
  //   throw std::invalid_argument(ss);
  // }

  // Store coordinates in vector from first particle and convert them
  // to wire and time tick coordinates.
  double xyz_origin[3];
  double ticks_origin[3];
  int channels_origin[3];

  // double 
  // fOriginIsInsideTPC = LArCVMaker::IsInsideTPC(part.Vx(),part.Vy(),part.Vz());
  // if (!fOriginIsInsideTPC){
  //   std::cout << "Event originates outside TPC!" << std::endl;
  //   std::cout << "Using nearest edge as origin point." << std::endl;
  // }
  LArCVMaker::AssignOrigin(part, xyz_origin);
  // std::cout << "Particle origin: (" << part.Vx() << ", " << part.Vy() << ", " << part.Vz() << ")" << std::endl;
  // std::cout << "Nearest edge origin: (" << xyz_origin[0] << ", " << xyz_origin[1] << ", " << xyz_origin[2] << ")" << std::endl;
  LArCVMaker::OriginToChannels(xyz_origin, channels_origin);
  // std::cout << "Nearest origin channels: (" << channels_origin[0] << ", " << channels_origin[1] << ", " << channels_origin[2] << ")" << std::endl;
  LArCVMaker::OriginToTicks(xyz_origin[0], ticks_origin);
  // std::cout << "Timeticks: (" << ticks_origin[0] << ", " << ticks_origin[1] << ", " << ticks_origin[2] << ")" << std::endl;

  // Loop over each wire and add it to the wire map
  for (std::vector<recob::Wire>::const_iterator it = wireHandle->begin(); it != wireHandle->end(); ++it) {
    const recob::Wire & wire = *it;
    fWireMap.insert(std::pair<int,std::vector<float> >(wire.Channel(),std::vector<float>(wire.Signal())));
  }

  // Find ROI starting from origin point
  int topEdgeTick[3], bottomEdgeTick[3];
  int rightEdgeChannel[3], leftEdgeChannel[3];
  int centerTick[3], centerChannel[3];
  int tickWidth[3], channelWidth[3];
  int maxTickWidth, maxChannelWidth;

  // Plane loop
  for (int it_plane = 0; it_plane < 3; ++it_plane) {
    topEdgeTick[it_plane] = ticks_origin[it_plane];
    bottomEdgeTick[it_plane] = ticks_origin[it_plane];
    rightEdgeChannel[it_plane] = channels_origin[it_plane];
    leftEdgeChannel[it_plane] = channels_origin[it_plane];

    // Diagnostic message
    // std::cout << std::endl;
    // std::cout << "PLANE " << it_plane << std::endl;
    // std::cout << "Origin: (" << channels_origin[it_plane] << ", " << ticks_origin[it_plane] << ")" << std::endl;

    // IncreasingChannel loop
    for (int it_channel = channels_origin[it_plane]; it_channel < fLastChannel[it_plane]; ++it_channel) {
      // IncreasingChannel_IncreasingTick loop
      for (int it_tick = ticks_origin[it_plane]; it_tick < fMaxTick; ++it_tick) {
        double pixel = 0.;
        if (fWireMap.find(it_channel) != fWireMap.end()){
          pixel = fWireMap[it_channel][it_tick];
        }
        if (pixel>fADCCut && it_tick>topEdgeTick[it_plane]) topEdgeTick[it_plane] = it_tick;
        if (pixel>fADCCut && it_channel>rightEdgeChannel[it_plane]) rightEdgeChannel[it_plane] = it_channel;
      } // EOIncreasingChannel_IncreasingTick loop

      // IncreasingChannel_DecreasingTick loop
      for (int it_tick = ticks_origin[it_plane]; it_tick > 0; --it_tick) {
        double pixel = 0.;
        if (fWireMap.find(it_channel) != fWireMap.end()){
          pixel = fWireMap[it_channel][it_tick];
        }
        if (pixel>fADCCut && it_tick<bottomEdgeTick[it_plane]) bottomEdgeTick[it_plane] = it_tick;
        if (pixel>fADCCut && it_channel>rightEdgeChannel[it_plane]) rightEdgeChannel[it_plane] = it_channel;
      } // EOIncreasingChannel_DecreasingTick loop
    } //EOIncreasingChannel loop

    // DecreasingChannel loop
    for (int it_channel = channels_origin[it_plane]; it_channel > fFirstChannel[it_plane]; --it_channel) {
      // DecreasingChannel_IncreasingTick loop
      for (int it_tick = ticks_origin[it_plane]; it_tick < fMaxTick; ++it_tick) {
        double pixel = 0.;
        if (fWireMap.find(it_channel) != fWireMap.end()){
          pixel = fWireMap[it_channel][it_tick];
        }
        if (pixel>fADCCut && it_tick>topEdgeTick[it_plane]) topEdgeTick[it_plane] = it_tick;
        if (pixel>fADCCut && it_channel<leftEdgeChannel[it_plane]) leftEdgeChannel[it_plane] = it_channel;
      } // EODecreasingChannel_IncreasingTick loop

      // Decreasingchannels_DecreasingTick loop
      for (int it_tick = ticks_origin[it_plane]; it_tick > 0; --it_tick) {
        double pixel = 0.;
        if (fWireMap.find(it_channel) != fWireMap.end()){
          pixel = fWireMap[it_channel][it_tick];
        }
        if (pixel>fADCCut && it_tick<bottomEdgeTick[it_plane]) bottomEdgeTick[it_plane] = it_tick;
        if (pixel>fADCCut && it_channel<leftEdgeChannel[it_plane]) leftEdgeChannel[it_plane] = it_channel;
      } // EODecreasingChannel_DecreasingTick loop
    } //EODecreasingChannel loop

    // Diagnostic message
    // std::cout << "LRChannels: [" << leftEdgeChannel[it_plane] << ", " << rightEdgeChannel[it_plane] << "]" << std::endl;
    // std::cout << "UDTicks: [" << bottomEdgeTick[it_plane] << ", " << topEdgeTick[it_plane] << "]" << std::endl;

    tickWidth[it_plane] = topEdgeTick[it_plane] - bottomEdgeTick[it_plane];
    channelWidth[it_plane] = rightEdgeChannel[it_plane] - leftEdgeChannel[it_plane];
    centerTick[it_plane] = (topEdgeTick[it_plane] + bottomEdgeTick[it_plane])/2;
    centerChannel[it_plane] = (rightEdgeChannel[it_plane] + leftEdgeChannel[it_plane])/2;

    // Mark event as empty if there is no image in it, so it doesn't get recorded.
    if (channelWidth[it_plane]==0 || tickWidth[it_plane]==0) isEmptyEvent = true;

    // Diagnostic message
    // std::cout << "Minimum size ROI dimension: " << channelWidth[it_plane] << "x" << tickWidth[it_plane] << std::endl;
  } //EOPlane loop

  // Determine optimal box size for all planes.
  // All planes should have same size, size must also be multiple of setting size
  maxTickWidth = std::max({tickWidth[0],tickWidth[1],tickWidth[2]});
  maxChannelWidth = std::max({channelWidth[0],channelWidth[1],channelWidth[2]});
  // Find now the maximum multiple across all channels (must also be same for x and y if you want to keep aspect ratio)
  int maxTickMult = ceil(maxTickWidth/float(fImageSizeY));
  int maxChannelMult = ceil(maxChannelWidth/float(fImageSizeX));
  int maxTotMult = std::max({maxTickMult,maxChannelMult});
  if (fConstrainProportions) {
    maxTickMult = maxTotMult;
    maxChannelMult = maxTotMult;
  }

  int tickOptWidth = fImageSizeY*maxTickMult;
  int channelOptWidth = fImageSizeX*maxChannelMult;

  if (tickOptWidth!=0 && channelOptWidth!=0){
    printf("RESOLUTIONS:\n");
    for(int i=0;i!=3;i++) {
      printf("Plane %i : [%ix%i]\n", i, channelWidth[i], tickWidth[i]);
    }
    printf("\n");
    printf("Original common resolution: [%ix%i]\n", maxChannelWidth, maxTickWidth);
    printf("Expanded resolution: [%ix%i]\n", channelOptWidth, tickOptWidth);
    printf("Compression factor: [%i:%i]\n", maxChannelMult, maxTickMult);
    printf("Final resolution: [%ix%i]\n", fImageSizeX, fImageSizeY);
    printf("\n");
  }

  //Diagnostic message
  // std::cout << std::endl;
  // std::cout << "Optimal dimensions: (" << channelOptWidth << "x" << tickOptWidth << ")" << std::endl;

  // Solve awkward problems with optimal box lying outside the detector
  int realTopEdgeTick[3], realBottomEdgeTick[3];
  int realRightEdgeChannel[3], realLeftEdgeChannel[3];
  LArCVMaker::FitBoxInDetector(centerTick, centerChannel, realTopEdgeTick, realBottomEdgeTick, realRightEdgeChannel, realLeftEdgeChannel, tickOptWidth, channelOptWidth);

  // Diagnostic message
  // std::cout << std::endl << std::endl;
  // std::cout << "ROI RESIZED" << std::endl << std::endl;
  // for (int it_plane = 0; it_plane < 3; ++it_plane) {
  //   std::cout << std::endl << "PLANE " << it_plane << std::endl;
  //   std::cout << "Origin: (" << centerChannel[it_plane] << ", " << centerTick[it_plane] << ")" << std::endl;
  //   std::cout << "LRChannels: [" << realLeftEdgeChannel[it_plane] << ", " << realRightEdgeChannel[it_plane] << "]" << std::endl;
  //   std::cout << "UDTicks: [" << realBottomEdgeTick[it_plane] << ", " << realTopEdgeTick[it_plane] << "]" << std::endl;
  //   std::cout << "ROI dimension: " << realRightEdgeChannel[it_plane]-realLeftEdgeChannel[it_plane] << "x" << realTopEdgeTick[it_plane]-realBottomEdgeTick[it_plane] << std::endl;
  // }

  // Get handle on larcv image
  if (isEmptyEvent==false){
    auto images = (larcv::EventImage2D*)(fMgr.get_data(larcv::kProductImage2D, "tpc"));

    // Create images from the wire map
    for (int it_plane = 0; it_plane < 3; ++it_plane) {
      larcv::Image2D image(tickOptWidth,channelOptWidth);
      for (int it_channel = realLeftEdgeChannel[it_plane]; it_channel < realRightEdgeChannel[it_plane]; ++it_channel) {
        for (int it_tick = realBottomEdgeTick[it_plane]; it_tick < realTopEdgeTick[it_plane]; ++it_tick) {
          int xPixel = it_channel - realLeftEdgeChannel[it_plane];
          int yPixel = it_tick - realBottomEdgeTick[it_plane];
          if (fWireMap.find(it_channel) != fWireMap.end()) image.set_pixel(yPixel,xPixel,fWireMap[it_channel][it_tick]);
          else image.set_pixel(yPixel,xPixel,0);
        }
      }
      image.compress(fImageSizeY,fImageSizeX);
      images->Emplace(std::move(image));
    }

    auto roi = (larcv::EventROI*)(fMgr.get_data(larcv::kProductROI, "tpc"));
    roi->Emplace(larcv::ROI((larcv::ROIType_t)fEventType));

    fMgr.save_entry();
  }
  else {
    std::cout <<  "Empty image not saved [" << fMgr.event_id().run() << ", " << fMgr.event_id().subrun() << ", " << fMgr.event_id().event() << "]" << std::endl;
  }
} //EOFunction LArCVMaker::analyze

DEFINE_ART_MODULE(LArCVMaker)

} // namespace larhsn

