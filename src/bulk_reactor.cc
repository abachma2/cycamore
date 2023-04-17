#include "bulk_reactor.h"

using cyclus::Material;
using cyclus::Composition;
using cyclus::toolkit::ResBuf;
using cyclus::toolkit::MatVec;
using cyclus::KeyError;
using cyclus::ValueError;
using cyclus::Request;

namespace cycamore {

BulkReactor::BulkReactor(cyclus::Context* ctx)
    : cyclus::Facility(ctx),
      core_mass(0),
      discharge_mass(0),
      spent_mass(0),
      cycle_time(0),
      refuel_time(0),
      cycle_step(0),
      power_cap(0),
      power_name("power"),
      discharged(false),
      latitude(0.0),
      longitude(0.0),
      coordinates(latitude, longitude) {}


#pragma cyclus def clone cycamore::BulkReactor

#pragma cyclus def schema cycamore::BulkReactor

#pragma cyclus def annotations cycamore::BulkReactor

#pragma cyclus def infiletodb cycamore::BulkReactor

#pragma cyclus def snapshot cycamore::BulkReactor

#pragma cyclus def snapshotinv cycamore::BulkReactor

#pragma cyclus def initinv cycamore::BulkReactor

void BulkReactor::InitFrom(BulkReactor* m) {
  #pragma cyclus impl initfromcopy cycamore::BulkReactor
  cyclus::toolkit::CommodityProducer::Copy(m);
}

void BulkReactor::InitFrom(cyclus::QueryableBackend* b) {
  #pragma cyclus impl initfromdb cycamore::BulkReactor

  namespace tk = cyclus::toolkit;
  tk::CommodityProducer::Add(tk::Commodity(power_name),
                             tk::CommodInfo(power_cap, power_cap));

  for (int i = 0; i < side_products.size(); i++) {
    tk::CommodityProducer::Add(tk::Commodity(side_products[i]),
                               tk::CommodInfo(side_product_quantity[i],
                                              side_product_quantity[i]));
  }
}

void BulkReactor::EnterNotify() {
  cyclus::Facility::EnterNotify();


  // If the user ommitted fuel_prefs, we set it to zeros for each fuel
  // type.  Without this segfaults could occur - yuck.
  if (fuel_prefs.size() == 0) {
    for (int i = 0; i < fuel_outcommods.size(); i++) {
      fuel_prefs.push_back(cyclus::kDefaultPref);
    }
  }

  // Test if any side products have been defined.
  if (side_products.size() == 0){
    hybrid_ = false;
  }

  // input consistency checking:
  int n = recipe_change_times.size();
  std::stringstream ss;
  if (recipe_change_commods.size() != n) {
    ss << "prototype '" << prototype() << "' has "
       << recipe_change_commods.size()
       << " recipe_change_commods vals, expected " << n << "\n";
  }
  if (recipe_change_in.size() != n) {
    ss << "prototype '" << prototype() << "' has " << recipe_change_in.size()
       << " recipe_change_in vals, expected " << n << "\n";
  }
  if (recipe_change_out.size() != n) {
    ss << "prototype '" << prototype() << "' has " << recipe_change_out.size()
       << " recipe_change_out vals, expected " << n << "\n";
  }

  n = pref_change_times.size();
  if (pref_change_commods.size() != n) {
    ss << "prototype '" << prototype() << "' has " << pref_change_commods.size()
       << " pref_change_commods vals, expected " << n << "\n";
  }
  if (pref_change_values.size() != n) {
    ss << "prototype '" << prototype() << "' has " << pref_change_values.size()
       << " pref_change_values vals, expected " << n << "\n";
  }

  if (ss.str().size() > 0) {
    throw cyclus::ValueError(ss.str());
  }
  RecordPosition();
}

bool BulkReactor::CheckDecommissionCondition() {
  return core.count() == 0 && spent.count() == 0;
}

void BulkReactor::Tick() {
  // The following code must go in the Tick so they fire on the time step
  // following the cycle_step update - allowing for the all reactor events to
  // occur and be recorded on the "beginning" of a time step.  Another reason
  // they
  // can't go at the beginnin of the Tock is so that resource exchange has a
  // chance to occur after the discharge on this same time step.
  if (retired()) {
    Record("RETIRED", "");

    if (context()->time() == exit_time() + 1) { // only need to transmute once
      if (decom_transmute_all == true) {
        Discharge_Transmute(core.quantity());
      }
      else {
        Discharge_Transmute(core.quantity(), core.quantity()/2);
      }
    }
    // in case a cycle lands exactly on our last time step, we will need to
    // burn a batch from fresh inventory on this time step.  When retired,
    // this batch also needs to be discharged to spent fuel inventory.
    if(CheckDecommissionCondition()) {
      Decommission();
    }
    return;
  }

  if (cycle_step == cycle_time) {
    Record("CYCLE_END", "");
  }

  if (cycle_step >= cycle_time && !discharged) {
    discharged = Discharge_Transmute(discharge_mass);
  }
  if (cycle_step >= cycle_time) {
    Load();
  }

  int t = context()->time();

  // update preferences
  for (int i = 0; i < pref_change_times.size(); i++) {
    int change_t = pref_change_times[i];
    if (t != change_t) {
      continue;
    }

    std::string incommod = pref_change_commods[i];
    for (int j = 0; j < fuel_incommods.size(); j++) {
      if (fuel_incommods[j] == incommod) {
        fuel_prefs[j] = pref_change_values[i];
        break;
      }
    }
  }

  // update recipes
  for (int i = 0; i < recipe_change_times.size(); i++) {
    int change_t = recipe_change_times[i];
    if (t != change_t) {
      continue;
    }

    std::string incommod = recipe_change_commods[i];
    for (int j = 0; j < fuel_incommods.size(); j++) {
      if (fuel_incommods[j] == incommod) {
        fuel_inrecipes[j] = recipe_change_in[i];
        fuel_outrecipes[j] = recipe_change_out[i];
        break;
      }
    }
  }
}

std::set<cyclus::RequestPortfolio<Material>::Ptr> BulkReactor::GetMatlRequests() {
  using cyclus::RequestPortfolio;

  std::set<RequestPortfolio<Material>::Ptr> ports;
  Material::Ptr m;

  
  // second min expression reduces assembles to amount needed until
  // retirement if it is near.
  double order_mass = core_mass - core.quantity();
  

  if (order_mass == 0) {
    return ports;
  } else if (retired()) {
    return ports;
  }

  RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
  std::vector<Request<Material>*> mreqs;
  for (int j = 0 ; j < fuel_incommods.size(); j++){
    std::string commod = fuel_incommods[j];
    double pref = fuel_prefs[j];
  
    Composition::Ptr recipe = context()->GetRecipe(fuel_inrecipes[j]);
    m = Material::CreateUntracked(order_mass, recipe);
    
    Request<Material>* r = port->AddRequest(m, this, commod, pref, true);
    mreqs.push_back(r);
    }

  std::vector<double>::iterator result;
  result = std::max_element(fuel_prefs.begin(), fuel_prefs.end());
  int max_index = std::distance(fuel_prefs.begin(), result);
  
  cyclus::toolkit::RecordTimeSeries<double>("demand"+fuel_incommods[max_index], this,
                                              order_mass);
  port->AddMutualReqs(mreqs);
  ports.insert(port);
  

  return ports;
}

void BulkReactor::GetMatlTrades(
    const std::vector<cyclus::Trade<Material> >& trades,
    std::vector<std::pair<cyclus::Trade<Material>, Material::Ptr> >&
        responses) {
  using cyclus::Trade;

  
  for (int i = 0; i < trades.size(); i++) {
    double qty = trades[i].amt;
    Material::Ptr m = spent.Pop(qty);
    responses.push_back(std::make_pair(trades[i], m));
    res_indexes.erase(m->obj_id());
  }
}

void BulkReactor::AcceptMatlTrades(const std::vector<
    std::pair<cyclus::Trade<Material>, Material::Ptr> >& responses) {
  std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                        cyclus::Material::Ptr> >::const_iterator trade;
  
  
  std::stringstream ss;

  for (trade = responses.begin(); trade != responses.end(); ++trade) {
    std::string commod = trade->first.request->commodity();
    Material::Ptr m = trade->second;
    index_res(m, commod);

    if (core.quantity() < core_mass) {
      core.Push(m);
    } else {
      fresh.Push(m);
    }
  }
}

std::set<cyclus::BidPortfolio<Material>::Ptr> BulkReactor::GetMatlBids(
    cyclus::CommodMap<Material>::type& commod_requests) {
    using cyclus::BidPortfolio;
  std::set<BidPortfolio<Material>::Ptr> ports;
  
  bool gotmats = false;
  std::map<std::string, MatVec> all_mats;

  if (uniq_outcommods_.empty()) {
    for (int i = 0; i < fuel_outcommods.size(); i++) {
      uniq_outcommods_.insert(fuel_outcommods[i]);
      
    }
  }


  
  std::set<std::string>::iterator it;
  for (it = uniq_outcommods_.begin(); it != uniq_outcommods_.end(); ++it) {
  

    std::string commod = *it;
    std::vector<Request<Material>*>& reqs = commod_requests[commod];
  
    if (reqs.size() == 0) {
      continue;
    } else if (!gotmats) {
      all_mats = PeekSpent();
    }

    MatVec mats = all_mats[commod];
    if (mats.size() == 0) {
      continue;
    }

    BidPortfolio<Material>::Ptr port(new BidPortfolio<Material>());

    for (int j = 0; j < reqs.size(); j++) {
      Request<Material>* req = reqs[j];
      double tot_bid = 0;
      for (int k=0 ; k < mats.size(); k++){
        Material::Ptr m = mats[k];
        tot_bid += m->quantity();
        port->AddBid(req, m, this, true);
        if (tot_bid >= req->target()->quantity()){
            break;
        }
      }
    }

    double tot_qty = 0;
    for (int j=0; j<mats.size(); j++){
        tot_qty += mats[j]->quantity();
    }

    cyclus::CapacityConstraint<Material> cc(tot_qty);
    port->AddConstraint(cc);
    ports.insert(port);
  }

  return ports;
}


std::map<std::string, MatVec> BulkReactor::PeekSpent(){
    std::map<std::string, MatVec> mapped;
    MatVec mats = spent.PopN(spent.count());
    spent.Push(mats);
    for (int i=0; i<mats.size(); i++){
        std::string commod = fuel_outcommod(mats[i]);
        mapped[commod].push_back(mats[i]);
    }
    
    return mapped;
}

void BulkReactor::Tock() {
  
  if (retired()) {
    return;
  }
  
  // Check that irradiation and refueling periods are over, that 
  // the core is full and that fuel was successfully discharged in this refueling time.
  // If this is the case, then a new cycle will be initiated.
  if (cycle_step >= cycle_time + refuel_time && core.quantity() == core_mass && discharged == true) {
    discharged = false;
    cycle_step = 0;
  }

  if (cycle_step == 0 && core.quantity() == core_mass) {
    Record("CYCLE_START", "");
  }

  if (cycle_step >= 0 && cycle_step < cycle_time &&
      core.quantity() == core_mass) {
    cyclus::toolkit::RecordTimeSeries<cyclus::toolkit::POWER>(this, power_cap);
    cyclus::toolkit::RecordTimeSeries<double>("supplyPOWER", this, power_cap);
    RecordSideProduct(true);
  } else {
    cyclus::toolkit::RecordTimeSeries<cyclus::toolkit::POWER>(this, 0);
    cyclus::toolkit::RecordTimeSeries<double>("supplyPOWER", this, 0);
    RecordSideProduct(false);
  }

  // "if" prevents starting cycle after initial deployment until core is full
  // even though cycle_step is its initial zero.
  if (cycle_step > 0 || core.quantity() == core_mass) {
    cycle_step++;
  }
}


bool BulkReactor::Discharge_Transmute(double qty) {

  if (spent.capacity() - spent.quantity() < qty) {
    Record("DISCHARGE", "failed");
    return false;  // not enough room in spent buffer
  }
  Material::Ptr m = core.Pop(qty);
  m->Transmute(context()->GetRecipe(fuel_outrecipe(m)));
  spent.Push(m);
  return true;
}

bool BulkReactor::Discharge_Transmute(double qty, double transmute_qty){

  if (spent.capacity() - spent.quantity() < discharge_mass) {
    Record("DISCHARGE", "failed");
    return false;  // not enough room in spent buffer
  }
  Material::Ptr m = core.Pop(transmute_qty);
  m->Transmute(context()->GetRecipe(fuel_outrecipe(m)));
  Material::Ptr m2 = core.Pop(qty-transmute_qty);
  spent.Push(m);
  spent.Push(m2);
}

void BulkReactor::Load() {
  double load_qty = std::min(core.capacity() - core.quantity(), fresh.quantity());
  if (load_qty == 0) {
    return;
  }
  core.Push(fresh.Pop(load_qty));
}


std::string BulkReactor::fuel_incommod(Material::Ptr m) {
  int i = res_indexes[m->obj_id()];
  if (i >= fuel_incommods.size()) {
    throw KeyError("cycamore::Reactor - no incommod for material object");
  }
  return fuel_incommods[i];
}

std::string BulkReactor::fuel_outcommod(Material::Ptr m) {
  int i = res_indexes[m->obj_id()];
  if (i >= fuel_outcommods.size()) {
    throw KeyError("cycamore::Reactor - no outcommod for material object");
  }
  return fuel_outcommods[i];
}

std::string BulkReactor::fuel_inrecipe(Material::Ptr m) {
  int i = res_indexes[m->obj_id()];
  if (i >= fuel_inrecipes.size()) {
    throw KeyError("cycamore::Reactor - no inrecipe for material object");
  }
  return fuel_inrecipes[i];
}

std::string BulkReactor::fuel_outrecipe(Material::Ptr m) {
  int i = res_indexes[m->obj_id()];
  if (i >= fuel_outrecipes.size()) {
    throw KeyError("cycamore::BulkReactor - no outrecipe for material object");
  }
  return fuel_outrecipes[i];
}

double BulkReactor::fuel_pref(Material::Ptr m) {
  int i = res_indexes[m->obj_id()];
  if (i >= fuel_prefs.size()) {
    return 0;
  }
  return fuel_prefs[i];
}

void BulkReactor::index_res(cyclus::Resource::Ptr m, std::string incommod) {
  for (int i = 0; i < fuel_incommods.size(); i++) {
    if (fuel_incommods[i] == incommod) {
      res_indexes[m->obj_id()] = i;
      return;
    }
  }
  throw ValueError(
      "cycamore::BulkReactor - received unsupported incommod material");
}


void BulkReactor::RecordSideProduct(bool produce){
  if (hybrid_){
    double value;
    for (int i = 0; i < side_products.size(); i++) {
      if (produce){
          value = side_product_quantity[i];
      }
      else {
          value = 0;
      }

      context()
          ->NewDatum("ReactorSideProducts")
          ->AddVal("AgentId", id())
          ->AddVal("Time", context()->time())
          ->AddVal("Product", side_products[i])
          ->AddVal("Value", value)
          ->Record();
    }
  }
}

void BulkReactor::Record(std::string name, std::string val) {
  context()
      ->NewDatum("ReactorEvents")
      ->AddVal("AgentId", id())
      ->AddVal("Time", context()->time())
      ->AddVal("Event", name)
      ->AddVal("Value", val)
      ->Record();
}

void BulkReactor::RecordPosition() {
  std::string specification = this->spec();
  context()
      ->NewDatum("AgentPosition")
      ->AddVal("Spec", specification)
      ->AddVal("Prototype", this->prototype())
      ->AddVal("AgentId", id())
      ->AddVal("Latitude", latitude)
      ->AddVal("Longitude", longitude)
      ->Record();
}

extern "C" cyclus::Agent* ConstructBulkReactor(cyclus::Context* ctx) {
  return new BulkReactor(ctx);
}

}  // namespace cycamore
