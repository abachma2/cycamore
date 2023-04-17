#include <gtest/gtest.h>

#include <sstream>

#include "cyclus.h"

using pyne::nucname::id;
using cyclus::Composition;
using cyclus::Material;
using cyclus::QueryResult;
using cyclus::Cond;
using cyclus::toolkit::MatQuery;

namespace cycamore {
namespace bulk_reactortests {

Composition::Ptr c_uox() {
  cyclus::CompMap m;
  m[id("u235")] = 0.04;
  m[id("u238")] = 0.96;
  return Composition::CreateFromMass(m);
};

Composition::Ptr c_mox() {
  cyclus::CompMap m;
  m[id("u235")] = .7;
  m[id("u238")] = 100;
  m[id("pu239")] = 3.3;
  return Composition::CreateFromMass(m);
};

Composition::Ptr c_spentuox() {
  cyclus::CompMap m;
  m[id("u235")] =  .8;
  m[id("u238")] =  100;
  m[id("pu239")] = 1;
  return Composition::CreateFromMass(m);
};

Composition::Ptr c_spentmox() {
  cyclus::CompMap m;
  m[id("u235")] =  .2;
  m[id("u238")] =  100;
  m[id("pu239")] = .9;
  return Composition::CreateFromMass(m);
};

Composition::Ptr c_water() {
  cyclus::CompMap m;
  m[id("O16")] =  1;
  m[id("H1")] =  2;
  return Composition::CreateFromAtom(m);
};

// Test that with a zero refuel_time and a zero capacity fresh fuel buffer
// (the default), fuel can be ordered and the cycle started with no time step
// delay.
TEST(BulkReactorTests, JustInTimeOrdering) {
  std::string config =
     "  <fuel_inrecipes>  <val>lwr_fresh</val>  </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>lwr_spent</val>  </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>enriched_u</val> </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>      </fuel_outcommods>  "
     "  <fuel_prefs>      <val>1.0</val>        </fuel_prefs>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <core_mass>300</core_mass>  "
     "  <discharge_mass>10</discharge_mass>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("enriched_u").Finalize();
  sim.AddRecipe("lwr_fresh", c_uox());
  sim.AddRecipe("lwr_spent", c_spentuox());
  int id = sim.Run();

  QueryResult qr = sim.db().Query("Transactions", NULL);
  EXPECT_EQ(simdur, qr.rows.size()) << "failed to order+run on fresh fuel inside 1 time step";
}

// tests that the correct number of assemblies are popped from the core each
// cycle.
TEST(BulkReactorTests, BatchSizes) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <core_mass>300</core_mass>  "
     "  <discharge_mass>10</discharge_mass>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", c_uox());
  sim.AddRecipe("spentuox", c_spentuox());
  int id = sim.Run();

  QueryResult qr = sim.db().Query("Transactions", NULL);
  // 7 for initial core, 3 per time step for each new batch for remainder
  EXPECT_EQ(1+1*(simdur-1), qr.rows.size());
}

// tests that the refueling period between cycle end and start of the next
// cycle is honored.
TEST(BulkReactorTests, RefuelTimes) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>4</cycle_time>  "
     "  <refuel_time>3</refuel_time>  "
     "  <core_mass>300</core_mass>  "
     "  <discharge_mass>10</discharge_mass>  ";
  int simdur = 49;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", c_uox());
  sim.AddRecipe("spentuox", c_spentuox());
  int id = sim.Run();

  QueryResult qr = sim.db().Query("Transactions", NULL);
  int cyclet = 4;
  int refuelt = 3;
  int n_assem_want = simdur/(cyclet+refuelt)+1; // +1 for initial core
  EXPECT_EQ(n_assem_want, qr.rows.size());
}


// tests that a reactor decommissions on time without producing
// power at the end of its lifetime.
TEST(BulkReactorTests, DecomTimes) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>2</cycle_time>  "
     "  <refuel_time>2</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <core_mass>300</core_mass>  "
     "  <discharge_mass>10</discharge_mass>  ";
  int simdur = 12;
  int lifetime = 7;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur, lifetime);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", c_uox());
  sim.AddRecipe("spentuox", c_spentuox());
  int id = sim.Run();

  // operating for 2+2 months and shutdown for 2+1
  int on_time = 4;
  std::vector<Cond> conds;
  conds.push_back(Cond("Value", "==", 1000));
  QueryResult qr = sim.db().Query("TimeSeriesPower", &conds);
  EXPECT_EQ(on_time, qr.rows.size());

  int off_time = 3;
  conds.clear();
  conds.push_back(Cond("Value", "==", 0));
  qr = sim.db().Query("TimeSeriesPower", &conds);
  EXPECT_EQ(off_time, qr.rows.size());
}


// Tests if a reactor produces power at the time of its decommission
// given a refuel_time of zero.
TEST(BulkReactorTests, DecomZeroRefuel) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>2</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <assem_size>1</assem_size>  "
     "  <core_mass>300</core_mass>  "
     "  <discharge_mass>10</discharge_mass>  ";

  int simdur = 8;
  int lifetime = 6;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur, lifetime);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", c_uox());
  sim.AddRecipe("spentuox", c_spentuox());
  int id = sim.Run();

  // operating for 2+2 months and shutdown for 2+1
  int on_time = 6;
  std::vector<Cond> conds;
  conds.push_back(Cond("Value", "==", 1000));
  QueryResult qr = sim.db().Query("TimeSeriesPower", &conds);
  EXPECT_EQ(on_time, qr.rows.size());
}

// tests that new fuel is ordered immediately following cycle end - at the
// start of the refueling period - not before and not after. - thie is subtly
// different than RefuelTimes test and is not a duplicate of it.
TEST(BulkReactorTests, OrderAtRefuelStart) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>4</cycle_time>  "
     "  <refuel_time>3</refuel_time>  "
     "  <core_mass>300</core_mass>  "
     "  <discharge_mass>10</discharge_mass>  ";
  int simdur = 7;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", c_uox());
  sim.AddRecipe("spentuox", c_spentuox());
  int id = sim.Run();

  QueryResult qr = sim.db().Query("Transactions", NULL);
  int cyclet = 4;
  int refuelt = 3;
  int n_assem_want = simdur/(cyclet+refuelt)+1; // +1 for initial core
  EXPECT_EQ(n_assem_want, qr.rows.size());
}

// tests that the reactor halts operation when it has no more room in its
// spent fuel inventory buffer.
TEST(BulkReactorTests, FullSpentInventory) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <core_mass>300</core_mass>  "
     "  <spent_mass>30</spent_mass> "
     "  <discharge_mass>10</discharge_mass>  ";

  int simdur = 10;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", c_uox());
  sim.AddRecipe("spentuox", c_spentuox());
  int id = sim.Run();

  QueryResult qr = sim.db().Query("Transactions", NULL);
  int n_assem_spent = 3;

  // +1 is for the assembly in the core + the three in spent
  EXPECT_EQ(n_assem_spent+1, qr.rows.size());
}


// tests that the reactor cycle is delayed as expected when it is unable to
// acquire fuel in time for the next cycle start.  This checks that after a
// cycle is delayed past an original scheduled start time, as soon as enough fuel is
// received, a new cycle pattern is established starting from the delayed
// start time.
TEST(BulkReactorTests, FuelShortage) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>7</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <core_mass>300</core_mass>  "
     "  <discharge_mass>10</discharge_mass>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("uox").lifetime(1).Finalize(); // provide initial full batch
  sim.AddSource("uox").start(9).lifetime(1).capacity(5).Finalize(); // provide partial batch post cycle-end
  sim.AddSource("uox").start(15).Finalize(); // provide remainder of batch much later
  sim.AddRecipe("uox", c_uox());
  sim.AddRecipe("spentuox", c_spentuox());
  int id = sim.Run();

  // check that we never got a full refueled batch during refuel period
  std::vector<Cond> conds;
  conds.push_back(Cond("Time", "<", 15));
  QueryResult qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(5, qr.rows.size());

  // after being delayed past original scheduled start of new cycle, we got
  // final assembly for core.
  conds.clear();
  conds.push_back(Cond("Time", "==", 15));
  qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(1, qr.rows.size());

  // all during the next (delayed) cycle we shouldn't be requesting any new fuel
  conds.clear();
  conds.push_back(Cond("Time", "<", 21));
  qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(6, qr.rows.size());

  // as soon as this delayed cycle ends, we should be requesting/getting 3 new batches
  conds.clear();
  conds.push_back(Cond("Time", "==", 22));
  qr = sim.db().Query("Transactions", &conds);
  EXPECT_EQ(3, qr.rows.size());
}

// tests that discharged fuel is transmuted properly immediately at cycle end.
TEST(BulkReactorTests, DischargedFuelTransmute) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>4</cycle_time>  "
     "  <refuel_time>3</refuel_time>  "
     "  <core_mass>300</core_mass>  "
     "  <discharge_mass>10</discharge_mass>  ";

  int simdur = 7;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddSink("waste").Finalize();
  sim.AddRecipe("uox", c_uox());
  Composition::Ptr spentuox = c_spentuox();
  sim.AddRecipe("spentuox", spentuox);
  int id = sim.Run();

  std::vector<Cond> conds;
  conds.push_back(Cond("SenderId", "==", id));
  int resid = sim.db().Query("Transactions", &conds).GetVal<int>("ResourceId");
  Material::Ptr m = sim.GetMaterial(resid);
  MatQuery mq(m);
  EXPECT_EQ(spentuox->id(), m->comp()->id());
  EXPECT_TRUE(mq.mass(942390000) > 0) << "transmuted spent fuel doesn't have Pu239";
}


TEST(BulkReactorTests, PositionInitialize) {
  std::string config =
     "  <fuel_inrecipes>  <val>lwr_fresh</val>  </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>lwr_spent</val>  </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>enriched_u</val> </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>      </fuel_outcommods>  "
     "  <fuel_prefs>      <val>1.0</val>        </fuel_prefs>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <core_mass>10</core_mass>  "
     "  <discharge_mass>1</discharge_mass>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("enriched_u").Finalize();
  sim.AddRecipe("lwr_fresh", c_uox());
  sim.AddRecipe("lwr_spent", c_spentuox());
  int id = sim.Run();

  QueryResult qr = sim.db().Query("AgentPosition", NULL);
  EXPECT_EQ(qr.GetVal<double>("Latitude"), 0.0);
  EXPECT_EQ(qr.GetVal<double>("Longitude"), 0.0);
}

TEST(BulkReactorTests, PositionInitialize2) {
  std::string config =
     "  <fuel_inrecipes>  <val>lwr_fresh</val>  </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>lwr_spent</val>  </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>enriched_u</val> </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>      </fuel_outcommods>  "
     "  <fuel_prefs>      <val>1.0</val>        </fuel_prefs>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>0</refuel_time>  "
     "  <core_mass>10</core_mass>  "
     "  <discharge_mass>1</discharge_mass>  ";

  int simdur = 50;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("enriched_u").Finalize();
  sim.AddRecipe("lwr_fresh", c_uox());
  sim.AddRecipe("lwr_spent", c_spentuox());
  int id = sim.Run();

  QueryResult qr = sim.db().Query("AgentPosition", NULL);
  EXPECT_EQ(qr.GetVal<double>("Latitude"), 30.0);
  EXPECT_EQ(qr.GetVal<double>("Longitude"), 30.0);
}

TEST(BulkReactorTests, ByProduct) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>1</refuel_time>  "
     "  <core_mass>10</core_mass>  "
     "  <discharge_mass>1</discharge_mass>  "
     ""
     "  <side_products> <val>process_heat</val> </side_products>"
     "  <side_product_quantity> <val>10</val> </side_product_quantity>";

  int simdur = 10;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", c_uox());
  sim.AddRecipe("spentuox", c_spentuox());
  int id = sim.Run();

  std::vector<Cond> conds;
  // test if it produces side products only when reactor is running
  int quantity = 10;
  conds.push_back(Cond("Value", "==", quantity));
  QueryResult qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(5, qr.rows.size());

  // test if it doesn't produce side products when reactor is refueling
  conds.clear();
  conds.push_back(Cond("Value", "==", 0));
  qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(5, qr.rows.size());
}

TEST(BulkReactorTests, MultipleByProduct) {
  std::string config =
     "  <fuel_inrecipes>  <val>uox</val>      </fuel_inrecipes>  "
     "  <fuel_outrecipes> <val>spentuox</val> </fuel_outrecipes>  "
     "  <fuel_incommods>  <val>uox</val>      </fuel_incommods>  "
     "  <fuel_outcommods> <val>waste</val>    </fuel_outcommods>  "
     ""
     "  <cycle_time>1</cycle_time>  "
     "  <refuel_time>1</refuel_time>  "
     "  <core_mass>10</core_mass>  "
     "  <discharge_mass>1</discharge_mass>  "
     ""
     "  <side_products> <val>process_heat</val> <val>water</val> </side_products>"
     "  <side_product_quantity> <val>10</val> <val>100</val> </side_product_quantity>";

  int simdur = 10;
  cyclus::MockSim sim(cyclus::AgentSpec(":cycamore:BulkReactor"), config, simdur);
  sim.AddSource("uox").Finalize();
  sim.AddRecipe("uox", c_uox());
  sim.AddRecipe("spentuox", c_spentuox());
  int id = sim.Run();


  std::vector<Cond> conds;
  // test if it produces heat when reactor is running
  int quantity = 10;
  conds.push_back(Cond("Product", "==", std::string("process_heat")));
  conds.push_back(Cond("Value", "==", quantity));
  QueryResult qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(5, qr.rows.size());

  // test if it produces water when reactor is running
  conds.clear();
  quantity = 100;
  conds.push_back(Cond("Product", "==", std::string("water")));
  conds.push_back(Cond("Value", "==", quantity));
  qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(5, qr.rows.size());

  conds.clear();
  conds.push_back(Cond("Value", "==", 0));
  qr = sim.db().Query("ReactorSideProducts", &conds);
  EXPECT_EQ(10, qr.rows.size());

}

} // namespace bulk_reactortests
} // namespace cycamore

