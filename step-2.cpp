#define CUTOFF_RADIUS 1e-1      // cutoff radius for short-distance force


#include <iomanip>

#include "NBodySimulation.h"

/**
 * You can compile this file with
 *   make step-2-gcc   // Uses the GNU Compiler Collection.
 *   make asisigment-icpc  // Uses the Intel compiler.
 * and run it with
 *   ./step-2-gcc
 *   ./step-2-icpc
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */


/**
 * Lookup nearby particles in time O(1) (on average) * by keeping track of the 
 * cell membership of each particle.
 * 
 * - Grid (list of cells) is held in std::unordered_map that has average lookup
 *   time O(1)
 * - Cells are accessed by key of type std::tuple<int,int,int> representing
 *   their location in 3d
 * - Cells themselves are represented by std::unordered_set<int> for storing ids
 *   of their memebers
*/
#include <unordered_set>
#include <unordered_map>
#include <tuple>

struct Grid{
  typedef std::tuple<int,int,int> CellID;
  double cell_size;

  Grid() : cell_size(1.0) {};
  Grid(double cell_size) : cell_size(cell_size) {};

  struct hash{
    size_t operator()(const CellID& x) const {
      return std::get<0>(x) ^ std::get<1>(x) ^ std::get<2>(x);
    }
  };
  
  std::unordered_map<CellID, std::unordered_set<int>, hash> cells;

  CellID coordsToCellID(double x, double y, double z){
    CellID r = {
      std::floor(x/cell_size),
      std::floor(y/cell_size),
      std::floor(z/cell_size)
    };
    return  r;
  }
};
bool operator==(const Grid::CellID& lhs, const Grid::CellID& rhs) {
  return (
     std::get<0>(lhs) == std::get<0>(rhs) && 
     std::get<1>(lhs) == std::get<1>(rhs) &&
     std::get<2>(lhs) == std::get<2>(rhs)
    );
}

/**
 * O(N) simulation of molecular forces with cutoff radius
*/
class NBodySimulationMolecularForces : public NBodySimulation {
private:
  Grid grid;

public:
  void setUpGrid(double cell_size){
    grid = Grid(cell_size);
    for (int i = 0; i < NumberOfBodies; ++i){
      Grid::CellID cid = grid.coordsToCellID(xx[i],xy[i],xz[i]);
      grid.cells.emplace(cid, std::unordered_set<int>());
      grid.cells[cid].insert(i);
    }
  }

  void calc_force(int i, int j, double &fx, double &fy, double &fz, double &dst){
    double dx = xx[i] - xx[j];
    double dy = xy[i] - xy[j];
    double dz = xz[i] - xz[j];
    double dst2 = dx*dx + dy*dy + dz*dz;
    dst = std::sqrt(dst2);

    if (dst > CUTOFF_RADIUS){
      fx = 0; fy = 0; fz = 0;
    }
    else{
      double a = 0.1/dst;   // 0.1/r
      double b = a*a; b*=b; // b := (0.1/r)^4
      double c = b*b*a;     // c := (0.1/r)^9
      double fmag = 10 * c * (b - 1);      
      fx = fmag*dx;
      fy = fmag*dy;
      fz = fmag*dz;
    }
  }

  void process_interactions(){

    // Clear acceleration data
    std::fill(ax, ax+NumberOfBodies, 0);
    std::fill(ay, ay+NumberOfBodies, 0);
    std::fill(az, az+NumberOfBodies, 0);

    /**
     * For each particle i:
     *   let c := cell containing i
     *   For each particle j in c or neighbour of c:
     *     calculate interaction between i and j
     *     update acceleration of i
     * 
     * Note: no symmetry exploited here, but for uniform distribution of 
     * particles this is O(N)
    */
    for (int i = 0; i < NumberOfBodies; ++i){
      int cx,cy,cz;
      std::tie(cx,cy,cz) = grid.coordsToCellID(xx[i], xy[i], xz[i]);

      Grid::CellID nid;

      // TODO: Check which octant of the cell is occupied by the particle, and 
      //       only check the neighbouring cells adjacent to this octant

      // iterate over all neighbours and self (27 iterations)
      for (std::get<0>(nid)=cx-1;std::get<0>(nid)<=cx+1;++std::get<0>(nid)){
        for (std::get<1>(nid)=cy-1;std::get<1>(nid)<=cy+1;++std::get<1>(nid)){
          for (std::get<2>(nid)=cz-1; std::get<2>(nid)<=cz+1;++std::get<2>(nid))
          {

            // If the neighbour cell does not exist in memory - go to next one
            if (grid.cells.find(nid) == grid.cells.end()) continue;

            for (int j : grid.cells[nid]){
              if (i == j) continue;
              double fx,fy,fz,dst;
              calc_force(i,j,fx,fy,fz,dst);
              minDx = std::min(minDx,dst);
              ay[i] += fx/m[i];
              ax[i] += fy/m[i];
              az[i] += fz/m[i];
            }
          } 
        }  
      }
    }
  }

  void updateBody(){
    timeStepCounter++;
    maxV   = 0.0;
    minDx  = std::numeric_limits<double>::max();

    for (int i = 0; i<NumberOfBodies; ++i){    
      // 1. Compute half an Euler time step for v
      // v(t + dt/2) = v(t) + dt/2 * a(t)
      vx[i] += timeStepSize/2 * ax[i];
      vy[i] += timeStepSize/2 * ay[i];
      vz[i] += timeStepSize/2 * az[i];
    
      // 2. Update positions (and cell membership)
      // x(t+dt) = d(t) + dt * v(t + dt/2)
      Grid::CellID cell_id = grid.coordsToCellID(xx[i], xy[i], xz[i]);
      xx[i] += timeStepSize * vx[i];
      xy[i] += timeStepSize * vy[i];
      xz[i] += timeStepSize * vz[i];

      Grid::CellID new_cell_id = grid.coordsToCellID(xx[i], xy[i], xz[i]);
      if (!(new_cell_id == cell_id)){
        grid.cells[cell_id].erase(i);
        grid.cells[new_cell_id].insert(i);

        // free up memory if the old cell is left empty
        if (grid.cells[cell_id].size() == 0) grid.cells.erase(cell_id);
      }
    }

    // 3. Calculate acceleration
    process_interactions();

    // 4. Update the velocities
    // v(t + dt) = v(t + dt/2) + dt/2 * a(t + dt)
    // Euler time step
    for (int i = 0; i<NumberOfBodies; ++i){    
      vx[i] += timeStepSize/2 * ax[i];
      vy[i] += timeStepSize/2 * ay[i];
      vz[i] += timeStepSize/2 * az[i];
      
      maxV = std::max(maxV, std::sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]));
    }

    t += timeStepSize;
  }
};


/**
 * Main routine.
 *
 * No major changes are needed in the assignment. You can add initialisation or
 * or remove input checking, if you feel the need to do so. But keep in mind
 * that you may not alter what the program writes to the standard output.
 */

int main (int argc, char** argv) {

  std::cout << std::setprecision(15);

  // Code that initialises and runs the simulation.
  NBodySimulationMolecularForces nbs;
  nbs.setUp(argc,argv);
  nbs.setUpGrid(CUTOFF_RADIUS);
  nbs.openParaviewVideoFile();
  nbs.takeSnapshot();

  while (!nbs.hasReachedEnd()) {
    nbs.updateBody();
    nbs.takeSnapshot();
  }

  nbs.printSummary();
  nbs.closeParaviewVideoFile();

  return 0;
}
