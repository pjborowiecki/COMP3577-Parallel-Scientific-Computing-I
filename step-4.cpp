#include <iomanip>

#include "NBodySimulationVectorised.cpp"

/**
 * You can compile this file with
 *   make step-4-gcc   // Uses the GNU Compiler Collection.
 *   make asisigment-icpc  // Uses the Intel compiler.
 * and run it with
 *   ./step-4-gcc
 *   ./step-4-icpc
 *
 * Results will be added to the `paraview-output` directory. In it you will find
 * a result.pvd file that you can open with ParaView. To see the points you will
 * need to look a the properties of result.pvd and select the representation
 * "Point Gaussian". Pressing play will play your time steps.
 */

class NBodySimulationParallelised : public NBodySimulationVectorised {

  /**
   * Due to data race, symmetry is not exploited.
   * This could be mitigated by, for example, having each thread work on its own
   * copy of the acceleration data, and then summing partial contributions from 
   * each thread to the global array.
   * But this would require O(N*<number of threads>) extra memory, which could 
   * become problematic for very large scale simulations
  */
  bool process_gravity_and_detect_collision()
  {
    std::fill(ax, ax+NumberOfBodies, 0);
    std::fill(ay, ay+NumberOfBodies, 0);
    std::fill(az, az+NumberOfBodies, 0);
    double m_minDx = std::numeric_limits<double>::max();
    double m_minC = std::numeric_limits<double>::max();

    #pragma omp parallel for reduction(min:m_minDx,m_minC)
    for (int i = 0; i < NumberOfBodies; ++i){
      double axi(0),ayi(0),azi(0);
      double xxi(xx[i]), xyi(xy[i]), xzi(xz[i]), mi(m[i]);

      double t_minDx = std::numeric_limits<double>::max();
      double t_minC = std::numeric_limits<double>::max();

      #pragma omp simd reduction(+:axi,ayi,azi) reduction(min:t_minDx,t_minC)
      for (int j=i-1; j>=0; --j){
        double dx = xx[j]-xxi;
        double dy = xy[j]-xyi;
        double dz = xz[j]-xzi;
        double dst2 = dx*dx + dy*dy + dz*dz;
        double dst = std::sqrt(dst2);
        double dst3 = dst2 * dst;
        
        double gx = dx/dst3;
        double gy = dy/dst3;
        double gz = dz/dst3;

        axi += gx*m[j];
        ayi += gy*m[j];
        azi += gz*m[j];
        t_minC  = std::min(t_minC, dst/(mi + m[j]));
        t_minDx = std::min(t_minDx, dst);
      }

      #pragma omp simd reduction(+:axi,ayi,azi) reduction(min:t_minDx,t_minC)
      for (int j=i+1; j<NumberOfBodies; ++j){
        double dx = xx[j]-xxi;
        double dy = xy[j]-xyi;
        double dz = xz[j]-xzi;
        double dst2 = dx*dx + dy*dy + dz*dz;
        double dst = std::sqrt(dst2);
        double dst3 = dst2 * dst;
        
        double gx = dx/dst3;
        double gy = dy/dst3;
        double gz = dz/dst3;

        axi += gx*m[j];
        ayi += gy*m[j];
        azi += gz*m[j];
        t_minC  = std::min(t_minC, dst/(mi + m[j]));
        t_minDx = std::min(t_minDx, dst);
      }

      ax[i] += axi;
      ay[i] += ayi;
      az[i] += azi;
      m_minC  = std::min(m_minC, t_minC);
      m_minDx = std::min(m_minDx, t_minDx);
    }
    
    minDx = m_minDx;
    return m_minC <= C;
  }
  
public:
  void updateBody () {
    timeStepCounter++;
    maxV   = 0.0;
    minDx  = std::numeric_limits<double>::max();

    #pragma omp parallel for simd
    for (int i = 0; i<NumberOfBodies; ++i){    
      vx[i] += timeStepSize/2 * ax[i];
      vy[i] += timeStepSize/2 * ay[i];
      vz[i] += timeStepSize/2 * az[i];

      xx[i] += timeStepSize   * vx[i];
      xy[i] += timeStepSize   * vy[i];
      xz[i] += timeStepSize   * vz[i];
    }

    if (process_gravity_and_detect_collision())
    {
      process_collisions();
      process_gravity_and_detect_collision();
    }

    double m_maxV = 0;
    #pragma omp parallel for simd reduction(max:m_maxV)
    for (int i = 0; i<NumberOfBodies; ++i){    
      vx[i] += timeStepSize/2 * ax[i];
      vy[i] += timeStepSize/2 * ay[i];
      vz[i] += timeStepSize/2 * az[i];
      
      m_maxV = std::max(maxV, std::sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]));
    }
    
    maxV = m_maxV;
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
  NBodySimulationParallelised nbs;
  nbs.setUp(argc,argv);
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
