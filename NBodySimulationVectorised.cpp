#include "NBodySimulation.h"

class NBodySimulationVectorised : public NBodySimulation {
protected:

  /**
   * Both g++ 9.4.0 and icpc 2021.8.0 vectorise the inner loop
   * Not much savings to be made elsewhere
  */
  bool process_gravity_and_detect_collision()
  {
    std::fill(ax, ax+NumberOfBodies, 0);
    std::fill(ay, ay+NumberOfBodies, 0);
    std::fill(az, az+NumberOfBodies, 0);
    double m_minDx = std::numeric_limits<double>::max();
    double m_minC = std::numeric_limits<double>::max();

    for (int i = 0; i<NumberOfBodies; ++i){
      double axi(0),ayi(0),azi(0);
      double xxi(xx[i]), xyi(xy[i]), xzi(xz[i]), mi(m[i]);
     #pragma omp simd \
        reduction(+:axi,ayi,azi) \
        reduction(min:m_minDx,m_minC)
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
        ax[j] -= gx*mi;
        ay[j] -= gy*mi;
        az[j] -= gz*mi;

        m_minDx = std::min(m_minDx, dst);
        m_minC = std::min(m_minC, dst/(mi + m[j]));
      }

      if (m_minC < C) return true;
      ax[i] += axi;
      ay[i] += ayi;
      az[i] += azi;
    }

    minDx = m_minDx;
    return false;
  }


public:
  void updateBody () {
    timeStepCounter++;
    maxV   = 0.0;
    minDx  = std::numeric_limits<double>::max();

    #pragma omp simd
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
    #pragma omp simd reduction(max:m_maxV)
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
