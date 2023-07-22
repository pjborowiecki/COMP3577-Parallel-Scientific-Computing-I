#include "NBodySimulation.h"

NBodySimulation::NBodySimulation () :
  t(0), tFinal(0), tPlot(0), tPlotDelta(0), NumberOfBodies(0),
  xx(nullptr), xy(nullptr), xz(nullptr),
  vx(nullptr), vy(nullptr), vz(nullptr),
  ax(nullptr), ay(nullptr), az(nullptr),
  timeStepSize(0), maxV(0), minDx(0), videoFile(nullptr),
  snapshotCounter(0), timeStepCounter(0) {};

NBodySimulation::~NBodySimulation () {
  if (xx != nullptr) free(xx);
  if (xy != nullptr) free(xy);
  if (xz != nullptr) free(xz);
  if (vx != nullptr) free(vx);
  if (vy != nullptr) free(vy);
  if (vz != nullptr) free(vz);
  if (ax != nullptr) free(ax);
  if (ay != nullptr) free(ay);
  if (az != nullptr) free(az);
  if (m  != nullptr) free(m);
}

void NBodySimulation::checkInput(int argc, char** argv) {
    if (argc==1) {
    std::cerr << "usage: " << std::string(argv[0])
              << " plot-time final-time dt objects" << std::endl
              << " Details:" << std::endl
              << " ----------------------------------" << std::endl
              << "  plot-time:       interval after how many time units to plot."
                 " Use 0 to switch off plotting" << std::endl
              << "  final-time:      simulated time (greater 0)" << std::endl
              << "  dt:              time step size (greater 0)" << std::endl
              << "  objects:         any number of bodies, specified by position, velocity, mass" << std::endl
              << std::endl
              << "Examples of arguments:" << std::endl
              << "+ One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ One body spiralling around the other" << std::endl
              << "    0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0" << std::endl
              << "+ Three-body setup from first lecture" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0" << std::endl
              << "+ Five-body setup" << std::endl
              << "    0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0" << std::endl
              << std::endl;

    throw -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each body is given by seven entries"
                 " (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments"
                 " (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    throw -2;
  }
}

void NBodySimulation::setUp (int argc, char** argv) {

  checkInput(argc, argv);

  NumberOfBodies = (argc-4) / 7;
  C = 1e-2/NumberOfBodies;

  xx = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));
  xy = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));
  xz = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));
  vx = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));
  vy = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));
  vz = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));
  ax = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));
  ay = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));
  az = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));
  m  = static_cast<double*>(aligned_alloc(64, NumberOfBodies * sizeof(double)));

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  //maxV = 0.0;
  for (int i=0; i<NumberOfBodies; i++) {
    xx[i] = std::stof(argv[readArgument]); readArgument++;
    xy[i] = std::stof(argv[readArgument]); readArgument++;
    xz[i] = std::stof(argv[readArgument]); readArgument++;

    vx[i] = std::stof(argv[readArgument]); readArgument++;
    vy[i] = std::stof(argv[readArgument]); readArgument++;
    vz[i] = std::stof(argv[readArgument]); readArgument++;
    //maxV = std::max(maxV, std::sqrt(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]));

    m[i] = std::stof(argv[readArgument]); readArgument++;
    if (m[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies"
            << std::endl;

  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta
              << " time units" << std::endl;
    tPlot = 0.0;
  }
}

void NBodySimulation::handle_collision(int i, int j)
{
  // Update position, velocity and mass of body i
  double Minv = 1.0/(m[i]+m[j]);
  xx[i] = (m[i]*xx[i] + m[j]*xx[j])*Minv;
  xy[i] = (m[i]*xy[i] + m[j]*xy[j])*Minv;
  xz[i] = (m[i]*xz[i] + m[j]*xz[j])*Minv;
  
  vx[i] = (m[i]*vx[i] + m[j]*vx[j])*Minv;
  vy[i] = (m[i]*vy[i] + m[j]*vy[j])*Minv;
  vz[i] = (m[i]*vz[i] + m[j]*vz[j])*Minv;

  m[i] += m[j];


  // Remove body j by replacing it with the last one in the array
  xx[j] = xx[NumberOfBodies-1];
  xy[j] = xy[NumberOfBodies-1];
  vx[j] = xz[NumberOfBodies-1];
  
  xz[j] = vx[NumberOfBodies-1];
  vx[j] = vy[NumberOfBodies-1];
  vz[j] = vz[NumberOfBodies-1];

  m[j] = m[NumberOfBodies-1]; 
  

  // "Remove" the last body
  NumberOfBodies--;
}

/**
 * On a rare occasion of two bodies colliding, we need to find out which ones
 * and deal with them. We could remember the collision data while computing 
 * force interactions, but this makes vectorisation and parallelisation 
 * harder, while potential time savings are not worth the messy code.
*/
void NBodySimulation::process_collisions()
{
  for (int i = 0; i < NumberOfBodies; ++i){
    for (int j = i+1; j < NumberOfBodies; ++j){
      double dx = xx[j]-xx[i];
      double dy = xy[j]-xy[i];
      double dz = xz[j]-xz[i];
      double dst = sqrt(dx*dx + dy*dy + dz*dz);
      
      if (dst / (m[i] + m[j]) <= C){
        handle_collision(i,j);
      }
    }
  }
}

/**
 * Calculating distance is the most expensive, so we can do it only once 
 * and use it for both gravity and collision detection 
*/
bool NBodySimulation::process_gravity_and_detect_collision()
{
  // Clear acceleration data
  std::fill(ax, ax+NumberOfBodies, 0);
  std::fill(ay, ay+NumberOfBodies, 0);
  std::fill(az, az+NumberOfBodies, 0);

  for (int i = 0; i<NumberOfBodies; ++i){
    double axi(0),ayi(0),azi(0);
    double xxi(xx[i]), xyi(xy[i]), xzi(xz[i]), mi(m[i]);

    for (int j=i+1; j<NumberOfBodies; ++j){
      double dx = xx[j]-xxi;
      double dy = xy[j]-xyi;
      double dz = xz[j]-xzi;
      double dst2 = dx*dx + dy*dy + dz*dz;
      double dst = std::sqrt(dst2);

      // If there is a collision, then we can stop right there
      // - forces can be calculated after all collisions are handled
      if (dst/(mi + m[j]) <= C) return true; 

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

      minDx = std::min(minDx, dst);
    }

    ax[i] += axi;
    ay[i] += ayi;
    az[i] += azi;
  }

  return false;
}

void NBodySimulation::updateBody () {

  timeStepCounter++;
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  for (int i = 0; i<NumberOfBodies; ++i){    
    // 1. Compute half an Euler time step for v
    // v(t + dt/2) = v(t) + dt/2 * a(t)
    vx[i] += timeStepSize/2 * ax[i];
    vy[i] += timeStepSize/2 * ay[i];
    vz[i] += timeStepSize/2 * az[i];
  
    // 2. Update positions
    // x(t+dt) = d(t) + dt * v(t + dt/2)
    xx[i] += timeStepSize * vx[i];
    xy[i] += timeStepSize * vy[i];
    xz[i] += timeStepSize * vz[i];
  }

  // 3. Calculate acceleration
  if (process_gravity_and_detect_collision())
  {
    // if there are collisions - process them and recalculate acceleration
    process_collisions();
    process_gravity_and_detect_collision();
  }

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



/**
 * Check if simulation has been completed.
 */
bool NBodySimulation::hasReachedEnd () {
  return t > tFinal;
}

void NBodySimulation::takeSnapshot () {
  if (t >= tPlot) {
    printParaviewSnapshot();
    printSnapshotSummary();
    tPlot += tPlotDelta;
  }
}


void NBodySimulation::openParaviewVideoFile () {
  videoFile.open("paraview-output/result.pvd");
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\""
    " version=\"0.1\""
    " byte_order=\"LittleEndian\""
    " compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}

void NBodySimulation::closeParaviewVideoFile () {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
  videoFile.close();
}

void NBodySimulation::printParaviewSnapshot () {
  static int counter = -1;
  counter++;
  std::stringstream filename, filename_nofolder;
  filename << "paraview-output/result-" << counter <<  ".vtp";
  filename_nofolder << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\""
    " NumberOfComponents=\"3\""
    " format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << xx[i] << " " << xy[i] << " " << xz[i] << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  out.close();

  videoFile << "<DataSet timestep=\"" << counter
            << "\" group=\"\" part=\"0\" file=\"" << filename_nofolder.str()
            << "\"/>" << std::endl;
}

void NBodySimulation::printSnapshotSummary () {
  std::cout << "plot next snapshot"
            << ",\t time step=" << timeStepCounter
            << ",\t t="         << t
            << ",\t dt="        << timeStepSize
            << ",\t v_max="     << maxV
            << ",\t dx_min="    << minDx
            << std::endl;
}

void NBodySimulation::printSummary () {
  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: "
            << xx[0] << ", " << xy[0] << ", " << xz[0] << std::endl;
}
