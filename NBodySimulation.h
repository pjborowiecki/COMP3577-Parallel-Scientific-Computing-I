#include <cmath>

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

class NBodySimulation {
public:
// protected:
  double t;
  double tFinal;
  double tPlot;
  double tPlotDelta;

  int NumberOfBodies;

  /**
   * Storing this way is better for cache
  */
  double* xx __attribute__((aligned(64)));
  double* xy __attribute__((aligned(64)));
  double* xz __attribute__((aligned(64)));

  double* vx __attribute__((aligned(64)));
  double* vy __attribute__((aligned(64)));
  double* vz __attribute__((aligned(64)));

  double* ax __attribute__((aligned(64)));
  double* ay __attribute__((aligned(64)));
  double* az __attribute__((aligned(64)));
  
  double* m  __attribute__((aligned(64)));

  // C = 10^(-2)/NumberOfBodies
  double C;

  /**
   * Global time step size used.
   */
  double timeStepSize;

  /**
   * Maximum velocity of all particles.
   */
  double maxV;

  /**
   * Minimum distance between two elements.
   */
  double minDx;

  /**
   * Stream for video output file.
   */
  std::ofstream videoFile;

  /**
   * Output counters.
   */
  int snapshotCounter;
  int timeStepCounter;


// public:
  NBodySimulation ();
  ~NBodySimulation ();

  /**
   * Check that the number command line parameters is correct.
   */
  void checkInput(int argc, char** argv);

  /**
   * Set up scenario from the command line.
   *
   * If you need additional helper data structures, you can initialise them
   * here. Alternatively, you can introduce a totally new function to initialise
   * additional data fields and call this new function from main after setUp().
   * Either way is fine.
   *
   * The semantics of this operations are not to be changed in the assignment.
   */
  void setUp (int argc, char** argv);

  bool process_gravity_and_detect_collision();
  void process_collisions();
  void handle_collision(int i, int j);
  
  /**
   * Implement timestepping scheme and force updates.
   */
  void updateBody ();

  /**
   * Check if the last time step has been reached (simulation is completed).
   *
   * This operation is not to be changed in the assignment.
   */
  bool hasReachedEnd ();

  /**
   * Take simulations snapshopts and print summary to standard output.
   *
   * This operation is not to be changed in the assignment.
   */
  void takeSnapshot ();

  /**
   * Handle Paraview output.
   *
   * These operations are not to be changed in the assignment.
   *
   * The file format is documented at
   * http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
   */
  void openParaviewVideoFile ();
  void closeParaviewVideoFile ();
  void printParaviewSnapshot ();

  /**
   * Handle terminal output.
   *
   * These operations are not to be changed in the assignment.
   */
  void printSnapshotSummary ();
  void printSummary ();

};
