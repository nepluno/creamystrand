/**
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "StrandSimulatorApp.hh"

#include <boost/filesystem.hpp>
#include <chrono>
#include <thread>

#include "ProblemStepper.hh"
#include "Render/Camera.hh"
#include "Render/ViewController.hh"

// problems :
#include "Problems/XMLReader.hh"
#include "StrandSim/Utils/TextLog.hh"

std::vector<ProblemStepper*> problems;
int g_problem_idx = -1;

int g_window_width = 960;
int g_window_height = 540;

bool g_single_step = false;
bool g_paused = true;
bool g_exit = false;
bool g_save = false;

time_t g_rawtime;
struct tm*
    g_timeinfo;  // for time-stamped simulation capture & output directory

std::string g_outputdirectory;  // directory to save to
int g_dumpcoord = 0;            // dump coords
int g_dumpbinary = 0;
int g_render = 1;

int g_last_frame_num = -1;       // last frame # that was output
int g_last_checkpoint_num = -1;  // last frame # that was output
int g_current_frame = 0;
int g_current_checkpoint = 0;

ProblemStepper* g_ps;
ViewController controller;
void doSimulation();
void output();

void createProblemVector() {
  problems.push_back(new XMLReader());  // default problem loader
}

void setOptions() {}

void getOptions() {}

void display() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();

  controller.ApplyCamera();

  g_ps->render();

  glutSwapBuffers();
  glPopMatrix();
}

void reshape(int w, int h) {
  g_window_width = w;
  g_window_height = h;

  Camera& c = controller.getCamera();
  c.setPerspective(60, 1.0);
  const Scalar radius = controller.getBoundingRadius();
  c.setZClipping(0.01 * radius, 3 * radius);
  c.setViewport(w, h);

  glutPostRedisplay();
}

void output() {
  if (g_dumpcoord) {
    int steps_per_frame = g_dumpcoord;
    const int step = floor(g_ps->getTime() / g_ps->getDt() + 0.5);

    const int last_frame_num = step / steps_per_frame;
    std::cout << "step: " << step << " steps / frame: " << steps_per_frame
              << " last step: " << last_frame_num << std::endl;
    if (g_last_frame_num != last_frame_num) {
      g_last_frame_num = last_frame_num;

      if (!boost::filesystem::exists(g_outputdirectory)) {
        boost::filesystem::create_directory(g_outputdirectory);
        boost::filesystem::permissions(
            g_outputdirectory,
            boost::filesystem::add_perms | boost::filesystem::others_all |
                boost::filesystem::owner_all | boost::filesystem::group_all);
      }

      const int file_width = 20;

      g_ps->dumpData(g_outputdirectory, g_current_frame, file_width);
      ++g_current_frame;
    }
  }

  if (g_dumpbinary) {
    int step_per_checkpoint = g_dumpbinary;
    const int step = floor(g_ps->getTime() / g_ps->getDt() + 0.5);

    const int last_checkpoint_num = step / step_per_checkpoint;
    std::cout << "step: " << step
              << " steps / checkpoint: " << step_per_checkpoint
              << " last checkpoint: " << last_checkpoint_num << std::endl;
    if (g_last_checkpoint_num != last_checkpoint_num) {
      g_last_checkpoint_num = last_checkpoint_num;

      if (!boost::filesystem::exists(g_outputdirectory)) {
        boost::filesystem::create_directory(g_outputdirectory);
        boost::filesystem::permissions(
            g_outputdirectory,
            boost::filesystem::add_perms | boost::filesystem::others_all |
                boost::filesystem::owner_all | boost::filesystem::group_all);
      }

      const int file_width = 20;

      ++g_current_checkpoint;

      g_ps->dumpBinaryCheckpoint(g_outputdirectory, g_current_frame,
                                 g_current_checkpoint, file_width);
    }
  }
}

void idle() {
  // If in render mode, update the display
  if (g_render) {
#ifdef __APPLE__
    static int mojave_counter = 0;
    if (mojave_counter < 2)
      glutReshapeWindow(g_window_width - 1 + mojave_counter++, g_window_height);
#endif
    glutPostRedisplay();
  }
}

void setLighting() {
  // Create a directional white light with a small ambient component
  glEnable(GL_LIGHT0);
  GLfloat white_ambient[] = {0.1, 0.1, 0.1, 1.0};
  glLightfv(GL_LIGHT0, GL_AMBIENT, white_ambient);
  GLfloat white_diffuse[] = {0.55, 0.55, 0.55, 1.0};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, white_diffuse);
  GLfloat upper_corner[] = {1.0, 1.0, 1.0, 0.0};
  glLightfv(GL_LIGHT0, GL_POSITION, upper_corner);

  // Create a much weaker direction light
  glEnable(GL_LIGHT1);
  GLfloat weak_white_diffuse[] = {0.3, 0.3, 0.3, 1.0};
  glLightfv(GL_LIGHT1, GL_DIFFUSE, weak_white_diffuse);
  GLfloat negative_z[] = {0.0, 0.0, 1.0, 0.0};
  glLightfv(GL_LIGHT1, GL_POSITION, negative_z);

  glShadeModel(GL_FLAT);
}

Vec3d calcSimCenter() {
  Vec3d center = Vec3d::Zero();

  // TODO:
  //    for( std::vector<RenderBase*>::size_type i = 0; i <
  //    renderable_objects.size(); ++i )
  //    {
  //        center += renderable_objects[i]->calculateObjectCenter();
  //    }
  //
  //    if( !renderable_objects.empty() ) center /= ( (double)
  //    renderable_objects.size() );

  g_ps->getCenter(center);

  return center;
}

Scalar calcViewRadius(const Vec3d& simCenter) {
  Scalar radius = 0.0;

  // TODO:
  //    for( std::vector<RenderBase*>::size_type i = 0; i <
  //    renderable_objects.size(); ++i )
  //    {
  //        const Vec3d center = renderable_objects[i]->calculateObjectCenter();
  //        const Scalar r =
  //        renderable_objects[i]->calculateObjectBoundingRadius( center );
  //        radius = std::max( radius, r + ( center - simCenter ).norm() );
  //    }

  g_ps->getRadius(radius, simCenter);

  if (radius == 0.0) radius = 1.0;

  return radius;
}

void centerObject() {
  const Vec3d simCenter = calcSimCenter();
  controller.setCenterMode(ViewController::CENTER_OBJECT);
  controller.setViewCenter(simCenter);

  const Scalar radius = calcViewRadius(simCenter);
  controller.setBoundingRadius(radius);
}

void initCamera() {
  controller.setViewDirection(Vec3d(0, 0, -2));
  centerObject();
}

void scaleMousePos(int x, int y, Scalar& xx, Scalar& yy) {
  int w, h;
  controller.getCamera().getViewport(&w, &h);

  xx = 2 * x / (Scalar)(w - 1) - 1.0;
  yy = 2 * (h - y - 1) / (Scalar)(h - 1) - 1.0;
}

void mouse(int button, int state, int x, int y) {
  const bool zooming =
      (button == GLUT_RIGHT_BUTTON) ||
      ((button == GLUT_LEFT_BUTTON) && (glutGetModifiers() & GLUT_ACTIVE_CTRL));
  const bool translating =
      (button == GLUT_LEFT_BUTTON) && (glutGetModifiers() & GLUT_ACTIVE_SHIFT);
  const bool rotating =
      (button == GLUT_LEFT_BUTTON) && (glutGetModifiers() == 0);

  Scalar xx, yy;
  scaleMousePos(x, y, xx, yy);
  if (state == GLUT_DOWN) {
    if (translating) {
      controller.beginTranslationDrag(xx, yy);
    }

    if (zooming) {
      controller.beginZoomDrag(xx, yy);
    }

    if (rotating) {
      controller.beginRotationDrag(xx, yy);
    }
  } else {
    controller.endTranslationDrag(xx, yy);
    controller.endRotationDrag(xx, yy);
    controller.endZoomDrag(xx, yy);
  }

  glutPostRedisplay();
}

void motion(int x, int y) {
  Scalar xx, yy;
  scaleMousePos(x, y, xx, yy);
  controller.updateDrag(xx, yy);
  glutPostRedisplay();
}

void menu(int id) {
  switch (id) {
    case 'q':
      g_paused = true;
      g_exit = true;
    std:
      exit(EXIT_SUCCESS);
      break;

    case ' ':
      g_paused = !g_paused;
      break;
    case 'w':
      break;
    case 's':
      g_single_step = true;
      g_paused = !g_paused;
      break;
    case 'a':
      break;
    case 'd':
      break;
    case 'f':
      if (g_paused) {
        g_ps->step();
        // if in render mode, update the display
        if (g_render) {
          glutPostRedisplay();
        }
        output();
      }
      break;
  }
}

void keyboard(unsigned char key, int, int) {
  std::cout << "KEY PRESSED: " << key << std::endl;
  menu(key);
}

void initializeOpenGL(int argc, char** argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_ALPHA);
  glutInitWindowPosition(500, 500);
  glutInitWindowSize(g_window_width, g_window_height);
  glutCreateWindow(argv[0]);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  // glClearColor(161.0/255.0,161.0/255.0,161.0/255.0,0.0);
  // glClearColor(0.0/255.0,0.0/255.0,0.0/255.0,0.0);
  glClearColor(1.0, 1.0, 1.0, 0.0);

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);

  //    initMenu();
  initCamera();

  setLighting();
  // SetMaterial();
}

void cleanup() {
  for (std::vector<ProblemStepper*>::size_type i = 1; i < problems.size();
       ++i) {
    assert(problems[i] != NULL);
    delete problems[i];
  }

  if (g_log != NULL) {
    delete g_log;
    g_log = NULL;
  }
}

int parseCommandLine(int argc, char** argv) {
  ProblemConstraint<int> allowedProblems(0, (int)(problems.size() - 1));
  ProblemConstraint<int> statLoggingLevels(0, 2);

  try {
    TCLAP::CmdLine cmd(
        "A Multi-Scale Model for Coupling Strands with Shear-Dependent Liquid");

    // args to modify behavior of r
    TCLAP::ValueArg<std::string> file("f", "file", "XML file for a problem",
                                      true, "", "string", cmd);

    TCLAP::ValueArg<int> output("o", "outputfile",
                                "Between # steps several PLY files are written "
                                "to save simulation state to",
                                false, 0, "integer", cmd);

    TCLAP::ValueArg<int> checkpoint("c", "checkpoint",
                                    "Between # steps several binary files are "
                                    "written to cache simulation state to",
                                    false, 0, "integer", cmd);

    TCLAP::ValueArg<int> display(
        "d", "display",
        "Run the simulation with display enabled if 1, without if 0", false, 1,
        "integer", cmd);

    // args to start logging
    TCLAP::ValueArg<int> statlog(
        "l", "statlog",
        "Log runtime stats: 0: no stats, 1: timestep stats, 2: all stats",
        false, 0, &statLoggingLevels, cmd);

    // Require one and only one of print, run, options, or generate modes
    cmd.parse(argc, argv);

    g_render = display.getValue();
    g_dumpcoord = output.getValue();
    g_dumpbinary = checkpoint.getValue();

    int idx = 0;
    g_ps = problems[idx];
    g_problem_idx = idx;
    setOptions();

    if (file.isSet()) {
      if (g_ps->LoadOptions(file.getValue()) == -1) {
        return -1;
      }
    } else {
      std::cout << "No XML file specified. Please specify XML file with -f"
                << std::endl;
      return -1;
    }

    if (statlog.isSet()) {
      int statLogLevel = statlog.getValue();
      g_ps->setStatLogging(statLogLevel);
    }
  } catch (TCLAP::ArgException& e) {
    std::cerr << "ERROR: " << e.argId() << std::endl
              << "       " << e.error() << std::endl;
  }

  return 0;
}

void printCommandLineSplashScreen() {
  assert(g_problem_idx >= 0);
  assert(g_problem_idx < (int)problems.size());

  std::cout << "\033[32;1m";

#ifdef DEBUG
  std::cout << "# build mode: DEBUG" << std::endl;
#else
  std::cout << "# build mode: RELEASE" << std::endl;
#endif

  std::cout << std::setfill('0');
  std::cout << "# timestamp: " << (1900 + g_timeinfo->tm_year) << "/"
            << std::setw(2) << (1 + g_timeinfo->tm_mon) << "/";
  std::cout << (g_timeinfo->tm_mday) << ", " << std::setw(2)
            << (g_timeinfo->tm_hour) << ":" << std::setw(2)
            << (g_timeinfo->tm_min);
  std::cout << ":" << std::setw(2) << (g_timeinfo->tm_sec) << std::endl;

  std::cout << "# problem: " << problems[g_problem_idx]->ProblemName()
            << std::endl;
  std::cout << "\033[m";
}

std::string generateOutputDirName() {
  assert(g_timeinfo != NULL);
  if (!g_ps) return "";

  std::stringstream datestream;
  datestream.fill('0');
  datestream << "simulation_capture_" << g_ps->ProblemName() << "_";
  datestream << std::setw(4) << (1900 + g_timeinfo->tm_year) << "_"
             << std::setw(2) << (1 + g_timeinfo->tm_mon) << "_";
  datestream << std::setw(2) << (g_timeinfo->tm_mday) << "_" << std::setw(2)
             << (g_timeinfo->tm_hour) << "_" << std::setw(2)
             << (g_timeinfo->tm_min) << "_";
  datestream << std::setw(2) << (g_timeinfo->tm_sec) << "/";

  return datestream.str();
}

int main(int argc, char** argv) {
  Eigen::initParallel();
  Eigen::setNbThreads(0);
  // make tme-stamped directory name
  time(&g_rawtime);
  g_timeinfo = localtime(&g_rawtime);

  createProblemVector();
  atexit(cleanup);

  if (parseCommandLine(argc, argv) < 0) {
    return -1;
  }

  g_ps->setup(g_current_frame, g_current_checkpoint);

  const int step = floor(g_ps->getTime() / g_ps->getDt() + 0.5);

  int step_per_checkpoint = g_dumpbinary;
  g_last_checkpoint_num = step / std::max(1, step_per_checkpoint);

  int steps_per_frame = g_dumpcoord;
  g_last_frame_num = step / std::max(1, steps_per_frame);

  g_outputdirectory = generateOutputDirName();

  g_ps->setOutputDirectory(g_outputdirectory);

  if (g_dumpcoord) {
    // dump copy of options used for the record
#ifdef WIN32
    const std::string base_path(_getcwd(NULL, 0));
    // Windows will eliminate "." at trailing of a folder name.
    // Use SHCreateDirectoryEx allows to create this such folder.
    const std::string g_outputdirectory_ = base_path + "\\" + g_outputdirectory;
    int sh_err_code =
        SHCreateDirectoryEx(NULL, g_outputdirectory_.c_str(), NULL);
    if (sh_err_code) {
      std::cout << "Create directory " << g_outputdirectory_
                << " failed, results will not be dumped. Error code"
                << sh_err_code << std::endl;
    }
#else
    if (!boost::filesystem::exists(g_outputdirectory)) {
      try {
        boost::filesystem::create_directory(g_outputdirectory);
      } catch (boost::filesystem::filesystem_error& e) {
        const boost::system::error_code errorCode = e.code();
        std::cout
            << "boost::filesystem::create_directories failed with error code: "
            << errorCode.message();
      }
      boost::filesystem::permissions(
          g_outputdirectory,
          boost::filesystem::add_perms | boost::filesystem::others_all |
              boost::filesystem::owner_all | boost::filesystem::group_all);
    }
#endif  // !WIN32

    const int file_width = 20;

    g_ps->dumpData(g_outputdirectory, g_current_frame, file_width);
    ++g_current_frame;

    g_ps->PrintAdditionalSettings(g_outputdirectory);
  }

  printCommandLineSplashScreen();

  if (g_render) {
    initializeOpenGL(argc, argv);
    std::thread(doSimulation).detach();
    glutMainLoop();
  } else {
    g_paused = false;
    for (;;) doSimulation();
  }

  return 0;
}

void doSimulation() {
  // If the simulation isn't paused, take a timestep
  while (!g_exit) {
    if (!g_paused) {
      if (!g_ps->step()) {
        g_paused = true;
      }
      output();

      if (g_single_step) {
        g_paused = true;
        g_single_step = false;
      }
    } else {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
  }
}
