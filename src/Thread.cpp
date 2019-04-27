// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Thread.cpp
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the implementation of the controller
//  class.
//					The purpose of the controller is to manage the
//simulation and 					its associated events. This includes controlling the simulation
//					time, managing events, controlling the
//active/inactive work- 					stations, and collecting/saving/printing the results of
//the 					simulation.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  06/02/2009	v1.02	Minor optimizations and bug fixes.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#include "Thread.h"

#include <sstream>

#ifndef NO_ALLEGRO

#include "allegro5/allegro.h"

extern ALLEGRO_BITMAP* mainbuf;
extern ALLEGRO_BITMAP* progbarbmp;
extern ALLEGRO_BITMAP* detailinfo;
extern ALLEGRO_BITMAP* popup;

#endif

extern Thread* threadZero;
extern Thread** threads;

extern size_t threadCount;

extern std::vector<AlgorithmToRun*> algParams;

const double Thread::TEN_HOURS = 10.0 * 60.0 * 60.0;
const double Thread::SPEED_OF_LIGHT = double(299792458);

///////////////////////////////////////////////////////////////////
//
// Function Name:	Thread
// Description:		Default constructor with no arguments, terminates
//					the program as command line arguments
//are 					required.
//
///////////////////////////////////////////////////////////////////
Thread::Thread()
    : CurrentActiveWorkstations(0),
      CurrentProbeStyle(ProbeStyle::NUMBER_OF_PROBE_STYLES),
      CurrentQualityAware(false),
      CurrentWavelengthAlgorithm(
          WavelengthAlgorithm::NUMBER_OF_WAVELENGTH_ALGORITHMS),
      globalTime(0.0),
      logger(nullptr),
      maxRunCount(0),
      maxSpans(0),
      minDuration(0.0),
      numOfWavelengths(0),
      numberOfConnections(0),
      numberOfEdges(0),
      numberOfRouters(0),
      numberOfWorkstations(0),
      order_init(false),
      qualityParams(),
      queue(nullptr),
      randomSeed(0),
      rm(nullptr),
      runCount(0),
      stats(),
      workstationOrder(nullptr),
      CurrentRoutingAlgorithm(RoutingAlgorithm::NUMBER_OF_ROUTING_ALGORITHMS) {
  threadZero->recordEvent(std::string("Unable to initialize the controller "
                                      "without command line arguments.\n"),
                          true, controllerIndex);
  exit(ERROR_THREAD_INIT);
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	Thread
// Description:		Constructor that takes the command line arguments
//					and initializes the controller.
//
///////////////////////////////////////////////////////////////////
Thread::Thread(size_t ci, int argc, const char* argv[], bool isLPS)
    : CurrentActiveWorkstations(0),
      CurrentProbeStyle(ProbeStyle::NUMBER_OF_PROBE_STYLES),
      CurrentQualityAware(false),
      CurrentWavelengthAlgorithm(
          WavelengthAlgorithm::NUMBER_OF_WAVELENGTH_ALGORITHMS),
      globalTime(0.0),
      logger(nullptr),
      maxRunCount(0),
      maxSpans(0),
      minDuration(0.0),
      numOfWavelengths(0),
      numberOfConnections(0),
      numberOfEdges(0),
      numberOfRouters(0),
      numberOfWorkstations(0),
      order_init(false),
      qualityParams(),
      queue(nullptr),
      randomSeed(0),
      rm(nullptr),
      runCount(0),
      stats(),
      workstationOrder(nullptr),
      CurrentRoutingAlgorithm(RoutingAlgorithm::NUMBER_OF_ROUTING_ALGORITHMS) {
  isLoadPrevious = isLPS;

  if (!isLoadPrevious) {
    threads[ci] = this;
  }

  if (ci == 0)
    threadZero = this;
  else
    setMinDuration(static_cast<size_t>(ceil(
        threadZero->maxSpans * threadZero->getQualityParams().QFactor_factor)));

#ifndef NO_ALLEGRO
  sprintf(topology, "%s", argv[1]);
  maxMaxUsage = 0.0;
  paintDir = false;
  paintSpans = false;
  paintRealTime = false;
  terminateProgram = false;
  numUsagesToDisplay = 0;  // initialize this to zero, but should later be set
                           // to usageList.size();
  usageHistStartTime =
      0;  // initialize this to zero. First hist graph should start at t = 0
  progBarLengthpx = 600;
  paintedConns = 0;
#endif

  controllerIndex = ci;
  setGlobalTime(0.0);

  std::string topology = argv[1];
  std::string wavelengths = argv[2];
  std::string seed = argv[3];
  std::string threads = argv[4];
  std::string iterations = argv[5];
  std::string probes = argv[6];

  generator = std::default_random_engine(std::stoi(seed));

  generateZeroToOne = std::uniform_real_distribution<double>(0, 1);

  setTopology(topology);

  if (controllerIndex == 0 && isLoadPrevious == false) {
    logger = new MessageLogger(topology, wavelengths, seed, probes);

    std::string quality =
        "input/Quality-" + topology + "-" + wavelengths + ".txt";
    setQualityParameters(quality);
  }
  else {
	logger = nullptr;
  }

  std::string topologyfile = "input/Topology-" + topology + ".txt";
  setTopologyParameters(topologyfile);

  std::string workstation =
      "input/Workstation-" + topology + "-" + wavelengths + ".txt";
  setWorkstationParameters(workstation);

  if (controllerIndex == 0 && isLoadPrevious == false) {
    qualityParams.max_probes = atoi(argv[6]);

    std::string algorithm = "input/Algorithm.txt";
    setAlgorithmParameters(algorithm, static_cast<size_t>(std::stoi(iterations)));

    numberOfConnections = static_cast<size_t>(TEN_HOURS / threadZero->getQualityParams().arrival_interval);
  }

  if (isLoadPrevious == false) {
    queue = new EventQueue();

    randomSeed = atoi(argv[3]);

    order_init = false;
    workstationOrder = new size_t[getNumberOfWorkstations()];
  }

  return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	~Thread
// Description:		Closes the Event Logger and deletes all of the
//					memory allocated by the controller
//
///////////////////////////////////////////////////////////////////
Thread::~Thread() {
  if (isLoadPrevious == false) {
    delete queue;
  }

  if (controllerIndex == 0 && isLoadPrevious == false) {
    delete[] qualityParams.ASE_perEDFA;

    delete logger;
    delete rm;
  }

  for (size_t r = 0; r < routers.size(); ++r) delete routers[r];

  for (size_t w = 0; w < workstations.size(); ++w) delete workstations[w];

  routers.clear();
  workstations.clear();

  if (isLoadPrevious == false) {
    delete[] workstationOrder;
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	initPriorityQueue
// Description:		Creates a new Event Queue and initalizes several
//					events in it.
//
///////////////////////////////////////////////////////////////////
void Thread::initPriorityQueue() {
  Event* activate = new Event;
  activate->e_type = ACTIVATE_WORKSTATIONS;
  activate->e_time = 0.0;
  activate->e_data = nullptr;

  Event* deactivate = new Event;
  deactivate->e_type = DEACTIVATE_WORKSTATIONS;
  deactivate->e_time = std::numeric_limits<double>::infinity();
  deactivate->e_data = nullptr;

  queue->addEvent(*deactivate);
  queue->addEvent(*activate);

  delete deactivate;
  delete activate;

  if (CurrentRoutingAlgorithm == PABR || CurrentRoutingAlgorithm == LORA) {
    Event* event = new Event();

    event->e_type = UPDATE_USAGE;
    event->e_time = 0.0;
    event->e_data = nullptr;

    queue->addEvent(*event);

    delete event;
  }

#ifndef NO_ALLEGRO
  Event* event = new Event();

  event->e_type = UPDATE_GUI;
  event->e_time = 0.0;
  event->e_data = nullptr;

  queue->addEvent(*event);

  delete event;
#endif
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	initResourceManager
// Description:		Creates a new Resource Manager and initalizes several
//					events in it.
//
///////////////////////////////////////////////////////////////////
void Thread::initResourceManager() {
  if (controllerIndex == 0)
    rm = new ResourceManager();
  else
    rm = nullptr;

  for (size_t r = 0; r < getNumberOfRouters(); ++r) {
    if (threadZero->getQualityParams().dest_dist != UNIFORM) {
      getRouterAt(r)->generateProbabilities();
    }
  }
}
///////////////////////////////////////////////////////////////////
//
// Function Name:	initThread
// Description:		Runs the simulation for the program until there
//					are no more events to run.
//
///////////////////////////////////////////////////////////////////
void Thread::initThread(AlgorithmToRun* alg) {
  CurrentRoutingAlgorithm = alg->ra;
  CurrentWavelengthAlgorithm = alg->wa;
  CurrentProbeStyle = alg->ps;
  CurrentQualityAware = alg->qa;
  CurrentActiveWorkstations = alg->workstations;

#ifndef NO_ALLEGRO
  rectfill(mainbuf, 0, 50 * controllerIndex + 85 - 1, SCREEN_W,
           50 * (controllerIndex + 1) + 85 - 1, makecol(0, 0, 0));
  char buffer[200];
  sprintf(buffer, "THREAD %d: RA = %s, WA = %s, PS = %s, QA = %d, N = %d",
          controllerIndex,
          threadZero->getRoutingAlgorithmName(CurrentRoutingAlgorithm)->c_str(),
          threadZero->getWavelengthAlgorithmName(CurrentWavelengthAlgorithm)
              ->c_str(),
          threadZero->getProbeStyleName(CurrentProbeStyle)->c_str(),
          (int)CurrentQualityAware, CurrentActiveWorkstations);
  textprintf_ex(mainbuf, font, 40, 50 * controllerIndex + 85,
                makecol(0, 255, 0), -1, buffer);
#endif

  initPriorityQueue();

  delete alg;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	runThread
// Description:		Runs the simulation for the program until there
//					are no more events to run.
//
///////////////////////////////////////////////////////////////////
int Thread::runThread(AlgorithmToRun* alg) {
  initThread(alg);

#ifndef NO_ALLEGRO
  rectfill(mainbuf, 49, controllerIndex * 50 + 99, 651,
           controllerIndex * 50 + 99 + 26, makecol(0, 0, 255));
#endif

  while (queue->getSize() > 0) {
#ifndef NO_ALLEGRO
    if (terminateProgram == true) {
      threadZero->recordEvent(
          "ERROR: User clicked on the close button, exiting simulation.", true,
          controllerIndex);
      exit(ERROR_USER_CLOSED);
    }
#endif

    Event event = queue->getNextEvent();

    setGlobalTime(event.e_time);

#ifndef NO_ALLEGRO  // PROGRESS BAR

    int prog = stats.ConnectionRequests * multFactor;
    if (prog >= 1) {
      if ((stats.ConnectionRequests - paintedConns) > (ConnsPerPx * 5) &&
          stats.ConnectionRequests > 0) {
        stretch_blit(progbarbmp, mainbuf, 0, 0, progbarbmp->w, progbarbmp->h,
                     50, controllerIndex * 50 + 100, prog,
                     progbarbmp->h);  // new
        masked_blit(mainbuf, screen, 0, 0, 0, 0, SCREEN_W, SCREEN_H);
        paintedConns = stats.ConnectionRequests;
      }
    }
#endif

    switch (event.e_type) {
      case ACTIVATE_WORKSTATIONS:
        activate_workstations();
        break;
      case DEACTIVATE_WORKSTATIONS:
        deactivate_workstations();
        break;
      case UPDATE_USAGE:
        update_link_usage();
        break;
#ifndef NO_ALLEGRO
      case UPDATE_GUI:
        update_gui();
        break;
#endif
      case CONNECTION_REQUEST:
        connection_request(static_cast<ConnectionRequestEvent*>(event.e_data));
        delete static_cast<ConnectionRequestEvent*>(event.e_data);
        break;
      case COLLISION_NOTIFICATION:
        collision_notification(
            static_cast<CollisionNotificationEvent*>(event.e_data));
        break;
      case CREATE_CONNECTION_PROBE:
        create_connection_probe(
            static_cast<CreateConnectionProbeEvent*>(event.e_data));
        break;
      case CREATE_CONNECTION_CONFIRMATION:
        create_connection_confirmation(
            static_cast<CreateConnectionConfirmationEvent*>(event.e_data));
        break;
      case DESTROY_CONNECTION_PROBE:
        destroy_connection_probe(
            static_cast<DestroyConnectionProbeEvent*>(event.e_data));
        break;
      default:
        threadZero->recordEvent(
            "ERROR: Unknown event type, exiting simulation.", true,
            controllerIndex);
        exit(ERROR_THREAD_EVENT_TYPE);
        break;
    }
  }

  return 0;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	activateWorkstations
// Description:		Activates the specified number of workstations
//					randomly
//
///////////////////////////////////////////////////////////////////
void Thread::activate_workstations() {
#ifndef NO_ALLEGRO
  ConnsPerPx =
      (threadZero->getNumberOfConnections() * getCurrentActiveWorkstations()) /
      progBarLengthpx;
  multFactor = 600.0 / double(getCurrentActiveWorkstations() *
                              threadZero->getNumberOfConnections());
#endif

  stats.ConnectionRequests = 0;
  stats.ConnectionSuccesses = 0;
  stats.CollisionFailures = 0;
  stats.NoPathFailures = 0;
  stats.QualityFailures = 0;
  stats.DroppedFailures = 0;
  stats.ProbeSentCount = 0;

  stats.totalSetupDelay = 0.0;
  stats.totalHopCount = 0;
  stats.totalSpanCount = 0;

  stats.aseNoiseTotal = 0.0;
  stats.fwmNoiseTotal = 0.0;
  stats.xpmNoiseTotal = 0.0;
  stats.raRunTime = 0.0;

  // Random generator for destination router
  generateRandomRouter =
      std::uniform_int_distribution<size_t>(0, getNumberOfRouters() - 1);

  // Random generator for duration time
  generateRandomDuration = std::exponential_distribution<double>(
      1.0 / double(threadZero->getQualityParams().duration));

  // Random generator for arrival interval
  generateArrivalInterval = std::exponential_distribution<double>(
      1.0 / double(threadZero->getQualityParams().arrival_interval));

  if (CurrentRoutingAlgorithm == SHORTEST_PATH) {
    threadZero->getResourceManager()->initSPMatrix();
  } else {
    bool SP_active = false;

    for (size_t t = 0; t < threadCount; ++t) {
      if (threads[t]->getCurrentRoutingAlgorithm() == SHORTEST_PATH) {
        SP_active = true;
        break;
      }
    }

    if (SP_active == false) threadZero->getResourceManager()->freeSPMatrix();
  }

  if (CurrentRoutingAlgorithm == ADAPTIVE_QoS) {
    for (size_t r = 0; r < getNumberOfRouters(); ++r) {
      getRouterAt(r)->resetFailures();
    }
  }

  if (order_init == false) {
    std::uniform_int_distribution<size_t> generateRandomWorkstation =
        std::uniform_int_distribution<size_t>(0, getNumberOfWorkstations() - 1);

    bool* init = new bool[getNumberOfWorkstations()];

    for (size_t w = 0; w < getNumberOfWorkstations(); ++w) {
      init[w] = false;
      workstationOrder[w] = 0;
    }

    size_t numberFound = 0;

    while (numberFound < getNumberOfWorkstations()) {
      size_t wkstn = generateRandomWorkstation(generator);

      if (init[wkstn] == false) {
        workstationOrder[numberFound] = wkstn;
        init[wkstn] = true;

        numberFound++;
      }
    }

    delete[] init;

    order_init = true;
  }

  std::ostringstream buffer;
  buffer << "Activate " << getCurrentActiveWorkstations() << "workstations.";
  threadZero->recordEvent(buffer.str(), false, controllerIndex);

  for (size_t w = 0; w < getCurrentActiveWorkstations(); ++w) {
    getWorkstationAt(workstationOrder[w])->setActive(true);
#ifndef NO_ALLEGRO
    getRouterAt(getWorkstationAt(workstationOrder[w])->getParentRouterIndex())
        ->incNumWorkstations();
#endif
    buffer.clear();
    buffer << "\tActivate workstation " << workstationOrder[w];
    threadZero->recordEvent(buffer.str(), false, controllerIndex);
  }

  if (CurrentRoutingAlgorithm == PABR || CurrentRoutingAlgorithm == LORA) {
    for (size_t r = 0; r < getNumberOfRouters(); ++r) {
      getRouterAt(r)->resetUsage();
    }
  } else if (CurrentRoutingAlgorithm == Q_MEASUREMENT ||
             CurrentRoutingAlgorithm == ADAPTIVE_QoS) {
    for (size_t r = 0; r < getNumberOfRouters(); ++r) {
      getRouterAt(r)->resetQMDegredation();
    }
  }

  for (size_t w = 0; w < getNumberOfWorkstations(); ++w)
    if (getWorkstationAt(w)->getActive() == true)
      generateTrafficEvent(w * threadZero->getNumberOfConnections());
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	deactivateWorkstations
// Description:		Deactivates all of the workstations and prints
//					statistics on their performance.
//
///////////////////////////////////////////////////////////////////
void Thread::deactivate_workstations() {
  threadZero->getLogger()->LockResultsMutex();

  std::ostringstream algorithm;
  algorithm << "**ALGORITHM = "
            << threadZero->getRoutingAlgorithmName(CurrentRoutingAlgorithm)
            << "-"
            << threadZero->getWavelengthAlgorithmName(
                   CurrentWavelengthAlgorithm)
            << ", WORKS = " << getCurrentActiveWorkstations()
            << ", PROBE = " << threadZero->getProbeStyleName(CurrentProbeStyle)
            << ", QA = " << getCurrentQualityAware();
  threadZero->recordEvent(algorithm.str(), true, controllerIndex);

#ifndef NO_ALLEGRO
  strcpy(routing,
         threadZero->getRoutingAlgorithmName(CurrentRoutingAlgorithm)->c_str());
  strcpy(wavelength,
         threadZero->getWavelengthAlgorithmName(CurrentWavelengthAlgorithm)
             ->c_str());
  strcpy(probing, threadZero->getProbeStyleName(CurrentProbeStyle)->c_str());
  qualityAware = getCurrentQualityAware();

  overallBlocking =
      double(stats.ConnectionRequests - stats.ConnectionSuccesses) /
      double(stats.ConnectionRequests);
  collisions =
      double(stats.CollisionFailures) / double(stats.ConnectionRequests);
  badquality = double(stats.QualityFailures) / double(stats.ConnectionRequests);
  nonresource = double(stats.NoPathFailures) / double(stats.ConnectionRequests);
  avgprobesperrequest =
      double(stats.ProbeSentCount) / double(stats.ConnectionRequests);
  avgrequestdelaytime =
      stats.totalSetupDelay / double(stats.ConnectionSuccesses);
  avgconnhopcount =
      double(stats.totalHopCount) / double(stats.ConnectionSuccesses);
  avgconnspancount =
      double(stats.totalSpanCount) / double(stats.ConnectionSuccesses);
  avgASEnoise = stats.aseNoiseTotal / double(stats.ConnectionSuccesses);
  avgFWMnoise = stats.fwmNoiseTotal / double(stats.ConnectionSuccesses);
  avgXPMnoise = stats.xpmNoiseTotal / double(stats.ConnectionSuccesses);
#endif

  std::ostringstream overall;
  overall << "OVERALL BLOCKING ("
          << stats.ConnectionRequests - stats.ConnectionSuccesses << "/"
          << stats.ConnectionRequests << ") = "
          << double(stats.ConnectionRequests - stats.ConnectionSuccesses) /
                 double(stats.ConnectionRequests);
  threadZero->recordEvent(overall.str(), true, controllerIndex);

  std::ostringstream collisions;
  collisions << "COLLISIONS (" << stats.CollisionFailures << "/"
             << stats.ConnectionRequests << ") = "
             << double(stats.CollisionFailures) /
                    double(stats.ConnectionRequests);
  threadZero->recordEvent(collisions.str(), true, controllerIndex);

  std::ostringstream quality;
  quality << "BAD QUALITY (" << stats.QualityFailures << "/"
          << stats.ConnectionRequests << ") = "
          << double(stats.QualityFailures) / double(stats.ConnectionRequests);
  threadZero->recordEvent(quality.str(), true, controllerIndex);

  std::ostringstream resources;
  resources << "NON RESOURCES (" << stats.NoPathFailures << "/"
            << stats.ConnectionRequests << ") = "
            << double(stats.NoPathFailures) / double(stats.ConnectionRequests);
  threadZero->recordEvent(resources.str(), true, controllerIndex);

  std::ostringstream probes;
  probes << "AVERAGE PROBES PER REQUEST (" << stats.ProbeSentCount << "/"
         << stats.ConnectionRequests << ") = "
         << double(stats.ProbeSentCount) / double(stats.ConnectionRequests);
  threadZero->recordEvent(probes.str(), true, controllerIndex);

  std::ostringstream delay;
  delay << "AVERAGE REQUEST DELAY TIME (" << stats.totalSetupDelay << "/"
        << stats.ConnectionSuccesses
        << ") = " << stats.totalSetupDelay / double(stats.ConnectionSuccesses);
  threadZero->recordEvent(delay.str(), true, controllerIndex);

  std::ostringstream hopcount;
  hopcount << "AVERAGE CONNECTION HOP COUNT (" << stats.totalHopCount << "/"
           << stats.ConnectionSuccesses << ") = "
           << double(stats.totalHopCount) / double(stats.ConnectionSuccesses);
  threadZero->recordEvent(hopcount.str(), true, controllerIndex);

  std::ostringstream connection;
  connection << "AVERAGE CONNECTION SPAN COUNT (" << stats.totalSpanCount << "/"
             << stats.ConnectionSuccesses << ") = "
             << double(stats.totalSpanCount) /
                    double(stats.ConnectionSuccesses);
  threadZero->recordEvent(connection.str(), true, controllerIndex);

  std::ostringstream ase;
  ase << "AVERAGE ASE NOISE (" << stats.aseNoiseTotal << "/"
      << stats.ConnectionSuccesses
      << ") = " << stats.aseNoiseTotal / double(stats.ConnectionSuccesses);
  threadZero->recordEvent(ase.str(), true, controllerIndex);

  std::ostringstream fwm;
  fwm << "AVERAGE FWM NOISE (" << stats.fwmNoiseTotal << "/"
      << stats.ConnectionSuccesses
      << ") = " << stats.fwmNoiseTotal / double(stats.ConnectionSuccesses);
  threadZero->recordEvent(fwm.str(), true, controllerIndex);

  std::ostringstream xpm;
  xpm << "AVERAGE XPM NOISE (" << stats.xpmNoiseTotal << "/"
      << stats.ConnectionSuccesses
      << ") = " << stats.xpmNoiseTotal / double(stats.ConnectionSuccesses);
  threadZero->recordEvent(xpm.str(), true, controllerIndex);

  std::ostringstream runtime;
  runtime << "AVERAGE RA RUN TIME (" << stats.raRunTime << "/"
          << stats.ConnectionRequests
          << ") = " << stats.raRunTime / double(stats.ConnectionRequests);
  threadZero->recordEvent(runtime.str(), true, controllerIndex);

  if (threadZero->getQualityParams().q_factor_stats == true) {
    double worstInitQ = std::numeric_limits<double>::infinity();
    double bestInitQ = 0.0;
    double averageInitQ = 0.0;

    double worstAvgQ = std::numeric_limits<double>::infinity();
    double bestAvgQ = 0.0;
    double averageAvgQ = 0.0;

    double worstPerQ = std::numeric_limits<double>::infinity();
    double bestPerQ = 0.0;
    double averagePerQ = 0.0;

    double timeTotal = 0.0;
    double countTotal = 0.0;
    double droppedTotal = 0.0;

    for (size_t r = 0; r < threadZero->getNumberOfRouters(); ++r) {
      Router* router = getRouterAt(r);

      for (size_t e = 0; e < router->getNumberOfEdges(); ++e) {
        EdgeStats* stats = router->getEdgeByIndex(e)->getEdgeStats();

        averageInitQ += stats->totalInitalQFactor;
        averageAvgQ += stats->totalAverageQFactor;
        averagePerQ += stats->totalPercentQFactor;

        timeTotal += stats->totalTime;
        countTotal += stats->count;

        droppedTotal += stats->droppedConnections;

        if (stats->minInitalQFactor < worstInitQ)
          worstInitQ = stats->minInitalQFactor;
        else if (stats->maxInitalQFactor > bestInitQ)
          bestInitQ = stats->maxInitalQFactor;

        if (stats->minAverageQFactor < worstAvgQ)
          worstAvgQ = stats->minAverageQFactor;
        else if (stats->maxAverageQFactor > bestAvgQ)
          bestAvgQ = stats->maxAverageQFactor;

        if (stats->minPercentQFactor < worstPerQ)
          worstPerQ = stats->minPercentQFactor;
        else if (stats->maxPercentQFactor > bestPerQ)
          bestPerQ = stats->maxPercentQFactor;

        router->getEdgeByIndex(e)->resetEdgeStats();
      }
    }

    std::ostringstream dropped;
    dropped << "DROPPED CONNECTIONS (" << int(droppedTotal) << "/"
            << stats.ConnectionRequests
            << ") = " << droppedTotal / double(stats.ConnectionRequests);
    threadZero->recordEvent(dropped.str(), true, controllerIndex);

    std::ostringstream odropped;
    odropped << "OVERALL W/DROPPED ("
             << int(stats.ConnectionRequests - stats.ConnectionSuccesses +
                    droppedTotal)
             << "/" << stats.ConnectionRequests << ") = "
             << double(stats.ConnectionRequests - stats.ConnectionSuccesses +
                       droppedTotal) /
                    double(stats.ConnectionRequests);
    threadZero->recordEvent(odropped.str(), true, controllerIndex);

    std::ostringstream initialq;
    initialq << "INITIAL Q: MIN = " << worstInitQ << ", MAX = " << bestInitQ
             << ", AVG = " << averageInitQ / countTotal;
    threadZero->recordEvent(initialq.str(), true, controllerIndex);

    std::ostringstream averageq;
    averageq << "AVERAGE Q: MIN = " << worstAvgQ << ", MAX  = " << bestAvgQ
             << ", AVG = " << averageAvgQ / countTotal;
    threadZero->recordEvent(averageq.str(), true, controllerIndex);

    std::ostringstream timebelow;
    timebelow << "%% TIME Q BELOW : MIN = " << worstPerQ
              << ", MAX = " << bestPerQ
              << ", AVG = " << averagePerQ / timeTotal;
    threadZero->recordEvent(timebelow.str(), true, controllerIndex);
  }

  std::ostringstream stars;
  stars << std::string("**********************************************")
        << std::endl;
  threadZero->recordEvent(stars.str(), true, controllerIndex);

  threadZero->flushLog(true);

  threadZero->getLogger()->UnlockResultsMutex();

  for (size_t r1 = 0; r1 < getNumberOfRouters(); ++r1) {
    for (size_t r2 = 0; r2 < getNumberOfRouters(); ++r2) {
      Edge* edge = getRouterAt(r1)->getEdgeByDestination(r2);

      if (edge != 0) {
        for (size_t k = 0; k < threadZero->getNumberOfWavelengths();
             ++k) {
          if (edge->getStatus(k) != EDGE_FREE) {
            threadZero->recordEvent(
                "ERROR: Edge is still used when all edges should be free.",
                true, controllerIndex);
            exit(ERROR_EDGE_IS_USED);
          }
        }
      }
    }
  }

  threadZero->recordEvent("Deactivate workstations.", false, controllerIndex);

  for (size_t w = 0; w < getNumberOfWorkstations(); ++w) {
    if (getWorkstationAt(w)->getActive() == true) {
      std::ostringstream deactivate;
      deactivate << "\tDeactivate workstation " << w;
      threadZero->recordEvent(deactivate.str(), false, controllerIndex);

      getWorkstationAt(w)->setActive(false);
    }
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	generateTrafficEvent
// Description:		Generates connection request events to setup requests
// for new
//					connections.
//
///////////////////////////////////////////////////////////////////
void Thread::generateTrafficEvent(size_t session) {
  size_t workstation = session / threadZero->getNumberOfConnections();

  Event* tr = new Event();
  ConnectionRequestEvent* tr_data = new ConnectionRequestEvent();

  tr->e_time = getGlobalTime() + generateArrivalInterval(generator);
  tr->e_type = CONNECTION_REQUEST;
  tr->e_data = tr_data;

  tr_data->connectionDuration = generateRandomDuration(generator);

  if (tr_data->connectionDuration < threadZero->getMinDuration())
    tr_data->connectionDuration = threadZero->getMinDuration();

  tr_data->requestBeginTime = tr->e_time;
  tr_data->sourceRouterIndex =
      getWorkstationAt(workstation)->getParentRouterIndex();
  tr_data->destinationRouterIndex = tr_data->sourceRouterIndex;
  tr_data->session = session;
  tr_data->sequence = 0;
  tr_data->max_sequence = 0;
  tr_data->qualityFail = false;

  while (tr_data->sourceRouterIndex == tr_data->destinationRouterIndex) {
    if (threadZero->getQualityParams().dest_dist == UNIFORM)
      tr_data->destinationRouterIndex = (generateRandomRouter(generator));
    else if (threadZero->getQualityParams().dest_dist == DISTANCE ||
             threadZero->getQualityParams().dest_dist == INVERSE_DISTANCE)
      tr_data->destinationRouterIndex =
          getRouterAt(tr_data->sourceRouterIndex)
              ->generateDestination(generateZeroToOne(generator));
  }

#ifndef NO_ALLEGRO
  routers[tr_data->sourceRouterIndex]->incConnAttemptsFrom();
  routers[tr_data->destinationRouterIndex]->incConnAttemptsTo();
#endif

  queue->addEvent(*tr);

  delete tr;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	connection_request
// Description:		Handles a connection request event
//
///////////////////////////////////////////////////////////////////
void Thread::connection_request(ConnectionRequestEvent* cre) {
  if (getGlobalTime() < TEN_HOURS) {
    if ((cre->session + 1) % threadZero->getNumberOfConnections() != 0)
      generateTrafficEvent(cre->session + 1);
    else
      generateTrafficEvent(cre->session + 1 -
                           threadZero->getNumberOfConnections());
  }

  ++stats.ConnectionRequests;

  if (CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING ||
      CurrentRoutingAlgorithm == IMPAIRMENT_AWARE) {
    ++stats.ProbeSentCount;
  }

  long long int probesToSend = 0;
  long long int probeStart = 0;
  long long int probesSkipped = 0;

  if (CurrentProbeStyle == SINGLE)
    probesToSend = 1;
  else
    probesToSend = threadZero->getQualityParams().max_probes;

  kShortestPathReturn* kPath = getKShortestPaths(cre, probesToSend);

  CreateConnectionProbeEvent** probesList =
      calcProbesToSend(cre, kPath, probesToSend, probeStart, probesSkipped);

  sendProbes(cre, kPath, probesList, probesToSend, probeStart, probesSkipped);
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	create_connection_probe
// Description:		Handles the probe message for creating a new
//					connection.
//
///////////////////////////////////////////////////////////////////
void Thread::create_connection_probe(CreateConnectionProbeEvent* ccpe) {
  ++ccpe->numberOfHops;

  if (ccpe->numberOfHops == ccpe->connectionLength) {
    ccpe->atDestination = true;

    if (CurrentProbeStyle != PARALLEL || sendResponse(ccpe) == true) {
      double q_factor = 0.0;
      double xpm_noise = 0.0;
      double fwm_noise = 0.0;
      double ase_noise = 0.0;

      Event* event = new Event();
      CreateConnectionConfirmationEvent* ccce =
          new CreateConnectionConfirmationEvent;

      ccce->connectionDuration = ccpe->connectionDuration;
      ccce->requestBeginTime = ccpe->requestBeginTime;
      ccce->connectionLength = ccpe->connectionLength;
      ccce->connectionPath = ccpe->connectionPath;
      ccce->destinationRouterIndex = ccpe->destinationRouterIndex;
      ccce->numberOfHops = 0;
      ccce->sourceRouterIndex = ccpe->sourceRouterIndex;
      ccce->kPaths = ccpe->kPaths;
      ccce->session = ccpe->session;
      ccce->sequence = ccpe->sequence;
      ccce->max_sequence = ccpe->max_sequence;
      ccce->originalWavelength = ccpe->wavelength;
      ccce->probes = ccpe->probes;
      ccce->finalFailure = false;

      ccpe->decisionTime = getGlobalTime();

      // Again...some algorithms have to be treated differently because they
      // use a forward reservation scheme.
      if (CurrentRoutingAlgorithm == IMPAIRMENT_AWARE ||
          CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING) {
        double q_factor = 0.0;
        double xpm_noise = 0.0;
        double fwm_noise = 0.0;
        double ase_noise = 0.0;

        if (ccpe->wavelength >= 0) {
          for (size_t p = 0; p < ccpe->connectionLength; ++p) {
            if (ccpe->connectionPath[p]->getStatus(ccpe->wavelength) !=
                EDGE_FREE) {
              ccpe->wavelength = NO_PATH_FAILURE;
              break;
            }
          }

          if (ccpe->wavelength != NO_PATH_FAILURE) {
            q_factor = threadZero->getResourceManager()->estimate_Q(
                ccpe->wavelength, ccpe->connectionPath, ccpe->connectionLength,
                &xpm_noise, &fwm_noise, &ase_noise, controllerIndex);

            if (getCurrentQualityAware() == true &&
                q_factor < threadZero->getQualityParams().TH_Q)
              ccpe->wavelength = QUALITY_FAILURE;
          }
        }

        if (ccpe->wavelength >= 0 ||
            CurrentRoutingAlgorithm == IMPAIRMENT_AWARE) {
          threadZero->getResourceManager()->print_connection_info(
              ccpe, q_factor, ase_noise, fwm_noise, xpm_noise, controllerIndex);
        } else {
          time_t start;
          time_t end;

          time(&start);

          ccpe->wavelength =
              threadZero->getResourceManager()->choose_wavelength(
                  ccpe, controllerIndex);

          time(&end);
          stats.raRunTime += difftime(end, start);
        }
      } else {
        time_t start;
        time_t end;

        time(&start);

        ccpe->wavelength = threadZero->getResourceManager()->choose_wavelength(
            ccpe, controllerIndex);

        time(&end);
        stats.raRunTime += difftime(end, start);
      }

      if (CurrentProbeStyle == PARALLEL &&
          ccpe->wavelength == NO_PATH_FAILURE) {
        for (long long int p = 0; p < ccce->max_sequence; ++p) {
          if (ccce->probes[p]->wavelength == QUALITY_FAILURE) {
            ccpe->wavelength = QUALITY_FAILURE;
            break;
          }
        }
      } else if (CurrentProbeStyle == SERIAL &&
                 ccpe->wavelength == NO_PATH_FAILURE) {
        if (ccpe->qualityFail == true) ccpe->wavelength = QUALITY_FAILURE;
      }

      ccce->wavelength = ccpe->wavelength;

      if ((ccpe->wavelength == QUALITY_FAILURE ||
           ccpe->wavelength == NO_PATH_FAILURE) &&
          otherResponse(ccpe) != -1) {
        long long int sequence = otherResponse(ccpe);

        --ccpe->probes[sequence]->numberOfHops;

        create_connection_probe(ccpe->probes[sequence]);

        delete ccce;
      } else if ((ccpe->wavelength == QUALITY_FAILURE ||
                  ccpe->wavelength == NO_PATH_FAILURE) &&
                 moreProbes(ccpe) == true) {
        // Don't send rejection, just wait for more probes
        delete ccce;
      } else {
        event->e_type = CREATE_CONNECTION_CONFIRMATION;
        event->e_time =
            getGlobalTime() +
            calculateDelay(ccpe->connectionPath[ccpe->numberOfHops - 1]
                               ->getNumberOfSpans());
        event->e_data = ccce;

        queue->addEvent(*event);
      }

      if (CurrentProbeStyle != PARALLEL) delete ccpe;

      delete event;
    }
  } else {
    // We are not at our destination yet, so continue with the probe.
    Event* event = new Event();

    event->e_type = CREATE_CONNECTION_PROBE;
    event->e_time =
        getGlobalTime() +
        calculateDelay(
            ccpe->connectionPath[ccpe->numberOfHops]->getNumberOfSpans());
    event->e_data = ccpe;

    queue->addEvent(*event);

    delete event;
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	destroy_connection_probe
// Description:		Handles the probe message for creating a new
//					connection.
//
///////////////////////////////////////////////////////////////////
void Thread::destroy_connection_probe(DestroyConnectionProbeEvent* dcpe) {
  if (dcpe->numberOfHops == 0) {
    if (CurrentRoutingAlgorithm == Q_MEASUREMENT ||
        CurrentRoutingAlgorithm == ADAPTIVE_QoS) {
      for (size_t p = 0; p < dcpe->connectionLength; ++p)
        dcpe->connectionPath[p]->removeEstablishedConnection(dcpe);

      updateQMDegredation(dcpe->connectionPath, dcpe->connectionLength,
                          dcpe->wavelength);
    } else if (threadZero->getQualityParams().q_factor_stats == true) {
      for (size_t p = 0; p < dcpe->connectionLength; ++p)
        dcpe->connectionPath[p]->removeEstablishedConnection(dcpe);

      updateQFactorStats(dcpe->connectionPath, dcpe->connectionLength,
                         dcpe->wavelength);
    }
  }

  Edge* edge = dcpe->connectionPath[dcpe->numberOfHops];

  if (edge->getStatus(dcpe->wavelength) == EDGE_USED) {
    edge->setFree(dcpe->wavelength);
  } else {
    size_t srcIndex = edge->getSourceIndex();
    size_t destIndex = edge->getDestinationIndex();

    std::ostringstream buffer;
    buffer << "ERROR: Something is wrong, an edge from " << srcIndex << " to "
           << destIndex << " on wave " << dcpe->wavelength
           << " shou;d be marked as used is actually free already.";
    threadZero->recordEvent(buffer.str(), true, controllerIndex);
    exit(ERROR_EDGE_IS_FREE);
  }

  ++dcpe->numberOfHops;

  if (dcpe->numberOfHops == dcpe->connectionLength) {
    std::string line;

    line.append("SETUP CONNECTION: Routers[" +
                std::to_string(dcpe->connectionLength) + "]:");

    for (size_t r = 0; r < dcpe->connectionLength; ++r) {
      line.append(std::to_string(dcpe->connectionPath[r]->getSourceIndex()));

      if (r < dcpe->connectionLength) line.append(",");
    }

    line.append(std::to_string(dcpe->connectionPath[dcpe->connectionLength - 1]
                                   ->getDestinationIndex()));

    line.append(" Wavelength = " + std::to_string(dcpe->wavelength));
    line.append(" Session = " + std::to_string(dcpe->session));

    threadZero->recordEvent(line, false, controllerIndex);

    if (CurrentProbeStyle == PARALLEL)
      clearResponses(dcpe->probes[dcpe->sequence]);

    delete[] dcpe->connectionPath;
    delete dcpe;
  } else {
    // We are not at our destination yet, so continue with the probe.
    Event* event = new Event();

    event->e_type = DESTROY_CONNECTION_PROBE;
    event->e_time =
        getGlobalTime() +
        calculateDelay(
            dcpe->connectionPath[dcpe->numberOfHops]->getNumberOfSpans());
    event->e_data = dcpe;

    queue->addEvent(*event);

    delete event;
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	create_connection_confirmation
// Description:		Handles the probe message for creating a new
//					connection.
//
///////////////////////////////////////////////////////////////////
void Thread::create_connection_confirmation(
    CreateConnectionConfirmationEvent* ccce) {
#ifndef NO_ALLEGRO
  int indx;
#endif

  Edge* edge =
      ccce->connectionPath[ccce->connectionLength - ccce->numberOfHops - 1];

  if (ccce->wavelength >= 0 && edge->getStatus(ccce->wavelength) == EDGE_FREE) {
    // No collision yet, so reserve the link.
    edge->setUsed(ccce->session, ccce->wavelength);
  } else if (ccce->wavelength >= 0 &&
             edge->getStatus(ccce->wavelength) == EDGE_USED) {
    // Oops...a collision. We need to send a response upstream to notify the
    // routers to release the resources that were scheduled for this event.
    std::ostringstream buffer;
    buffer << "COLLISION: Session = " << ccce->session
           << ", Sequence = " << ccce->sequence;
    threadZero->recordEvent(buffer.str(), false, controllerIndex);

    Event* event = new Event();
    CollisionNotificationEvent* cne = new CollisionNotificationEvent();

    event->e_type = COLLISION_NOTIFICATION;
    event->e_time = getGlobalTime() + calculateDelay(edge->getNumberOfSpans());
    event->e_data = cne;

    cne->sourceRouterIndex = edge->getSourceIndex();
    cne->destinationRouterIndex = ccce->destinationRouterIndex;
    cne->wavelength = ccce->wavelength;
    cne->numberOfHops = ccce->connectionLength - ccce->numberOfHops - 1;
    cne->connectionLength = ccce->connectionLength;
    cne->connectionPath = new Edge*[cne->connectionLength];
    cne->session = ccce->session;
    cne->sequence = ccce->sequence;
    cne->probes = ccce->probes;

    if (CurrentProbeStyle == PARALLEL)
      cne->max_sequence = ccce->probes[ccce->sequence]->max_sequence;
    else
      cne->max_sequence = 1;

    for (size_t p = 0; p < ccce->connectionLength; ++p)
      cne->connectionPath[p] = ccce->connectionPath[p];

    cne->finalFailure = true;

    if (CurrentProbeStyle == PARALLEL) {
      for (long long int p = 0; p < ccce->probes[ccce->sequence]->max_sequence;
           ++p) {
        if (p != ccce->sequence) {
          if (ccce->probes[p]->atDestination == false) {
            cne->finalFailure = false;
            break;
          } else {
            if (ccce->probes[p]->decisionTime == 0.0 ||
                ccce->probes[p]->decisionTime >
                    ccce->probes[ccce->sequence]->decisionTime) {
              cne->finalFailure = false;
              break;
            }
          }
        }
      }
    }

    ccce->finalFailure = cne->finalFailure;

    queue->addEvent(*event);

    delete event;

    ccce->wavelength = COLLISION_FAILURE;
  }

  ++ccce->numberOfHops;

  if (ccce->connectionLength > ccce->numberOfHops) {
    // Forward the confirmation downstream.
    Event* event = new Event();

    event->e_type = CREATE_CONNECTION_CONFIRMATION;
    event->e_time = getGlobalTime() + calculateDelay(edge->getNumberOfSpans());
    event->e_data = ccce;

    queue->addEvent(*event);

    delete event;
  } else if (ccce->connectionLength == ccce->numberOfHops) {
    if (ccce->wavelength >= 0) {
      ++stats.ConnectionSuccesses;

      if (CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING ||
          CurrentRoutingAlgorithm == IMPAIRMENT_AWARE) {
        stats.totalSetupDelay +=
            (threadZero->getResourceManager()
                 ->span_distance[ccce->sourceRouterIndex *
                                 threadZero->getNumberOfRouters()] *
             threadZero->getQualityParams().L * 1000) /
            (SPEED_OF_LIGHT / threadZero->getQualityParams().refractive_index);
      } else {
        stats.totalSetupDelay += getGlobalTime() - ccce->requestBeginTime;
      }

#ifndef NO_ALLEGRO
      indx = ccce->destinationRouterIndex;  // this code keeps track of how many
                                            // times each router is the
                                            // destination for a connection.
      routers[indx]->incConnSuccessesTo();
      indx = ccce->sourceRouterIndex;
      routers[indx]->incConnSuccessesFrom();  // i think this stuff works.
#endif

      // Yeah, we have made it back to the destination and everything worked
      // great.
      Event* event = new Event();
      DestroyConnectionProbeEvent* dcpe = new DestroyConnectionProbeEvent();

      event->e_type = DESTROY_CONNECTION_PROBE;
      event->e_time = getGlobalTime() + ccce->connectionDuration;
      event->e_data = dcpe;

      dcpe->connectionLength = ccce->connectionLength;
      dcpe->connectionPath = ccce->connectionPath;
      dcpe->numberOfHops = 0;
      dcpe->wavelength = ccce->wavelength;
      dcpe->session = ccce->session;
      dcpe->sequence = ccce->sequence;
      dcpe->probes = ccce->probes;

      queue->addEvent(*event);

      delete event;

      if (threadZero->getQualityParams().q_factor_stats == true ||
          CurrentRoutingAlgorithm == Q_MEASUREMENT ||
          CurrentRoutingAlgorithm == ADAPTIVE_QoS) {
        EstablishedConnection* ec = new EstablishedConnection();

        ec->connectionLength = ccce->connectionLength;
        ec->connectionPath = ccce->connectionPath;
        ec->wavelength = ccce->wavelength;
        ec->connectionStartTime = getGlobalTime();
        ec->connectionEndTime = getGlobalTime() + ccce->connectionDuration;

        if (threadZero->getQualityParams().q_factor_stats == true) {
          ec->QFactors = new std::vector<double>;
          ec->QTimes = new std::vector<double>;
        }

        for (size_t p = 0; p < ec->connectionLength; ++p) {
          ec->connectionPath[p]->insertEstablishedConnection(ec);
        }
      }

      stats.totalHopCount += ccce->connectionLength;

      for (size_t p = 0; p < ccce->connectionLength; ++p) {
        stats.totalSpanCount += ccce->connectionPath[p]->getNumberOfSpans();
      }

      if (CurrentRoutingAlgorithm == Q_MEASUREMENT ||
          CurrentRoutingAlgorithm == ADAPTIVE_QoS) {
        updateQMDegredation(ccce->connectionPath, ccce->connectionLength,
                            ccce->wavelength);
      } else if (threadZero->getQualityParams().q_factor_stats == true) {
        updateQFactorStats(ccce->connectionPath, ccce->connectionLength,
                           ccce->wavelength);
      }

      if (CurrentProbeStyle == SERIAL &&
          CurrentRoutingAlgorithm != SHORTEST_PATH) {
        delete[] ccce->kPaths->pathcost;
        delete[] ccce->kPaths->pathinfo;
        delete[] ccce->kPaths->pathlen;

        delete ccce->kPaths;
      }

      delete ccce;
    } else if (CurrentProbeStyle == SERIAL &&
               ccce->sequence < ccce->max_sequence - 1) {
      long long int probesToSend = 0;
      long long int probeStart = 0;
      long long int probesSkipped = ccce->sequence + 1;

      ConnectionRequestEvent* cre = new ConnectionRequestEvent;

      cre->connectionDuration = ccce->connectionDuration;
      cre->destinationRouterIndex = ccce->destinationRouterIndex;
      cre->requestBeginTime = ccce->requestBeginTime;
      cre->session = ccce->session;
      cre->sequence = ccce->sequence + 1;
      cre->sourceRouterIndex = ccce->sourceRouterIndex;
      cre->wavelength = ccce->originalWavelength;
      cre->max_sequence = ccce->max_sequence;

      if (ccce->wavelength == QUALITY_FAILURE)
        cre->qualityFail = true;
      else
        cre->qualityFail = false;

      CreateConnectionProbeEvent** probesList = calcProbesToSend(
          cre, ccce->kPaths, probesToSend, probeStart, probesSkipped);

      sendProbes(cre, ccce->kPaths, probesList, probesToSend, probeStart,
                 probesSkipped);

      delete[] ccce->connectionPath;
      delete ccce;

      delete cre;
    } else if (ccce->wavelength == COLLISION_FAILURE) {
      if (CurrentProbeStyle == PARALLEL && ccce->finalFailure == true) {
        clearResponses(ccce->probes[ccce->sequence]);
        ++stats.CollisionFailures;
      } else if (CurrentProbeStyle != PARALLEL) {
        ++stats.CollisionFailures;

        if (CurrentProbeStyle == SERIAL &&
            CurrentRoutingAlgorithm != SHORTEST_PATH) {
          delete[] ccce->kPaths->pathcost;
          delete[] ccce->kPaths->pathinfo;
          delete[] ccce->kPaths->pathlen;

          delete ccce->kPaths;
        }

        delete[] ccce->connectionPath;
      }

      delete ccce;
    } else if (ccce->wavelength == QUALITY_FAILURE) {
      if (CurrentRoutingAlgorithm == ADAPTIVE_QoS) {
        getRouterAt(ccce->sourceRouterIndex)->incrementQualityFailures();
      }

      if (CurrentProbeStyle == PARALLEL) {
        clearResponses(ccce->probes[ccce->sequence]);
      }

      ++stats.QualityFailures;

      if (CurrentProbeStyle == SERIAL &&
          CurrentRoutingAlgorithm != SHORTEST_PATH) {
        delete[] ccce->kPaths->pathcost;
        delete[] ccce->kPaths->pathinfo;
        delete[] ccce->kPaths->pathlen;

        delete ccce->kPaths;
      }

      delete[] ccce->connectionPath;
      delete ccce;
    } else if (ccce->wavelength == NO_PATH_FAILURE) {
      if (CurrentRoutingAlgorithm == ADAPTIVE_QoS) {
        getRouterAt(ccce->sourceRouterIndex)->incrementWaveFailures();
      }

      if (CurrentProbeStyle == PARALLEL) {
        clearResponses(ccce->probes[ccce->sequence]);
      }

      ++stats.NoPathFailures;

      if (CurrentProbeStyle == SERIAL &&
          CurrentRoutingAlgorithm != SHORTEST_PATH) {
        delete[] ccce->kPaths->pathcost;
        delete[] ccce->kPaths->pathinfo;
        delete[] ccce->kPaths->pathlen;

        delete ccce->kPaths;
      }

      delete[] ccce->connectionPath;
      delete ccce;
    }
  } else {
    threadZero->recordEvent(
        "ERROR: Invalid connection confirmation, something has really gone "
        "wrong here.",
        true, controllerIndex);
    exit(ERROR_INVALID_CONFIRMATION);
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	collision_notification
// Description:		Handles the probe message for creating a new
//					connection.
//
///////////////////////////////////////////////////////////////////
void Thread::collision_notification(CollisionNotificationEvent* cne) {
  if (cne->connectionPath[cne->numberOfHops]->getSourceIndex() !=
      cne->sourceRouterIndex) {
    if (cne->connectionPath[cne->numberOfHops]->getStatus(cne->wavelength) ==
        EDGE_USED) {
      cne->connectionPath[cne->numberOfHops]->setFree(cne->wavelength);
    } else {
      threadZero->recordEvent(
          "ERROR: Something has gone wrong here...this edge should be used and "
          "it is free.",
          true, controllerIndex);
      exit(ERROR_EDGE_IS_FREE);
    }
  }

  ++cne->numberOfHops;

  if (cne->numberOfHops == cne->connectionLength) {
    if (cne->max_sequence > 1) {
      if (cne->finalFailure == false) {
        cne->probes[cne->sequence]->wavelength = COLLISION_FAILURE;

        if (otherResponse(cne->probes[cne->sequence]) != -1) {
          long long int sequence = otherResponse(cne->probes[cne->sequence]);

          --cne->probes[sequence]->numberOfHops;

          create_connection_probe(cne->probes[sequence]);
        }
      }
    }

    delete[] cne->connectionPath;
    delete cne;
  } else {
    Event* event = new Event();

    event->e_type = COLLISION_NOTIFICATION;
    event->e_time =
        getGlobalTime() +
        calculateDelay(
            cne->connectionPath[cne->numberOfHops]->getNumberOfSpans());
    event->e_data = cne;

    queue->addEvent(*event);

    delete event;
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	setQualityParameters
// Description:		Opens up the Quality Paremeters file that was
//					specified by the command line arguement
//
///////////////////////////////////////////////////////////////////
void Thread::setQualityParameters(const std::string& f) {
  // Default setting is uniform. Can be modifed using the parameter file.
  qualityParams.dest_dist = UNIFORM;

  std::string log = "Reading Quality Parameters from " + f + " file.";
  threadZero->recordEvent(log, true, 0);

  std::ifstream inFile(f);

  if (!inFile.is_open()) {
    std::cerr << "Error opening quality file: " << f << std::endl;
    exit(ERROR_QUALITY_FILE);
  }

  std::string line;

  while (std::getline(inFile, line)) {
    std::string param;
    std::string value;

    std::vector<std::string> entire_line_tokens = split(line, '=');

    if (entire_line_tokens.size() < 2) {
      std::cerr << "Invalid algorithm line: " << line << std::endl;
      exit(ERROR_QUALITY_FILE);
    }

    std::vector<std::string> value_tokens = split(entire_line_tokens[1], ' ');

    if (value_tokens.size() < 2) {
      std::cerr << "Invalid algorithm line: " << line << std::endl;
      exit(ERROR_QUALITY_FILE);
    }

    param = entire_line_tokens[0];
    value = value_tokens[0];

    if (param == "arrival_interval") {
      qualityParams.arrival_interval = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tarrival_interval = " << qualityParams.arrival_interval;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "duration") {
      qualityParams.duration = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tduration = " << qualityParams.duration;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "nonlinear_halfwin") {
      qualityParams.nonlinear_halfwin = std::stoi(value);
      std::ostringstream buffer;
      buffer << "\tnonlinear_halfwin = " << qualityParams.nonlinear_halfwin;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "halfwavelength") {
      qualityParams.halfwavelength = std::stoi(value);
      std::ostringstream halfbuffer;
      halfbuffer << "\thalfwavelength = " << qualityParams.halfwavelength;
      threadZero->recordEvent(halfbuffer.str(), true, 0);

      numOfWavelengths = 2 * qualityParams.halfwavelength + 1;
      std::ostringstream buffer;
      buffer << "\tsys_wavelength = " << threadZero->getNumberOfWavelengths();
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "fc") {
      qualityParams.fc = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tfc = " << qualityParams.fc;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "f_step") {
      qualityParams.f_step = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tf_step = " << qualityParams.f_step;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "channel_power") {
      qualityParams.channel_power = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tchannel_power = " << qualityParams.channel_power;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "L") {
      qualityParams.L = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tL = " << qualityParams.L;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "alphaDB") {
      qualityParams.alphaDB = std::stod(value);
      qualityParams.alpha =
          double(qualityParams.alphaDB * 0.1 / log10(exp(1.0)));
      std::ostringstream buffer;
      buffer << "\talphaDB = " << qualityParams.alphaDB;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "D") {
      qualityParams.D = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tD = " << qualityParams.D;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "S") {
      qualityParams.S = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tS = " << qualityParams.S;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "gamma") {
      qualityParams.gamma = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tgamma = " << qualityParams.gamma;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "QFactor_factor") {
      qualityParams.QFactor_factor = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tQFactor_factor = " << qualityParams.QFactor_factor;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "EDFA_Noise_Figure") {
      qualityParams.EDFA_Noise_Figure = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tEDFA_Noise_Figure = " << qualityParams.EDFA_Noise_Figure;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "EDFA_Gain") {
      qualityParams.EDFA_Gain = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tEDFA_Gain = " << qualityParams.EDFA_Gain;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "B_w") {
      qualityParams.B_w = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tB_w = " << qualityParams.B_w;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "usage_update_interval") {
      qualityParams.usage_update_interval = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tusage_update_interval = "
             << qualityParams.usage_update_interval;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "beta") {
      qualityParams.beta = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tbeta = " << qualityParams.beta;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "gui_update_interval") {
      qualityParams.gui_update_interval = std::stoi(value);
      std::ostringstream buffer;
      buffer << "\tgui_update_interval = " << qualityParams.gui_update_interval;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "refractive_index") {
      qualityParams.refractive_index = std::stod(value);
      std::ostringstream buffer;
      buffer << "\trefractive_index = " << qualityParams.refractive_index;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "q_factor_stats") {
      if (std::stoi(value) == 1)
        qualityParams.q_factor_stats = true;
      else if (std::stoi(value) == 0)
        qualityParams.q_factor_stats = false;
      else {
        std::ostringstream buffer;
        buffer << "Unexpected value input for q_factor_stats.";
        threadZero->recordEvent(buffer.str(), true, 0);
        qualityParams.q_factor_stats = false;
      }

      std::ostringstream buffer;
      buffer << "\tq_factor_stats = " << qualityParams.q_factor_stats;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "detailed_log") {
      if (std::stoi(value) == 1)
        qualityParams.detailed_log = true;
      else if (std::stoi(value) == 0)
        qualityParams.detailed_log = false;
      else {
        std::ostringstream buffer;
        buffer << "Unexpected buffer input for detailed_log.";
        threadZero->recordEvent(buffer.str(), true, 0);
        qualityParams.detailed_log = false;
      }

      std::ostringstream buffer;
      buffer << "\tdetailed_log = " << qualityParams.detailed_log;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "dest_dist") {
      if (std::stoi(value) == 1)
        qualityParams.dest_dist = UNIFORM;
      else if (std::stoi(value) == 2)
        qualityParams.dest_dist = DISTANCE;
      else if (std::stoi(value) == 3)
        qualityParams.dest_dist = INVERSE_DISTANCE;
      else {
        std::ostringstream buffer;
        buffer << "Unexpected value input for dest_dist.";
        threadZero->recordEvent(buffer.str(), true, 0);
        qualityParams.detailed_log = false;
      }

      std::ostringstream buffer;
      buffer << "\tdest_dist = " << qualityParams.dest_dist;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "DP_alpha") {
      qualityParams.DP_alpha = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tDP_alpha = " << qualityParams.DP_alpha;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "ACO_ants") {
      qualityParams.ACO_ants = std::stoi(value);
      std::ostringstream buffer;
      buffer << "\tACO_ants = " << qualityParams.ACO_ants;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "ACO_alpha") {
      qualityParams.ACO_alpha = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tACO_alpha = " << qualityParams.ACO_alpha;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "ACO_beta") {
      qualityParams.ACO_beta = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tACO_beta = " << qualityParams.ACO_beta;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "ACO_rho") {
      qualityParams.ACO_rho = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tACO_rho = " << qualityParams.ACO_rho;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "MM_ACO_gamma") {
      qualityParams.MM_ACO_gamma = std::stod(value);
      std::ostringstream buffer;
      buffer << "\tMM_ACO_gamma = " << qualityParams.MM_ACO_gamma;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "MM_ACO_N_iter") {
      qualityParams.MM_ACO_N_iter = std::stoi(value);
      std::ostringstream buffer;
      buffer << "\tMM_ACO_N_iter = " << qualityParams.MM_ACO_N_iter;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else if (param == "MM_ACO_N_reset") {
      qualityParams.MM_ACO_N_reset = std::stoi(value);
      std::ostringstream buffer;
      buffer << "\tMM_ACO_N_reset = " << qualityParams.MM_ACO_N_reset;
      threadZero->recordEvent(buffer.str(), true, 0);
    } else {
      threadZero->recordEvent("ERROR: Invalid line in the input file!!!", true,
                              0);
      inFile.close();
      exit(ERROR_QUALITY_INPUT);
    }
  }

  qualityParams.ASE_perEDFA = new double[getNumberOfWavelengths()];

  for (size_t w = 0; w < getNumberOfWavelengths(); ++w) {
    double F_n = pow(10.0, (qualityParams.EDFA_Noise_Figure / 10.0));
    double h = 6.6260689633e-34;
    double f_c;

    if (w < qualityParams.halfwavelength) {
      f_c = qualityParams.fc -
            (qualityParams.halfwavelength - w) * qualityParams.f_step;
    } else if (w > qualityParams.halfwavelength) {
      f_c = qualityParams.fc +
            (w - qualityParams.halfwavelength) * qualityParams.f_step;
    } else {
      f_c = qualityParams.fc;
    }

    double G = pow(10.0, (qualityParams.EDFA_Gain / 10.0));
    double B_w = qualityParams.B_w;
    double P_ch = qualityParams.channel_power;

    qualityParams.ASE_perEDFA[w] =
        2.0 * P_ch * (F_n * h * f_c * (G - 1.0) * B_w);
  }

  inFile.close();
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	setTopologyParameters
// Description:		Opens up the Topology Paremeters file that was
//					specified by the command line arguement
//
///////////////////////////////////////////////////////////////////
void Thread::setTopologyParameters(const std::string& f) {
  std::ostringstream buffer;
  buffer << "Reading Topology Parameters from " << f << " file.";
  threadZero->recordEvent(buffer.str(), false, 0);

#ifndef NO_ALLEGRO
  strcpy(topoFile, f);
#endif

  std::ifstream inFile(f);

  if (!inFile.is_open()) {
    std::cerr << "Error opening toplogy file: " << f << std::endl;
    exit(ERROR_TOPOLOGY_FILE);
  }

  std::string line;

  numberOfRouters = 0;
  numberOfEdges = 0;

  while (std::getline(inFile, line)) {
    std::vector<std::string> tokens = split(line, '=');

    if (tokens.size() != 2) {
      std::string err = "Invalid line in topology file: " + line;
      threadZero->recordEvent(err, false, 0);
      inFile.close();
      exit(ERROR_TOPOLOGY_FILE);
    }

    std::string param = tokens[0];

    if (param == "Router") {
      Router* r = new Router;

      r->setIndex(numberOfRouters);

#ifndef NO_ALLEGRO
      std::vector<std::string> coordinates = split(tokens[1], ',');

      if (coordinates.size() != 2) {
        std::string err = "Invalid router line in topology file: " + line;
        threadZero->recordEvent(err, false, 0);
        inFile.close();
        exit(ERROR_TOPOLOGY_INPUT_ROUTERS)
      }

      r->setXPercent(std::stoi(coordinates[0]));
      r->setYPercent(std::stoi(coordinates[1]))
          :
#endif
            addRouter(r);

      ++numberOfRouters;
    } else if (param == "Edge") {
      std::vector<std::string> coordinates = split(tokens[1], ',');

      if (coordinates.size() != 3) {
        std::string err = "Invalid edge line in topology file: " + line;
        threadZero->recordEvent(err, false, 0);
        inFile.close();
        exit(ERROR_TOPOLOGY_INPUT_EDGES);
      }

      size_t from = std::stoi(coordinates[0]);
      size_t to = std::stoi(coordinates[1]);
      size_t spans = std::stoi(coordinates[2]);

      Edge* e1 = new Edge(from, to, spans);
      Edge* e2 = new Edge(to, from, spans);

      getRouterAt(from)->addEdge(e1);
      getRouterAt(to)->addEdge(e2);

      numberOfEdges += 2;
    } else {
      threadZero->recordEvent("ERROR: Invalid line in the input file!!!", true,
                              0);
      inFile.close();
      exit(ERROR_TOPOLOGY_FILE);
    }
  }

  std::ostringstream routers;
  routers << "\tCreated " << numberOfRouters << " routers.";
  threadZero->recordEvent(buffer.str(), false, 0);

  std::ostringstream edges;
  edges << "\tCreated " << numberOfEdges << " edges.";
  threadZero->recordEvent(buffer.str(), false, 0);

  inFile.close();
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	setWorkstationParameters
// Description:		Opens up the Workstation Paremeters file that was
//					specified by the command line arguement
//
///////////////////////////////////////////////////////////////////
void Thread::setWorkstationParameters(const std::string& f) {
  std::ostringstream reading;
  reading << "Reading Workstation Parameters from " << f << " file.";
  threadZero->recordEvent(reading.str(), false, 0);

#ifndef NO_ALLEGRO
  strcpy(wkstFile, f);
#endif

  try {
    std::ifstream inFile(f);

    if (!inFile.is_open()) {
      std::cerr << "Error opening workstation file: " << f << std::endl;
      exit(ERROR_WORKSTATION_FILE);
    }

    std::string line;

    if (!std::getline(inFile, line)) {
      std::cerr << "Error reading workstation file: " << f << std::endl;
      exit(ERROR_WORKSTATION_FILE);
    }

    std::vector<std::string> tokens = split(line, '=');

    if (tokens.size() != 2) {
      std::cerr << "Invalid algorithm line: " << line;
      exit(ERROR_WORKSTATION_FILE);
    }

    if (tokens[0] == "NumberOfWorkstations")
      numberOfWorkstations = std::stoi(tokens[1]);
    else {
      std::cerr << "ERROR: Invalid line in the input file!!!" << std::endl;
      exit(ERROR_WORKSTATION_INPUT_QUANTITY);
    }

    // If workstations are not specified, then just uniformly distribute them
    // amongst the routers.
    for (size_t n = workstations.size(); n < numberOfWorkstations; ++n) {
	  addWorkstation(new Workstation(n % getNumberOfRouters()));
    }

    inFile.close();
  } catch (...) {
    std::cerr << "Error reading/opening workstation file: " << f << std::endl;
    exit(ERROR_WORKSTATION_FILE);
  }

  std::ostringstream created;
  created << "\tCreated " << numberOfWorkstations << " workstations.";
  threadZero->recordEvent(created.str(), false, 0);
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	setAlgorithmParameters
// Description:		Opens up the Quality Parameters file that was
//					specified by the command line argument
//
///////////////////////////////////////////////////////////////////
void Thread::setAlgorithmParameters(const std::string& f, size_t iterationCount) {
  RoutingAlgorithmNames[SHORTEST_PATH] = std::string("SP");
  RoutingAlgorithmNames[PABR] = std::string("PABR");
  RoutingAlgorithmNames[LORA] = std::string("LORA");
  RoutingAlgorithmNames[IMPAIRMENT_AWARE] = std::string("IA");
  RoutingAlgorithmNames[Q_MEASUREMENT] = std::string("QM");
  RoutingAlgorithmNames[ADAPTIVE_QoS] = std::string("AQoS");
  RoutingAlgorithmNames[DYNAMIC_PROGRAMMING] = std::string("DP");
  RoutingAlgorithmNames[ACO] = std::string("ACO");
  RoutingAlgorithmNames[MAX_MIN_ACO] = std::string("MM-ACO");

  WavelengthAlgorithmNames[FIRST_FIT] = std::string("FF");
  WavelengthAlgorithmNames[FIRST_FIT_ORDERED] = std::string("FFwO");
  WavelengthAlgorithmNames[BEST_FIT] = std::string("BF");
  WavelengthAlgorithmNames[RANDOM_FIT] = std::string("RP");
  WavelengthAlgorithmNames[QUAL_FIRST_FIT] = std::string("Q-FF");
  WavelengthAlgorithmNames[QUAL_FIRST_FIT_ORDERED] = std::string("Q-FFwO");
  WavelengthAlgorithmNames[QUAL_RANDOM_FIT] = std::string("Q-RP");
  WavelengthAlgorithmNames[LEAST_QUALITY] = std::string("LQ");
  WavelengthAlgorithmNames[MOST_QUALITY] = std::string("MQ");
  WavelengthAlgorithmNames[MOST_USED] = std::string("MU");
  WavelengthAlgorithmNames[QUAL_MOST_USED] = std::string("Q-MU");

  ProbeStyleNames[SINGLE] = std::string("SINGLE");
  ProbeStyleNames[SERIAL] = std::string("SERIAL");
  ProbeStyleNames[PARALLEL] = std::string("PARALLEL");

  std::ostringstream log;
  log << "Reading Algorithm Parameters from " << f << " file.";
  threadZero->recordEvent(log.str(), false, 0);

  std::string line;

  std::ifstream inFile(f.c_str());

  while (std::getline(inFile, line)) {
    std::vector<std::string> tokens = split(line, ',');

    if (tokens.size() != 5) {
      std::string token = "Invalid algorithm line: " + line;
      threadZero->recordEvent(token, true, 0);
      continue;
    }

    std::string ra = tokens[0].replace(0, 3, "");
    std::string wa = tokens[1].replace(0, 3, "");
    std::string ps = tokens[2].replace(0, 3, "");
    std::string qa = tokens[3].replace(0, 3, "");
    std::string run = tokens[4].replace(0, 4, "");

    int i_qa = std::stoi(qa);
    int i_run = std::stoi(run);

    if (i_run) {
      CurrentRoutingAlgorithm = NUMBER_OF_ROUTING_ALGORITHMS;
      CurrentWavelengthAlgorithm = NUMBER_OF_WAVELENGTH_ALGORITHMS;
      CurrentProbeStyle = NUMBER_OF_PROBE_STYLES;
      CurrentQualityAware = true;

      for (size_t r = 0; r < NUMBER_OF_ROUTING_ALGORITHMS; ++r) {
        if (RoutingAlgorithmNames[r].compare(ra) == 0) {
          CurrentRoutingAlgorithm = static_cast<RoutingAlgorithm>(r);
          break;
        }
      }

      for (size_t w = 0; w < NUMBER_OF_WAVELENGTH_ALGORITHMS; ++w) {
        if (WavelengthAlgorithmNames[w].compare(wa) == 0) {
          CurrentWavelengthAlgorithm = static_cast<WavelengthAlgorithm>(w);
          break;
        }
      }

      for (size_t p = 0; p < NUMBER_OF_PROBE_STYLES; ++p) {
        if (ProbeStyleNames[p].compare(ps) == 0) {
          CurrentProbeStyle = static_cast<ProbeStyle>(p);
          break;
        }
      }

      if (i_qa) {
        CurrentQualityAware = true;
      } else {
        CurrentQualityAware = false;
      }

      if (CurrentRoutingAlgorithm == NUMBER_OF_ROUTING_ALGORITHMS ||
          CurrentWavelengthAlgorithm == NUMBER_OF_WAVELENGTH_ALGORITHMS ||
          CurrentProbeStyle == NUMBER_OF_PROBE_STYLES) {
        threadZero->recordEvent("ERROR: Invalid line in the input file!!!",
                                true, 0);
        exit(ERROR_ALGORITHM_INPUT);
      } else {
        size_t iterationWorkstationDelta = static_cast<size_t>(
            double(1.0) / double(iterationCount) *
            double(threadZero->getNumberOfWorkstations()));

        for (size_t i = 0; i < iterationCount; ++i) {
          AlgorithmToRun* ap = new AlgorithmToRun;

          ap->ra = CurrentRoutingAlgorithm;
          ap->wa = CurrentWavelengthAlgorithm;
          ap->ps = CurrentProbeStyle;
          ap->qa = CurrentQualityAware;
          ap->workstations = (i + 1) * iterationWorkstationDelta;

          algParams.push_back(ap);
        }
      }
    }
  }

  inFile.close();
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	update_link_usage
// Description:		Updates the link usage to be used for the
//					LORA and PABR algorithms
//
///////////////////////////////////////////////////////////////////
void Thread::update_link_usage() {
  time_t start;
  time_t end;

  time(&start);

  if (getGlobalTime() >= TEN_HOURS) {
    return;
  }

  for (size_t r = 0; r < getNumberOfRouters(); ++r) {
    getRouterAt(r)->updateUsage();
  }

  Event* event = new Event();

  event->e_type = UPDATE_USAGE;
  event->e_time =
      getGlobalTime() + threadZero->getQualityParams().usage_update_interval;
  event->e_data = nullptr;

  queue->addEvent(*event);

  delete event;

  time(&end);
  stats.raRunTime += difftime(end, start);
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	getKShortestPaths
// Description:		Calculates the k shortest paths based upon
//					the current routing algorithm.
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* Thread::getKShortestPaths(ConnectionRequestEvent* cre,
                                               size_t probesToSend) {
  kShortestPathReturn* kPath;

  time_t start;
  time_t end;

  time(&start);

  if (CurrentRoutingAlgorithm == PABR) {
    kPath = threadZero->getResourceManager()->calculate_PAR_path(
        cre->sourceRouterIndex, cre->destinationRouterIndex, probesToSend,
        controllerIndex);
  } else if (CurrentRoutingAlgorithm == SHORTEST_PATH) {
    kPath = threadZero->getResourceManager()->calculate_SP_path(
        cre->sourceRouterIndex, cre->destinationRouterIndex, probesToSend,
        controllerIndex);
  } else if (CurrentRoutingAlgorithm == LORA) {
    kPath = threadZero->getResourceManager()->calculate_LORA_path(
        cre->sourceRouterIndex, cre->destinationRouterIndex, probesToSend,
        controllerIndex);
  } else if (CurrentRoutingAlgorithm == IMPAIRMENT_AWARE) {
    kPath = threadZero->getResourceManager()->calculate_IA_path(
        cre->sourceRouterIndex, cre->destinationRouterIndex, controllerIndex);
  } else if (CurrentRoutingAlgorithm == Q_MEASUREMENT) {
    kPath = threadZero->getResourceManager()->calculate_QM_path(
        cre->sourceRouterIndex, cre->destinationRouterIndex, probesToSend,
        controllerIndex);
  } else if (CurrentRoutingAlgorithm == ADAPTIVE_QoS) {
    if (getRouterAt(cre->sourceRouterIndex)->getQualityFailures() >=
        getRouterAt(cre->sourceRouterIndex)->getWaveFailures()) {
      kPath = threadZero->getResourceManager()->calculate_QM_path(
          cre->sourceRouterIndex, cre->destinationRouterIndex, probesToSend,
          controllerIndex);
    } else {
      kPath = threadZero->getResourceManager()->calculate_AQoS_path(
          cre->sourceRouterIndex, cre->destinationRouterIndex, probesToSend,
          controllerIndex);
    }
  } else if (CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING) {
    kPath = threadZero->getResourceManager()->calculate_DP_path(
        cre->sourceRouterIndex, cre->destinationRouterIndex, probesToSend,
        controllerIndex);
  } else if (CurrentRoutingAlgorithm == ACO) {
    kPath = threadZero->getResourceManager()->calculate_ACO_path(
        cre->sourceRouterIndex, cre->destinationRouterIndex, probesToSend,
        controllerIndex);
  } else if (CurrentRoutingAlgorithm == MAX_MIN_ACO) {
    kPath = threadZero->getResourceManager()->calculate_MM_ACO_path(
        cre->sourceRouterIndex, cre->destinationRouterIndex, probesToSend,
        controllerIndex);
  } else {
    std::ostringstream buffer;
    buffer << "Invalid values for CurrentAlgorithm (" << CurrentRoutingAlgorithm
           << ") and CurrentWavelengthAlgorithm (" << CurrentWavelengthAlgorithm
           << ")";
    threadZero->recordEvent(buffer.str(), true, controllerIndex);
    exit(ERROR_ALGORITHM_INPUT);
  }

  time(&end);

  stats.raRunTime += difftime(end, start);

  return kPath;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calcProbesToSend
// Description:		Calculates which probes to send based upon the
//					kPaths and the Routing/Wavelength
//Algorithms
//
///////////////////////////////////////////////////////////////////
CreateConnectionProbeEvent** Thread::calcProbesToSend(
    ConnectionRequestEvent* cre, kShortestPathReturn* kPath,
    long long int& probesToSend, long long int& probeStart,
    long long int& probesSkipped) {
  CreateConnectionProbeEvent** probesList = nullptr;

  // We have to handle some algorithms differently because it uses a forward
  // reservation scheme, while the others use a backward reservation scheme.
  if (CurrentRoutingAlgorithm == IMPAIRMENT_AWARE ||
      CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING) {
    long long int probesTotal = 0;
    long long int probeFirst = -1;

    if (CurrentRoutingAlgorithm == IMPAIRMENT_AWARE) {
      for (long long int w = probesSkipped;
           w < static_cast<long long int>(threadZero->getNumberOfWavelengths());
           ++w) {
        if (kPath->pathcost[w] != std::numeric_limits<double>::infinity()) {
          ++probesTotal;

          if (probeFirst == -1) probeFirst = w;
        }
      }
    } else if (CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING) {
      size_t max_probes = threadZero->getQualityParams().max_probes;

      if (CurrentProbeStyle == SINGLE) max_probes = 1;

      for (size_t w = probesSkipped; w < max_probes; ++w) {
        if (kPath->pathcost[w] != std::numeric_limits<double>::infinity()) {
          ++probesTotal;

          if (probeFirst == -1) probeFirst = w;
        }
      }
    }

    if (probeFirst == -1) {
      probesToSend = 0;
    } else {
      if (CurrentProbeStyle == PARALLEL) {
        probesToSend = probesTotal;

        probesList = new CreateConnectionProbeEvent*[probesToSend];
      } else if (CurrentProbeStyle == SINGLE || CurrentProbeStyle == SERIAL) {
        probesToSend = 1;

        if (CurrentWavelengthAlgorithm == FIRST_FIT ||
            CurrentWavelengthAlgorithm == QUAL_FIRST_FIT) {
          probeStart = probeFirst;
        } else if (CurrentWavelengthAlgorithm == BEST_FIT) {
          time_t start;
          time_t end;

          time(&start);

          double minCost = std::numeric_limits<double>::infinity();

          for (size_t w = 0; w < threadZero->getNumberOfWavelengths();
               ++w) {
            if (kPath->pathcost[w] < minCost) {
              minCost = kPath->pathcost[w];
              probeStart = w;
            }
          }

          if (minCost == std::numeric_limits<double>::infinity())
            probesToSend = 0;

          time(&end);

          stats.raRunTime += difftime(end, start);
        } else {
          threadZero->recordEvent(
              "ERROR: Invalid Wavelength Algorithm for Impairment Aware.", true,
              controllerIndex);
          exit(ERROR_WAVELENGTH_ALGORITHM_IA);
        }
      }
    }
  } else if (CurrentProbeStyle == SINGLE) {
    probeStart = 0;

    if (kPath->pathcost[0] != std::numeric_limits<double>::infinity()) {
      probesToSend = 1;
    } else {
      probesToSend = 0;
    }
  } else if (CurrentProbeStyle == SERIAL) {
    probesToSend = 0;

    for (size_t a = static_cast<size_t>(probesSkipped);
         a < threadZero->getQualityParams().max_probes; ++a) {
      if (kPath->pathcost[a] != std::numeric_limits<double>::infinity()) {
        probeStart = a;
        probesToSend = 1;
        break;
      }
    }
  } else if (CurrentProbeStyle == PARALLEL) {
    probeStart = 0;
    probesToSend = 0;

    for (size_t p = 0; p < threadZero->getQualityParams().max_probes;
         ++p) {
      if (kPath->pathcost[p] != std::numeric_limits<double>::infinity())
        ++probesToSend;
    }

    if (probesToSend > 0) {
      probesList = new CreateConnectionProbeEvent*[probesToSend];
    }
  }

  return probesList;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	sendProbes
// Description:		Generates the probe events based upon the
//					input parameters.
//
///////////////////////////////////////////////////////////////////
void Thread::sendProbes(ConnectionRequestEvent* cre, kShortestPathReturn* kPath,
                        CreateConnectionProbeEvent** probesList,
                        long long int probesToSend, long long int probeStart,
                        long long int probesSkipped) {
  if (probesToSend == 0) {
    if (CurrentProbeStyle == SERIAL) {
      if (cre->wavelength >= 0 &&
          CurrentRoutingAlgorithm != DYNAMIC_PROGRAMMING) {
        ++stats.QualityFailures;
      } else {
        ++stats.NoPathFailures;
      }
    } else {
      ++stats.NoPathFailures;
    }
  } else {
    long long int newProbesSkipped = 0;

    for (long long int p = probeStart;
         p < probeStart + newProbesSkipped + probesToSend; ++p) {
      if (CurrentRoutingAlgorithm == IMPAIRMENT_AWARE ||
          CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING) {
        while (kPath->pathcost[p] == std::numeric_limits<double>::infinity()) {
          ++newProbesSkipped;
          ++p;
        }
      }

      CreateConnectionProbeEvent* probe = new CreateConnectionProbeEvent();

      probe->sourceRouterIndex = cre->sourceRouterIndex;
      probe->destinationRouterIndex = cre->destinationRouterIndex;
      probe->connectionDuration = cre->connectionDuration;
      probe->requestBeginTime = cre->requestBeginTime;
      probe->numberOfHops = 0;
      probe->connectionLength = kPath->pathlen[p] - 1;
      probe->connectionPath = new Edge*[probe->connectionLength];
      probe->session = cre->session;
      probe->qualityFail = cre->qualityFail;

      if (CurrentProbeStyle == SERIAL)
        probe->sequence = cre->sequence;
      else if (CurrentProbeStyle == PARALLEL)
        probe->sequence = cre->sequence++;
      else if (CurrentProbeStyle == SINGLE)
        probe->sequence = 0;

      probe->atDestination = false;
      probe->decisionTime = 0.0;
      probe->kPaths = kPath;

      if (CurrentProbeStyle == PARALLEL) {
        probe->probes = probesList;
        probesList[probe->sequence] = probe;
      }

      if (CurrentProbeStyle == SERIAL) {
        if (cre->max_sequence == 0) {
          if (CurrentRoutingAlgorithm != IMPAIRMENT_AWARE &&
              CurrentRoutingAlgorithm != DYNAMIC_PROGRAMMING) {
            for (size_t a = 0;
                 a < threadZero->getQualityParams().max_probes; ++a) {
              if (kPath->pathcost[a] != std::numeric_limits<double>::infinity())
                ++cre->max_sequence;
            }
          } else if (CurrentRoutingAlgorithm == IMPAIRMENT_AWARE) {
            for (size_t a = 0; a < threadZero->getNumberOfWavelengths();
                 ++a) {
              if (kPath->pathcost[a] != std::numeric_limits<double>::infinity())
                ++cre->max_sequence;
            }
          } else if (CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING) {
            for (size_t a = 0;
                 a < threadZero->getQualityParams().max_probes; ++a) {
              if (kPath->pathcost[a] != std::numeric_limits<double>::infinity())
                ++cre->max_sequence;
            }
          }
        }

        probe->max_sequence = cre->max_sequence;
      } else {
        probe->max_sequence = probesToSend;
      }

      if (CurrentRoutingAlgorithm == IMPAIRMENT_AWARE) {
        probe->wavelength = p;
      } else if (CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING) {
        probe->wavelength = static_cast<short>(kPath->pathcost[p]);
      } else {
        probe->wavelength = 0;
      }

      for (size_t r = 0; r < probe->connectionLength; ++r) {
        probe->connectionPath[r] =
            getRouterAt(
                kPath->pathinfo[p * (threadZero->getNumberOfRouters() - 1) + r])
                ->getEdgeByDestination(
                    kPath->pathinfo[p * (threadZero->getNumberOfRouters() - 1) +
                                    r + 1]);
      }

      if (CurrentRoutingAlgorithm != SHORTEST_PATH) {
        kPath->pathcost[p] = std::numeric_limits<double>::infinity();
        kPath->pathlen[p] = std::numeric_limits<int>::infinity();
      }

      Event* event = new Event();

      event->e_type = CREATE_CONNECTION_PROBE;
      event->e_time =
          getGlobalTime() +
          calculateDelay(probe->connectionPath[0]->getNumberOfSpans());
      event->e_data = probe;

      queue->addEvent(*event);

      delete event;

      if (CurrentRoutingAlgorithm != IMPAIRMENT_AWARE &&
          CurrentRoutingAlgorithm != DYNAMIC_PROGRAMMING) {
        ++stats.ProbeSentCount;
      }
    }
  }

  if ((CurrentProbeStyle != SERIAL &&
       CurrentRoutingAlgorithm != SHORTEST_PATH) ||
      probesToSend == 0) {
    delete[] kPath->pathcost;
    delete[] kPath->pathinfo;
    delete[] kPath->pathlen;

    delete kPath;
  }
}

#ifndef NO_ALLEGRO
///////////////////////////////////////////////////////////////////
//
// Function Name:	update_gui
// Description:		Updates the gui with the most recent usage levels
//
///////////////////////////////////////////////////////////////////
void Thread::update_gui() {
  if (getGlobalTime() >= TEN_HOURS) {
    return;
  }

  for (size_t r = 0; r < getNumberOfRouters(); ++r) {
    getRouterAt(r)->updateGUI();
  }

  Event* event = new Event();

  event->e_type = UPDATE_GUI;
  event->e_time =
      getGlobalTime() + threadZero->getQualityParams().gui_update_interval;
  event->e_data = nullptr;

  queue->addEvent(*event);

  delete event;
}
#endif

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculateDelay
// Description:		Calculates the propogation delay on the optical
//					link.
//
///////////////////////////////////////////////////////////////////
double Thread::calculateDelay(size_t spans) const {
  if (CurrentRoutingAlgorithm == IMPAIRMENT_AWARE ||
      CurrentRoutingAlgorithm == DYNAMIC_PROGRAMMING) {
    return 0.0;
  } else {
    return double(spans * threadZero->getQualityParams().L * 1000) /
           (SPEED_OF_LIGHT / threadZero->getQualityParams().refractive_index);
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	sendResponse
// Description:		Determines if the router needs to respond to
//					this probe. (i.e. Either no response has been
//sent 					yet or the sent responses were declined.)
//
///////////////////////////////////////////////////////////////////
bool Thread::sendResponse(CreateConnectionProbeEvent* probe) {
  for (long long int p = 0; p < probe->max_sequence; ++p) {
    if (p != probe->sequence) {
      if (probe->probes[p]->atDestination == true) {
        if (probe->probes[p]->decisionTime != 0 &&
            probe->probes[p]->wavelength >= 0) {
          return false;
        }
      }
    }
  }

  return true;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	clearResponses
// Description:		Clears the pending responses as the connection has
//					been setup succesfully.
//
///////////////////////////////////////////////////////////////////
void Thread::clearResponses(CreateConnectionProbeEvent* probe) {
  for (long long int p = 0; p < probe->max_sequence; ++p) {
    if (p != probe->sequence) {
      delete[] probe->probes[p]->connectionPath;
      delete probe->probes[p];
    }
  }

  delete[] probe->probes;
  delete probe;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	moreProbes
// Description:		Returns true if more probes will be arriving
//					at the destination.
//
///////////////////////////////////////////////////////////////////
bool Thread::moreProbes(CreateConnectionProbeEvent* probe) {
  if (CurrentProbeStyle != PARALLEL) return false;

  for (long long int p = 0; p < probe->max_sequence; ++p) {
    if (p != probe->sequence) {
      if (probe->probes[p]->atDestination == false) {
        return true;
      }
    }
  }

  return false;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	setMinDuration
// Description:		Sets the minimum duration of a connection to
//					three times the maximum minumum distance
//between 					two difference routers.
//
///////////////////////////////////////////////////////////////////
void Thread::setMinDuration(size_t spans) {
  // Routing algorithm should be a non centralized algorithm to ensure
  // the delay is not zero.
  RoutingAlgorithm origAlg = CurrentRoutingAlgorithm;
  CurrentRoutingAlgorithm = SHORTEST_PATH;

  minDuration = 3.0 * calculateDelay(spans);

  CurrentRoutingAlgorithm = origAlg;

  maxSpans = threadZero->getMaxSpans();
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	setQFactorMin
// Description:		Sets the value of the Q-factor threshold to
//					within a specified factor of the max min
//
///////////////////////////////////////////////////////////////////
void Thread::setQFactorMin(size_t spans) {
  qualityParams.TH_Q = static_cast<double>(
      qualityParams.QFactor_factor * 10.0 *
      log10(qualityParams.channel_power /
            sqrt(qualityParams.ASE_perEDFA[qualityParams.halfwavelength] *
                 double(spans))));

  double Q = qualityParams.TH_Q;
  maxSpans = spans;

  while (Q >= qualityParams.TH_Q) {
    ++maxSpans;

    Q = 10.0 *
        log10(qualityParams.channel_power /
              sqrt(qualityParams.ASE_perEDFA[qualityParams.halfwavelength] *
                   double(maxSpans + 1)));
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	otherResponse
// Description:		Determines if another probe stored could be used
//					to responsd to this probe.
//
///////////////////////////////////////////////////////////////////
long long int Thread::otherResponse(CreateConnectionProbeEvent* probe) {
  if (CurrentProbeStyle != PARALLEL) return -1;

  for (long long int p = 0; p < probe->max_sequence; ++p) {
    if (p != probe->sequence) {
      if (probe->probes[p]->atDestination == true) {
        if (probe->probes[p]->decisionTime == 0.0) {
          return p;
        }
      }
    }
  }

  return -1;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	updateQMDegredation
// Description:		Updates the link costs when a connection is either added
//					or dropped.
//
///////////////////////////////////////////////////////////////////
void Thread::updateQMDegredation(Edge** connectionPath, size_t connectionLength,
                                 long long int wavelength) {
  time_t start;
  time_t end;

  time(&start);

  for (size_t p = 0; p < connectionLength; ++p) {
    connectionPath[p]->updateQMDegredation(controllerIndex, wavelength);
  }

  time(&end);
  stats.raRunTime += difftime(end, start);
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	updateQFactorStats
// Description:		Updates the q-factor stats when a connection is either
// added
//					or dropped.
//
///////////////////////////////////////////////////////////////////
void Thread::updateQFactorStats(Edge** connectionPath, size_t connectionLength,
                                long long int wavelength) {
  for (size_t p = 0; p < connectionLength; ++p) {
    connectionPath[p]->updateQFactorStats(controllerIndex, wavelength);
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	split
// Description:		Tokenizes a string by specified delimiter
//
///////////////////////////////////////////////////////////////////
std::vector<std::string> Thread::split(const std::string& s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

#ifndef NO_ALLEGRO
///////////////////////////////////////////////////////////////////
//
// Function Name:	detailScreen()
// Description:		displays configuration info as well as overall
//					results for the given simulation.
//
///////////////////////////////////////////////////////////////////
void Thread::detailScreen() {
  int startX = ceil((double)SCREEN_W / 2.0 - (double)detailinfo->w / 2.0);
  int startY = 100;
  int numtopaint;
  int black = makecol(0, 0, 0);

  masked_blit(detailinfo, popup, 0, 0, startX, startY, SCREEN_W, SCREEN_H);
  textprintf_right_ex(popup, font, startX + 270, startY + 35, black, -1, "%s",
                      topology);
  textprintf_right_ex(popup, font, startX + 270, startY + 52, black, -1, "%s",
                      routing);
  textprintf_right_ex(popup, font, startX + 270, startY + 69, black, -1, "%s",
                      wavelength);
  textprintf_right_ex(popup, font, startX + 270, startY + 86, black, -1, "%s",
                      probing);
  if (qualityAware)
    textprintf_right_ex(popup, font, startX + 270, startY + 103, black, -1,
                        "Yes");
  else
    textprintf_right_ex(popup, font, startX + 270, startY + 103, black, -1,
                        "No");
  textprintf_right_ex(popup, font, startX + 270, startY + 130, black, -1, "%d",
                      CurrentActiveWorkstations);
  textprintf_right_ex(popup, font, startX + 270, startY + 147, black, -1, "%d",
                      numOfWavelengths);

  textprintf_right_ex(popup, font, startX + 270, startY + 202, black, -1,
                      "%2.3f", overallBlocking);
  textprintf_right_ex(popup, font, startX + 270, startY + 219, black, -1,
                      "%2.3f", collisions);
  textprintf_right_ex(popup, font, startX + 270, startY + 238, black, -1,
                      "%2.3f", badquality);
  textprintf_right_ex(popup, font, startX + 270, startY + 255, black, -1,
                      "%2.3f", nonresource);
  textprintf_right_ex(popup, font, startX + 270, startY + 285, black, -1,
                      "%2.3f", avgprobesperrequest);
  textprintf_right_ex(popup, font, startX + 270, startY + 303, black, -1,
                      "%2.3f", avgrequestdelaytime);
  textprintf_right_ex(popup, font, startX + 270, startY + 321, black, -1,
                      "%2.3f", avgconnhopcount);
  textprintf_right_ex(popup, font, startX + 270, startY + 339, black, -1,
                      "%2.3f", avgconnspancount);
  textprintf_right_ex(popup, font, startX + 270, startY + 357, black, -1,
                      "%2.3e", avgASEnoise);
  textprintf_right_ex(popup, font, startX + 270, startY + 375, black, -1,
                      "%2.3e", avgFWMnoise);
  textprintf_right_ex(popup, font, startX + 270, startY + 393, black, -1,
                      "%2.3e", avgXPMnoise);

  while (!(mouse_b & 1)) {
  }
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	saveThread
// Description:		Creates appropriately-named .bin file and stores
//					all data needed to explore history of thread
//later.
//
///////////////////////////////////////////////////////////////////
void Thread::saveThread(char* dir) {
  char allName[75];
  char fileName[40];
  std::ofstream myFile;

  // routing algorithm - wavelength algorithm - probe style - num workstations -
  // quality aware (1 if true) - maxprobes
  sprintf(fileName, "%s-%s-%s-%d-%d-%d",
          threadZero->getRoutingAlgorithmName(CurrentRoutingAlgorithm)->c_str(),
          threadZero->getWavelengthAlgorithmName(CurrentWavelengthAlgorithm)
              ->c_str(),
          threadZero->getProbeStyleName(CurrentProbeStyle)->c_str(),
          getCurrentActiveWorkstations(), (int)getCurrentQualityAware(),
          threadZero->getQualityParams().max_probes);

  sprintf(allName, "%s/index.idx", dir);
  myFile.open(allName, std::ios::out | std::ios::app | std::ios::binary);
  myFile << fileName << "|\n";
  myFile.close();

  sprintf(allName, "%s/%s.thd", dir, fileName);

  myFile.open(allName, std::ios::out | std::ios::app | std::ios::binary);
  myFile << topology << "|\n";
  myFile << routing << "|\n";
  myFile << wavelength << "|\n";
  myFile << probing << "|\n";            // need
  myFile << (int)qualityAware << "|\n";  // need

  myFile << overallBlocking << "|\n";
  myFile << collisions << "|\n";
  myFile << badquality << "|\n";
  myFile << nonresource << "|\n";
  myFile << avgprobesperrequest << "|\n";
  myFile << avgrequestdelaytime << "|\n";
  myFile << avgconnhopcount << "|\n";
  myFile << avgconnspancount << "|\n";
  myFile << avgASEnoise << "|\n";
  myFile << avgFWMnoise << "|\n";
  myFile << avgXPMnoise << "|\n";

  myFile << threadZero->getNumberOfWavelengths() << "|\n";
  myFile << routers[0]->getEdgeByIndex(0)->getUsageNums() << "|\n";
  myFile << getCurrentActiveWorkstations() << "|\n";
  myFile.close();

  for (size_t r = 0; r < routers.size(); ++r)
    routers[r]->saveData(allName);
}
#endif

///////////////////////////////////////////////////////////////////
//
// Function Name:	recordEvent
// Description:		Passes the event to the logger for recording
//
///////////////////////////////////////////////////////////////////
void Thread::recordEvent(const std::string& s, bool print, size_t ci) {
	if (isLoadPrevious == true)
		return;
	else if (controllerIndex == 0)
		logger->recordEvent(s, print, ci);
	else
		exit(ERROR_RECORD_EVENT);
};

///////////////////////////////////////////////////////////////////
//
// Function Name:	flushLog
// Description:		Passes the flush command to the logger
//
///////////////////////////////////////////////////////////////////
void Thread::flushLog(bool print) {
	if (isLoadPrevious == true)
		return;
	else if (controllerIndex == 0)
		logger->flushLog(print);
	else
		exit(ERROR_NO_FLUSH);
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	getBeta
// Description:		Returns the beta parameter
//
///////////////////////////////////////////////////////////////////
double Thread::getBeta() const
{
  return threadZero->getQualityParams().beta;
}
