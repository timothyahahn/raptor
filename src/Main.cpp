// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Main.cpp
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the implementation of main(), whose only
//					purpose is to pass the command line arguments to the controller.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#include "ErrorCodes.h"
#include "Thread.h"

#define HAVE_STRUCT_TIMESPEC
#include "pthread.h"

#include <iostream>
#include <vector>

#ifdef RUN_GUI

#include "AllegroWrapper.h"

extern void close_button_handler();
extern void flash();
extern void explore_history();

extern void switcher( int argc, const char* argv[] );

extern int color,color2,white,explore_time;

BITMAP *buffer;
BITMAP *pointer;
BITMAP *graph;
BITMAP *graphbuttons;
BITMAP *topbuttons;
BITMAP *edgespans;
BITMAP *routersbmp;
BITMAP *routerinfo;
BITMAP *editrouterinfo;
BITMAP *editedgeinfo;
BITMAP *edgesbmp;
BITMAP *mainbuf;
BITMAP *progbarbmp;
BITMAP *routerpic;
BITMAP *graphbackground;
BITMAP *topobackground;
BITMAP *topomenu;
BITMAP *popup;
BITMAP *detailinfo;

BITMAP *colorkey;

char foldName[25];
char folder[50];

#define SCRNWID 1300
#define SCRNHEI 700
#endif

//Hate to make a global instance of the threads, but the other classes
//need access to the threads information. They can't just include a 
//pointer to the class due to a circular reference.
Thread* threadZero = nullptr;
Thread** threads = nullptr;
unsigned int threadCount = 0;

std::vector<AlgorithmToRun*> algParams;

void *runThread(void* n);
void runSimulation(int argc, const char* argv[]);

pthread_mutex_t ScheduleMutex;

int main( int argc, const char* argv[] )
{
	if(argc != 7)
	{
		std::cerr << "Usage: " << argv[0] << " <Topology> <Wavelengths> <Random Seed> <Thread Count> <Iteration Count> <Probe Count>" << std::endl;

		std::cerr << std::endl;

		std::cerr << "Topology: NSF, Mesh, Mesh6x6, Mesh8x8, Mesh10x10" << std::endl;
		std::cerr << "Wavelengths: 21, 41, 81, 161, 321, 641, 1281" << std::endl;
		std::cerr << "Random Seed: <any valid unsigned int>" << std::endl;
		std::cerr << "Thread Count: <maximum number of threads to run, 1 to n>" << std::endl;
		std::cerr << "Iteration Count: <number of iterations, 1 to n>" << std::endl;
		std::cerr << "Probe Count: <number of probes, 1 to n>" << std::endl;

		return ERROR_INVALID_PARAMETERS;
	}

#ifdef RUN_GUI
	allegro_init();
	set_close_button_callback(close_button_handler);
	set_color_depth(16);

	install_keyboard();
	install_mouse();
	install_timer();
	set_window_title("RAPTOR (Route Assignment Program for Transparent Optical Routes)");
	graph = create_bitmap(SCRNWID,SCRNHEI);
	graphbuttons = create_bitmap(SCRNWID,16);
	topbuttons = create_bitmap(SCRNWID,31);
	buffer = create_bitmap(SCRNWID,SCRNHEI);
	mainbuf = create_bitmap(SCRNWID,SCRNHEI);
	routersbmp = create_bitmap(SCRNWID,SCRNHEI);
	edgesbmp = create_bitmap(SCRNWID,SCRNHEI);
	edgespans = create_bitmap(SCRNWID,SCRNHEI);
	popup = create_bitmap(SCRNWID,SCRNHEI);
	graphbackground = load_bitmap("bitmaps/graphbackground.bmp",nullptr);
	topobackground = load_bitmap("bitmaps/topobackground2.bmp",nullptr);
	progbarbmp = load_bitmap("bitmaps/progressbar.bmp",nullptr);
	routerinfo = load_bitmap("bitmaps/routerinfo.bmp",nullptr);
	editrouterinfo = load_bitmap("bitmaps/editrouterinfo.bmp",nullptr);
	editedgeinfo = load_bitmap("bitmaps/editedgeinfo.bmp",nullptr);
	detailinfo = load_bitmap("bitmaps/configurationinfo.bmp",nullptr);
	colorkey = load_bitmap("bitmaps/colorkey.bmp",nullptr);
	pointer = load_bitmap("bitmaps/pointer2.bmp",nullptr);
	topomenu = load_bitmap("bitmaps/topomenu.bmp",nullptr);

	flash();

	set_gfx_mode(GFX_AUTODETECT_WINDOWED,SCRNWID,SCRNHEI,0,0);

	color = makecol(0,0,0);
	color2 = makecol(0,255,0);
	white = makecol(255,255,255);
	set_mouse_sprite(pointer);
	set_mouse_sprite_focus(15,15);
	show_mouse(screen);
	set_display_switch_mode(SWITCH_BACKGROUND);

	bool doMenu = true;
	//flash();
	switcher(argc, argv);

	set_mouse_sprite(nullptr);
	destroy_bitmap(buffer);
	destroy_bitmap(graph);
	destroy_bitmap(pointer);
	destroy_bitmap(topobackground);
	destroy_bitmap(popup);
	destroy_bitmap(routerinfo);
	allegro_exit();
#else
	runSimulation(argc,argv);
#endif

	return 0;
}

#ifdef RUN_GUI
END_OF_MAIN()
#endif

void runSimulation(int argc, const char* argv[])
{
	int* threadZeroReturn = 0;
	int runCount = 0;

	while(threadZeroReturn == 0 || *threadZeroReturn == Thread::MORE_SIMULATIONS)
	{
		delete threadZeroReturn;

		std::vector<pthread_t*> pThreads;

		threadCount = atoi(argv[4]);

		threads = new Thread*[threadCount];

		unsigned int iterationCount = atoi(argv[5]);

#ifdef RUN_GUI
	rectfill(screen, 0, 0, SCREEN_W, 40, color);
	textprintf_ex(screen,font,20,15,color2,color,"Building XPM Database, please wait..."); 
#endif

		Thread* thread = new Thread(0,argc,argv,false,runCount);

		threadZero->initResourceManager();

#ifdef RUN_GUI
		rectfill(screen, 0, 0, SCREEN_W, 40, color);
#endif

		pthread_mutex_init(&ScheduleMutex,nullptr);

		if(threadCount > algParams.size())
			threadCount = static_cast<unsigned int>(algParams.size());

		for(unsigned int t = 1; t < threadCount; ++t)
		{
			Thread* thread = new Thread(t,argc,argv,false,runCount);
			pThreads.push_back(new pthread_t);

			thread->initResourceManager();
		}

		std::ostringstream buffer;
		buffer << "Created " << threadCount << " threads " << std::endl;
		threadZero->recordEvent(buffer.str(),true,0);

		for(unsigned int t = 1; t < threadCount; ++t)
		{
			pthread_create(pThreads[t-1],nullptr,runThread,new unsigned int(t));
		}

		threadZeroReturn = static_cast<int*>(runThread(new unsigned int(0)));

		for(unsigned int t = 1; t < threadCount; ++t)
		{
			pthread_join(*pThreads[t-1],nullptr);

#ifdef RUN_GUI
		textprintf_ex(screen,font,20,500,makecol(0,255,0),makecol(0,0,0),"Thread %5d",t-1);
#endif
		}

		for(unsigned int t = 0; t < threadCount; ++t)
		{
			delete threads[t];

			if(t != 0)
				delete pThreads[t-1];
		}

		delete[] threads;

		pThreads.clear();

		pthread_mutex_destroy(&ScheduleMutex);

		++runCount;
	}

	delete threadZeroReturn;

#ifdef RUN_GUI

	explore_time = 0;
	while(!key[KEY_ESC])
	{
	}
#endif
}

void *runThread(void* n)
{
	unsigned int* t_id = static_cast<unsigned int *>(n);

	int *retVal = new int;

	while(algParams.size() > 0)
	{
		pthread_mutex_lock(&ScheduleMutex);

		AlgorithmToRun* alg = algParams.back();
		algParams.pop_back();

		pthread_mutex_unlock(&ScheduleMutex);

		*retVal = threads[*t_id]->runThread(alg);

#ifdef RUN_GUI
		textprintf_ex(screen,font,20,SCREEN_H-30,color2,color,"%s",folder); 
		threads[*t_id]->saveThread(folder);
#endif
	}

	delete t_id;

	return retVal;
}
