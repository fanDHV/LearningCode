#include  <cassert>
#include  <memory.h>
#include  <iostream>
#include  <fstream>
#include  <stdio.h>
#include  <stdlib.h>
#include  <string>

using namespace std;

#include  "sarray.h"
#include  "Renderer.h"
//#include  "Leda.h"
#include "myCGAL.h"
#include "global-settings.h"

Renderer g_renderer;

// man, all these wrappers are slowly killing me

void
RendererDisplayWrapper ()
{
  g_renderer.display ();
}

void
RendererReshapeWrapper (int w, int h)
{
  g_renderer.reshape (w, h);
}

void
RendererKeyboardWrapper (unsigned char key, int x, int y)
{
  g_renderer.keyboard (key, x, y);
}

void
RendererMouseWrapper (int button, int state, int x, int y)
{
  g_renderer.mouse (button, state, x, y);
}

void
RendererFlyByWrapper ()
{
  g_renderer.flyBy();
}

void
RendererUpdateDebugDrawWrapper ()
{
  g_renderer.updateDebugDraw ();
}

void
RendererAlgorithmDisplayMenuWrapper (int value)
{
  g_renderer.algorithmDisplayMenu (value);
}

void
RendererAlgorithmAnimateMenuWrapper (int value)
{
  g_renderer.algorithmAnimateMenu (value);
}

void
RendererOtherDisplayMenuWrapper (int value)
{
  g_renderer.otherDisplayMenu (value);
}

void
RendererMainMenuWrapper (int value)
{
  g_renderer.mainMenu (value);
}

void
RendererHandleSpecialKeypressWrapper (int key, int x, int y)
{
  g_renderer.handleSpecialKeypress(key, x, y);
}


int main (int argc, char **argv)
{
  bool f_test = true;         // f_test = true -> read srcdest from file.in
  bool f_super_test = false;  // f_super_test = false -> run in graphical mode ...                               
  bool f_print, f_white;

  f_white = f_print = false;

  bool f_LiuWong = false;  // f_LiuWong = true: run LiuWong algorithm

  
  //%%%%%%%%% Write into file 2022-4-6
  ofstream ofs("result.txt", ios::app);
  if(!ofs)
  {
    cerr << "Error: file not opened." << endl;
    return 0;
  }



  for (int ind = 1; ind < argc; ind++)
    {
      string arg ((char *) (argv[ind]));

      if (arg == "-h")
        {
          // use hersh-suri
          g_renderer.useHershSuri();
          g_renderer.MainLoop( string ((char *) (argv[ind + 1])), 0.0,
                               string ((char *) (argv[ind + 2])),
                               argc, argv, f_test, f_print, f_super_test,
                               f_white);

          ofs << "\n Use hersh-suri" << endl;

          ind += 2;
          continue;
        }

      if (arg == "-t")
        {
          f_test = true;
          continue;
        }
      if (arg == "-s")
        {
          f_test = true;
          f_super_test = true;
          continue;
        }
      if (arg == "-p")
        {
          f_print = true;
          continue;
        }
      if (arg == "-w")
        {
          f_white = true;
          continue;
        }
      /******************************************/
      /* Hoai - 15/3/2013                       */
      /******************************************/
      //--------------------------------------------
      // maximum number of iterations
      if ( arg == "-i" )
        {
          MAX_ITERATIONS = atoi( argv[ ind + 1 ] );
          cout << "%% max iteration: " << MAX_ITERATIONS << endl;

          ofs << " input max_iterators is" << MAX_ITERATIONS << endl;

          ind ++;
          continue;
        }
      //--------------------------------------------
      // log paths
      if ( arg == "-l" )
        {
          PATH_LOG = true;
          continue;
        }
      //--------------------------------------------
      // epsilon (currently unused)
      if ( arg == "-e" )
        {
          EPSILON = atof( argv[ ind + 1 ] );

           ofs << " input Epsilon is " << EPSILON << endl;

          ind ++;
          continue;
        }
      //--------------------------------------------
      // \epsilon_{\alpha} for stopping condition (in degree)
      if ( arg == "-a" )
        {
          // input is in degree, so conversion to radian is required
          ALPHA_EPSILON = atof( argv[ ind + 1 ] ) / 180.0 * CGAL_PI;

          ofs << " input \epsilon_{\alpha} for stopping condition (in degree) is" << 
                ALPHA_EPSILON << endl;

          ind ++;
          continue;
        }
      //--------------------------------------------
      // vibrate ratio from current point to end points of current edge
      if ( arg == "-v" )
        {
          VIBRATE_RATIO = atof( argv[ ind + 1 ] );
          ind ++;
          continue;
        }
      //--------------------------------------------
      // vibrate to move far from current edge
      if ( arg == "-n" )
        {
          VIBRATE_NEDGES = atoi( argv[ ind + 1 ] );
          ind ++;
          continue;
        }
      //--------------------------------------------

      // trang (2015/08/30): for computing the global solution
      if ( arg == "-g" )
      {
          globalComputing_ = true;

          ofs << " Find Global Path (not local): True" <<  endl;

          ind ++;
          continue;
      }

      // run LiuWong algorithm
      //g_renderer.LiuWong (arg, atof ((char *) (argv[ind + 1])),
      //                     string ((char *) (argv[ind + 2])), f_LiuWong);

      //------------ trang: main loop of gl renderer ---------------

      g_renderer.MainLoop (arg, atof ((char *) (argv[ind + 1])),
                           string ((char *) (argv[ind + 2])),
                           argc, argv, f_test, f_print, f_super_test,
                           f_white);

      
      

      ind += 2;
      ofs << " input max_iterators is: " << MAX_ITERATIONS << endl;
      ofs << " input Epsilon is: " << EPSILON << endl;

    }



  //Close file
  ofs.close();

  return 0;
}

/* main.cc - end of file */
