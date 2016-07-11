#include "Tracer.h"
#include "stdio.h"


void main(int argc, char** argv)
{
  CTracer tracer;
  CScene scene;

  int xRes = 1024;  // Default resolution
  int yRes = 768;

  if(argc == 2) // There is input file in parameters
  {
    FILE* file = fopen(argv[1], "r");
    if(file)
    {
      int xResFromFile = 0;
      int yResFromFile = 0;
      if(fscanf(file, "%d %d\n", &xResFromFile, &yResFromFile) == 2)
      {
        xRes = xResFromFile;
        yRes = yResFromFile;
		if(xResFromFile < 512 || yResFromFile < 512){
			printf("Too small size of picture.\r\n");
			xRes = 1024;
			yRes = 768;
		}
      }
      else
        printf("Invalid config format! Using default parameters.\r\n");
	  double Mass = 0;
	  double coefRadDisk = 0;
	  double phi = 0;
	  CBlackHole hole;
	  if(fscanf(file, "%lf %lf %lf\n %d %d %d %d", &Mass, &coefRadDisk, &phi, &tracer.Antial, &tracer.AlphaD, &tracer.Sph) == 6)
	  {
		hole.m_mas = Mass;
		hole.m_rad = 2*G*Mass / pow(C, 2);
		hole.m_radDisk = coefRadDisk * hole.m_rad;
		/*if(xRes > yRes){
			hole.m_k = hole.m_radDisk/yRes*2;
		}else{
			hole.m_k = hole.m_radDisk/xRes*2;
		}*/
		hole.m_phi = phi*PI/180;
		scene.m_hole = hole;
	  }
	  else
        printf("Invalid config format! Using default parameters.\r\n");
      fclose(file);
    }
    else
      printf("Invalid config path! Using default parameters.\r\n");
  }
  else
    printf("No config! Using default parameters.\r\n");

  tracer.m_pScene = &scene;
  tracer.RenderImage(xRes, yRes);
  tracer.SaveImageToFile("Result.png");
  unsigned int end_time = clock();
  std::cout << end_time << std::endl;
  scanf("\n");
}