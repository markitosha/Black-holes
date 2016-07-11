#include "Tracer.h"
#include <windows.h>
#include "glm/gtx/string_cast.hpp"
#include <omp.h>

using namespace glm;

double KP = 0;

SRay CTracer::MakeRay(glm::uvec2 pixelPos, double xPoint, double yPoint)
{
  SRay ray;
  ray.m_start = m_camera.m_pos;
  vec3 right(m_camera.m_right.x*((pixelPos.x+xPoint)/m_camera.m_resolution.x - 0.5),
	  m_camera.m_right.y*((pixelPos.x+xPoint)/m_camera.m_resolution.x - 0.5),
	  m_camera.m_right.z*((pixelPos.x+xPoint)/m_camera.m_resolution.x - 0.5));
   vec3 up(m_camera.m_up.x*((pixelPos.y+yPoint)/m_camera.m_resolution.y - 0.5),
	  m_camera.m_up.y*((pixelPos.y+yPoint)/m_camera.m_resolution.y - 0.5),
	  m_camera.m_up.z*((pixelPos.y+yPoint)/m_camera.m_resolution.y - 0.5));
  ray.m_dir = m_camera.m_forward + right + up;
  return ray;
}


vec3 InSphere(vec3 ray_pos, vec3 ray_dir, vec3 sph_pos, double rad, vec3 new_pos)
{
	
	double A = pow(ray_dir.x,2)+pow(ray_dir.y,2)+pow(ray_dir.z,2);
	double B = 2*(ray_pos.x - sph_pos.x)*ray_dir.x+2*(ray_pos.y - sph_pos.y)*ray_dir.y+2*(ray_pos.z - sph_pos.z)*ray_dir.z;
	double CK = pow(ray_pos.x - sph_pos.x, 2) + pow(ray_pos.y - sph_pos.y, 2) + pow(ray_pos.z - sph_pos.z, 2) - pow(rad,2);
	double D = pow(B,2) - 4*A*CK;
	if(D < 0)
		return vec3(0,0,0);
	else{
		double t1 = (-B - sqrt(D)) / (2*A);
		double t2 = (-B + sqrt(D)) / (2*A);
		double min_t = min(t1, t2);
		double max_t = max(t1, t2);
		double t = (min_t >= 0) ? min_t : max_t;
		if(t<0)
			return vec3(0,0,0);
		double x = ray_pos.x + ray_dir.x*t;
		double y = ray_pos.y + ray_dir.y*t;
		double z = ray_pos.z + ray_dir.z*t;
		if (length(vec3(x,y,z) - ray_pos) < length(new_pos - ray_pos))
			return vec3(x,y,z);
		else
			return vec3(0,0,0);
	}
}

vec3 InDisk(vec3 pos, vec3 ray_dir, vec3 sph_pos, double rad, vec3 new_pos, double phi)
{
	glm::mat3 matrix(1, 0, 0, 0, cos(phi), -sin(phi), 0, sin(phi), cos(phi));
	glm::vec3 n = vec3(0,0,1)*matrix;
	glm::vec3 dir = new_pos - pos;
	double A = n.x*dir.x + n.y*dir.y + n.z*dir.z;
	double B = n.x*pos.x + n.y*pos.y + n.z*pos.z - n.x*sph_pos.x - n.y*sph_pos.y - n.z*sph_pos.z;
	if(A == 0){
		if(B == 0)
			return pos;
		else
			return vec3(0,0,0);
	}
	double t = - B / A;
	if (t<0)
		return vec3(0,0,0);
	else{
		double x = pos.x + dir.x*t;
		double y = pos.y + dir.y*t;
		double z = pos.z + dir.z*t;	
		if (length(vec3(x,y,z) - pos) < length(new_pos - pos) && length(vec3(x,y,z) - sph_pos) < rad)
			return vec3(x,y,z);
		else
			return vec3(0,0,0);
	}
}

vec3 FonColor(vec3 dir, CImage *fon)
{
  dir = normalize(dir);
  double phi = atan2(dir.x, dir.y); 
  double eta = asin(dir.z);
  int H = 0;
  H = fon -> GetHeight();
  phi = (phi + PI)/PI*H;
  if(phi > 2*H - 0.5)
	phi = 0;
  eta = (eta + PI/2)/PI*H;
  if(eta > H - 0.5)
	eta = 0;
  int pitch = fon -> GetPitch();
  unsigned char* imageBuffer = (unsigned char*)fon -> GetBits();
  int r = ((int)(eta))*pitch + ((int)(phi))*3;
  int g = ((int)(eta))*pitch + ((int)(phi))*3 + 1;
  int b = ((int)(eta))*pitch + ((int)(phi))*3 + 2;
  return vec3((double)(imageBuffer[r]) / 255, (double)(imageBuffer[g]) / 255, (double)(imageBuffer[b])/ 255);	
}

vec4 TextureColor(vec3 w_pos, double phi, double rad, CImage *texture, double dist)
{
	glm::mat3 matrix(1, 0, 0, 0, cos(phi), -sin(phi), 0, sin(phi), cos(phi));
	vec3 n(0, 0, 1); 
	n = n * matrix;
	n = (n + vec3(w_pos.x, w_pos.y, w_pos.z - dist));
	n = n*inverse(matrix);
	vec2 t_pos = vec2(n.x + rad, n.y + rad);
	t_pos = vec2(t_pos.x/(2*rad)*texture->GetWidth(), t_pos.y/(2*rad)*texture->GetHeight());
	auto pData = (unsigned char *)texture->GetBits();
    int pitch = texture->GetPitch();
	int r = ((int)(t_pos.x))*pitch + ((int)(t_pos.y))*4 + 2;
	int g = ((int)(t_pos.x))*pitch + ((int)(t_pos.y))*4 + 1;
	int b = ((int)(t_pos.x))*pitch + ((int)(t_pos.y))*4;
	int alpha = ((int)(t_pos.x))*pitch + ((int)(t_pos.y))*4 + 3;
	if (pitch < 0 && (alpha > -pitch || r > -pitch || g > -pitch || b > -pitch)){
		do{
			r += pitch;
			g += pitch;
			b += pitch;
			alpha += pitch;
		}while(alpha > -pitch || r > -pitch || g > -pitch || b > -pitch);
	}
	if (pData[alpha] == 0)
		return vec4(-1,-1,-1, 0);
	else{
		return vec4((double)pData[r] / 255, (double)pData[g] / 255, (double)pData[b] / 255, (double)pData[alpha] / 255);
	}
}

vec3 SphColor(double a, vec3 color)
{
	vec3 alpha(a,a,a);
	return color*alpha;
}

vec3 DiskColor(double a, vec3 color0, vec3 color1)
{
	vec3 alpha(a,a,a);
	vec3 alpha1(1-a,1-a,1-a);
	return color0*alpha + color1*alpha1;
}

vec3 AlphaFon(vec3 dir, CImage *fon, double a, vec3 color)
{
	vec3 alpha(a,a,a);
	vec3 alpha1(1-a,1-a,1-a);
	return color*alpha + FonColor(dir, fon)*alpha1;
}

glm::vec3 CTracer::TraceRay(SRay& ray, CImage *fon, CImage *texture)
{
  vec3 color(0, 0, 1);
  CPhoton photon;
  vec3 l_speed(C*KP, C*KP, C*KP);
  photon.m_pos = ray.m_dir;
  photon.m_dst = m_pScene->m_hole.m_pos - photon.m_pos;
  double k, alph = 0;
  double dt = 1;
  vec3 delta_t(dt, dt, dt);
  vec3 delta_t2 = vec3(pow(delta_t.x, 2)/2, pow(delta_t.y, 2)/2, pow(delta_t.z, 2)/2);
  photon.m_dir = normalize(ray.m_dir)*l_speed*delta_t;
  vec3 old_dir, old_pos;
  vec3 t_sph(0,0,0), t_pln(0,0,0), add_pos(0,0,0);
  vec4 t_color(0,0,0,0);
  vec3 old_color(0,0,0);
  double sph_rad1 = m_pScene->m_hole.m_rad / 2, sph_rad2 = m_pScene->m_hole.m_rad / 10;
  int i = 0;
  vec3 t_sph1, t_sph2, sph_pos1(-m_camera.m_right.x/2, m_camera.m_up.y/4, m_camera.m_forward.z * RHO/2);
  vec3 sph_pos2(m_camera.m_right.x/7, m_camera.m_up.y/15, m_camera.m_forward.z);
  do{
	  i++;
	  old_dir = photon.m_dir;
	  old_pos = photon.m_pos;
	  k = KP*KP*G*m_pScene->m_hole.m_mas / pow(length(photon.m_dst), 3);
	  photon.m_acc = vec3(photon.m_dst.x*k, photon.m_dst.y*k, photon.m_dst.z*k);
	  add_pos = photon.m_dir*delta_t*l_speed + photon.m_acc*delta_t2;
	  photon.m_pos = photon.m_pos + add_pos;
	  photon.m_dir = normalize(photon.m_dir + photon.m_acc*delta_t)*l_speed*delta_t;
	  photon.m_dst = m_pScene->m_hole.m_pos - photon.m_pos;
	  t_sph = InSphere(old_pos, old_dir, m_pScene->m_hole.m_pos, m_pScene->m_hole.m_rad, photon.m_pos);
	  if(Sph != 0){
		  t_sph1 = InSphere(old_pos, old_dir, sph_pos1, sph_rad1, photon.m_pos);
		  if(t_sph1 != vec3(0,0,0))
			  return vec3(1.0/255,11.0/255,103.0/255);
		  t_sph2 = InSphere(old_pos, old_dir, sph_pos2, sph_rad2, photon.m_pos);
		  if(t_sph2 != vec3(0,0,0))
			  return vec3(103.0/255,1.0/255,13.0/255);
	  }
	  t_pln = InDisk(old_pos, old_dir, m_pScene->m_hole.m_pos, m_pScene->m_hole.m_radDisk, photon.m_pos, m_pScene->m_hole.m_phi);
	  if(t_sph != vec3(0,0,0) && t_pln != vec3(0,0,0)){
		  if(length(t_sph) < length(t_pln)){
			  if(AlphaD == 0)
				return vec3(0,0,0);
			  else
				return SphColor(alph, old_color);
		  }else{
			  t_color = TextureColor(t_pln, m_pScene->m_hole.m_phi, m_pScene->m_hole.m_radDisk, texture, m_pScene->m_hole.m_pos.z);
			  if(t_color != vec4(-1,-1,-1,0)){
				if(AlphaD == 0){
					return vec3(t_color.x, t_color.y, t_color.z);
				}else if(alph == 0){
					alph = t_color.w;
					old_color = vec3(t_color.x, t_color.y, t_color.z);
				}else{
					return DiskColor(alph, old_color, vec3(t_color.x, t_color.y, t_color.z));
				}	
			  }
		  }
	  }else if(t_sph != vec3(0,0,0)){
		if(AlphaD == 0)
		  return vec3(0,0,0);
		else
		  return SphColor(alph, old_color);
	  }else if(t_pln != vec3(0,0,0)){
		  t_color = TextureColor(t_pln, m_pScene->m_hole.m_phi, m_pScene->m_hole.m_radDisk, texture, m_pScene->m_hole.m_pos.z);
		  if(t_color != vec4(-1,-1,-1,0)){
				if(AlphaD == 0){
					return vec3(t_color.x, t_color.y, t_color.z);
				}else if(alph == 0){
					alph = t_color.w;
					old_color = vec3(t_color.x, t_color.y, t_color.z);
				}else{
					return DiskColor(alph, old_color, vec3(t_color.x, t_color.y, t_color.z));
				}
		  }
	  }
  }while(photon.m_pos.z < 2*m_pScene->m_hole.m_pos.z && i < 10000);
  if(AlphaD == 0)
	return FonColor(photon.m_dir, fon);
  else
	return AlphaFon(photon.m_dir, fon, alph, old_color);
}

vec3 CTracer::Antialiacing(int j, int i, CImage *fon, CImage *texture)
{
	vec3 color(0,0,0);
	SRay ray = MakeRay(uvec2(j, i), 0.5, 0.5);
	color += TraceRay(ray, fon, texture);
	ray = MakeRay(uvec2(j, i), 0.25, 0.25);
	color += TraceRay(ray, fon, texture);
	ray = MakeRay(uvec2(j, i), 0.75, 0.25);
	color += TraceRay(ray, fon, texture);
	ray = MakeRay(uvec2(j, i), 0.25, 0.75);
	color += TraceRay(ray, fon, texture);
	ray = MakeRay(uvec2(j, i), 0.75, 0.75);
	color += TraceRay(ray, fon, texture);
	color = vec3(color.x / 5, color.y / 5, color.z / 5);
	return color;
}

void CTracer::RenderImage(int xRes, int yRes)
{
	CImage* fon = LoadImageFromFile("data/stars1.jpg");
	CImage* texture = LoadImageFromFile("data/disk_32.png");

  // Rendering
  m_camera.m_resolution = uvec2(xRes, yRes);
  m_camera.m_pixels.resize(xRes * yRes);
  m_camera.m_pos = vec3(0,0,0);
  m_camera.m_right = vec3(xRes, 0, 0); //TODO
  m_camera.m_up = vec3(0, yRes, 0);
  m_camera.m_forward = vec3(0, 0, xRes / 2 / tan(Angle/2));
//  m_camera.m_viewAngle = vec2(Angle, yRes / 2 / m_camera.m_forward.z);
  KP = m_camera.m_forward.z*2/(RHO*m_pScene -> m_hole.m_radDisk - m_camera.m_forward.z);
  m_pScene -> m_hole.m_pos = vec3(0, 0, KP*RHO*m_pScene -> m_hole.m_radDisk);
  m_pScene -> m_hole.m_mas = KP*m_pScene -> m_hole.m_mas;
  m_pScene -> m_hole.m_rad = KP*m_pScene -> m_hole.m_rad;
  m_pScene -> m_hole.m_radDisk = KP*m_pScene -> m_hole.m_radDisk;

#pragma omp parallel for
  for(int i = 0; i < yRes; i++){
	  std::cout<<i<<std::endl;
    for(int j = 0; j < xRes; j++)
    {
	  if(Antial == 0){
		SRay ray = MakeRay(uvec2(j, i), 0.5, 0.5);
		m_camera.m_pixels[i * xRes + j] = TraceRay(ray, fon, texture);
	  }else{
		m_camera.m_pixels[i * xRes + j] = Antialiacing(j, i, fon, texture);
	  }
    }
  }
}

void CTracer::SaveImageToFile(std::string fileName)
{
  CImage image;

  int width = m_camera.m_resolution.x;
  int height = m_camera.m_resolution.y;

  image.Create(width, height, 24);
    
	int pitch = image.GetPitch();
	unsigned char* imageBuffer = (unsigned char*)image.GetBits();

	if (pitch < 0)
	{
		imageBuffer += pitch * (height - 1);
		pitch =- pitch;
	}

	int i, j;
	int imageDisplacement = 0;
	int textureDisplacement = 0;

	for (i = 0; i < height; i++)
	{
    for (j = 0; j < width; j++)
    {
      vec3 color = m_camera.m_pixels[textureDisplacement + j];

      imageBuffer[imageDisplacement + j * 3] = clamp(color.b, 0.0f, 1.0f) * 255.0f;
      imageBuffer[imageDisplacement + j * 3 + 1] = clamp(color.g, 0.0f, 1.0f) * 255.0f;
      imageBuffer[imageDisplacement + j * 3 + 2] = clamp(color.r, 0.0f, 1.0f) * 255.0f;
    }

		imageDisplacement += pitch;
		textureDisplacement += width;
	}

  image.Save(fileName.c_str());
	image.Destroy();
}

CImage* CTracer::LoadImageFromFile(std::string fileName)
{
  CImage* pImage = new CImage;

  if(SUCCEEDED(pImage->Load(fileName.c_str())))
    return pImage;
  else
  {
    delete pImage;
    return NULL;
  }
}