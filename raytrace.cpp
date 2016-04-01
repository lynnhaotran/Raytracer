//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include "matm.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cfloat>
using namespace std;

#define NO_INTERSECTION -1.0f
#define NOT_VISIBLE -2.0f
#define SHADOW -3.0f

int g_width;
int g_height;

struct Ray {
    vec4 origin;
    vec4 dir;
};

//DONE: add structs for spheres, lights and anything else you may need.

struct Sphere {
	vec4 center;
	vec3 scale;
	vec4 color;
	float k_a;
	float k_d;
	float k_s;
	float k_r;
	float n;
	mat4 inverse_scale;
};

struct Light {
	vec4 position;
	vec4 intensity;
};

vector<Sphere> spheres;
vector<Light> lights;
vector<vec4> g_colors;

//background r,g,b
vec4 bg_colors;
//scene's ambient intensity
vec4 ambience;
//output file name
string output;

bool inside = false;
//float t;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;


// -------------------------------------------------------------------
// Input file parsing

vec3 toVec3(vec4 in)
{
	return vec3(in[0], in[1], in[2]);
}

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    //DONE: add parsing of NEAR, LEFT, RIGHT, BOTTOM, TOP, SPHERE, LIGHT, BACK, AMBIENT, OUTPUT.

	const int num_labels = 11;
	const string labels[] = { "NEAR",		//0
							"LEFT",			//1
							"RIGHT",		//2
							"BOTTOM",		//3
							"TOP",			//4
							"RES",			//5
							"SPHERE",		//6
							"LIGHT",		//7
							"BACK",			//8
							"AMBIENT",		//9
							"OUTPUT"		//10
	};
	unsigned label_id = find(labels, labels + num_labels, vs[0]) - labels;
	switch (label_id) {
		case 0:		g_near = toFloat(vs[1]);	break;		//NEAR
		case 1:		g_left = toFloat(vs[1]);	break;		//LEFT
		case 2:		g_right = toFloat(vs[1]);	break;		//RIGHT
		case 3:		g_bottom = toFloat(vs[1]);	break;		//BOTTOM
		case 4:		g_top = toFloat(vs[1]);		break;		//TOP
		case 5:												//RES
			g_width = (int)toFloat(vs[1]);
			g_height = (int)toFloat(vs[2]);
			g_colors.resize(g_width * g_height);
			break;
		case 6:												//SPHERE
		{
			Sphere new_sphere;
			new_sphere.center = toVec4(vs[2], vs[3], vs[4]);
			new_sphere.scale = vec3(toFloat(vs[5]), toFloat(vs[6]), toFloat(vs[7]));
			new_sphere.color = toVec4(vs[8], vs[9], vs[10]);
			new_sphere.k_a = toFloat(vs[11]);
			new_sphere.k_d = toFloat(vs[12]);
			new_sphere.k_s = toFloat(vs[13]);
			new_sphere.k_r = toFloat(vs[14]);
			new_sphere.n = toFloat(vs[15]);

			//store the inverse scale transform for later
			mat4 scale_matrix = Scale(new_sphere.scale);
			mat4 inverseScale;
			InvertMatrix(scale_matrix, inverseScale);
			new_sphere.inverse_scale = inverseScale;

			spheres.push_back(new_sphere);
		}
			break;
		case 7:												//LIGHT
		{
			Light new_light;
			new_light.position = toVec4(vs[2], vs[3], vs[4]);
			new_light.intensity = toVec4(vs[5], vs[6], vs[7]);
			lights.push_back(new_light);
		}
			break;
		case 8:												//BACK		
			bg_colors = toVec4(vs[1], vs[2], vs[3]);
			break;		
		case 9:												//AMBIENT
			ambience = toVec4(vs[1], vs[2], vs[3]);
			break;
		case 10:											//OUTPUT
		{
			output = vs[1];
			break;
		}
	}

}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}


// -------------------------------------------------------------------
// Intersection routine

// TODO: add your ray-sphere intersection routine here.

float intersect(const Ray& ray, const Sphere& sphere, float threshold) {
		
	float t_1, t_2;

		//Transform ray to find intersection with unit sphere
		vec4 direction = sphere.inverse_scale * ray.dir;
		vec4 ray_to_origin = sphere.inverse_scale * (ray.origin - sphere.center);  
		
		//Components of determinant in line-sphere formula
		//referenced from Wikipedia - "Line-Sphere Intersection"
		float a = dot(direction, direction);
		float b = dot(direction, ray_to_origin);
		float c = dot(ray_to_origin, ray_to_origin) - (1.0f);

		//Evaluate determinant
		float determinant = (b*b) - (a*c);

		//No intersection
		if (determinant < 0)
			return NO_INTERSECTION;

		t_1 = ((-1 * b) - sqrtf(determinant)) / a;
		t_2 = ((-1 * b) + sqrtf(determinant)) / a;

		if (t_1 < threshold && t_2 < threshold) {
			//Handle shadows
			if (t_1 > 0.001f || t_2 > 0.001f)
				return SHADOW;
			return NOT_VISIBLE;
		}
		else if (t_1 >= threshold && t_2 >= threshold) {
			inside = false;
			if (t_2 > t_1)
				return t_1;
			else
				return t_2;
		}
		else if (threshold == 1.0f) {
			//Cutting inside the sphere
			inside = true;
			if (t_1 < threshold && t_2 >= threshold)
				return t_2;
			else if (t_2 < threshold && t_1 >= threshold)
				return t_1;
		}
		else
			return FLT_MAX;

}

vec4 c_local(const Ray& ray, const Sphere& sphere, const Light& lightsource, vec4 intersection, vec4 normal) {

	vec4 pixel_diffuse = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	vec4 pixel_specular = vec4(0.0f, 0.0f, 0.0f, 1.0f);

	//Diffuse Lighting
	vec4 light = normalize(lightsource.position - intersection);
	float angle = dot(normal, light);

	//Specular Lighting
	vec4 v = normalize((ray.origin - intersection));
	vec4 r = (2 * normal * angle) - light;
	float spec = pow(dot(r, v), sphere.n);


	if (angle < 0) {}
	else
		pixel_specular = (lightsource.intensity * spec * sphere.k_s);

	if (angle < 0) {}
	else
		pixel_diffuse = (lightsource.intensity * angle * sphere.k_d) * sphere.color;

	return pixel_diffuse + pixel_specular;
}


// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray, int depth, float threshold)
{
	vec4 pixel_color = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	vec4 pixel_local = vec4(0.0f, 0.0f, 0.0f, 1.0f);
	vec4 pixel_refl = vec4(0.0f, 0.0f, 0.0f, 1.0f);

	//bool has_intersection = false;
	float closest_t = FLT_MAX;
	Sphere intersected_sphere;
	bool inside_sphere = false;

	//Find the closest intersection between the ray from the eye and a sphere
	for (std::vector<Sphere>::iterator sphere = spheres.begin(); sphere != spheres.end(); ++sphere) {

		float t = intersect(ray, *sphere, threshold);

		if (t == NO_INTERSECTION)
			continue;
		else {
			//has_intersection = true;
			if (t < closest_t) {
				closest_t = t;
				intersected_sphere = *sphere;
				inside_sphere = inside;
			}
		}
	}

	vec4 intersection = ray.origin + closest_t*ray.dir;
	vec4 normal = (intersection - intersected_sphere.center) / (1.0f); 	//distance from center of transformed sphere to intersection at surface
	normal = intersected_sphere.inverse_scale*normal;					//untransform that distance
	normal = normalize(transpose(intersected_sphere.inverse_scale) * normal); 	//multiply by the inverse transpose to get the actual normal

	//invert the normal if the intersection is inside the sphere
	if (inside_sphere)
		normal = -normal;

	//If within the viewing plane, calculate shadow rays and local illumination
	if (closest_t != FLT_MAX && closest_t >= threshold) {
		for (std::vector<Light>::iterator light_it = lights.begin(); light_it != lights.end(); ++light_it) {
			Ray light_ray;
			light_ray.origin = intersection;
			light_ray.dir = light_it->position - intersection;

			//Flag for object obstructing light
			bool shadow = false;

			for (std::vector<Sphere>::iterator shadow_sphere = spheres.begin(); shadow_sphere != spheres.end(); ++shadow_sphere) {
				if (intersect(light_ray, *shadow_sphere, threshold) == SHADOW) {
					shadow = true;
					break;
				}
			}
			if (!shadow)
				pixel_local += c_local(ray, intersected_sphere, *light_it, intersection, normal);
		}
		pixel_local += ambience * intersected_sphere.k_a * intersected_sphere.color;

		//Calculate contribution from reflected rays
		vec4 calc_refl = vec4(0.0f, 0.0f, 0.0f, 1.0f);

		Ray reflection;
		reflection.origin = intersection;
		reflection.dir = -2 * dot(normal, normalize(ray.dir)) * normal + normalize(ray.dir);

		if (depth != 3)
			calc_refl = trace(reflection, depth + 1, 0.0f);

		if (calc_refl.x != bg_colors.x || calc_refl.y != bg_colors.y || calc_refl.z != bg_colors.z || calc_refl.w != bg_colors.w)
			pixel_refl = calc_refl * intersected_sphere.k_r;
	}

	if (closest_t == FLT_MAX || closest_t == NOT_VISIBLE)
		return bg_colors;
	else {
		pixel_color = pixel_local + pixel_refl;
		return pixel_color;
	}
}

//return the direction from the origin to pixel (ix, iy).
vec4 getDir(int ix, int iy)
{
    vec4 dir;

	//interpolate to find x and y values
	float x = g_left*(1 - ( (float) ix / g_width)) + g_right*( (float) ix / g_width);
	float y = g_bottom*(1 - ((float)iy / g_height)) + g_top* ((float)iy / g_height);
	float z = -1 * g_near;

	dir = vec4(x, y, z, 0.0f);

    return dir;
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    vec4 color = trace(ray, 0, 1.0f);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, const char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    // Clamp values if out of range.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
			for (int i = 0; i < 3; i++) {
				if (((float*)g_colors[y*g_width + x])[i] > 1.0f)
					((float*)g_colors[y*g_width + x])[i] = 1;
				if (((float*)g_colors[y*g_width + x])[i] < 0.0f)
					((float*)g_colors[y*g_width + x])[i] = 0;
				buf[y*g_width * 3 + x * 3 + i] = (unsigned char)(((float*)g_colors[y*g_width + x])[i] * 255.9f);
			}
    
    //DONE: change file name based on input file name.

    savePPM(g_width, g_height, output.c_str(), buf);
    delete[] buf;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
	return 0;
}

