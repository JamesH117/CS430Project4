#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "structs.h"
#include "functions.h"

double* illuminate(int recursion_depth, double* color, int closest_object_index, double Ro[], double Rd[], double best_t);
double shoot_double(double* Ro, double* Rd);
int shoot_int(double* Ro, double* Rd);

int line = 1;
int num_camera;
int total_objects;
int list_i = 0;
int list_l = 0;
#define PI 3.141592653589793
obj_camera main_camera;
scene_object *obj_list;
scene_light *light_list;
pixels *pixel_buffer;

// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
  int c = fgetc(json);
#ifdef DEBUG
  printf("next_c: '%c'\n", c);
#endif
  if (c == '\n') {
    line += 1;
  }
  if (c == EOF) {
    fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
    exit(1);
  }
  return c;
}


// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
  int c = next_c(json);
  if (c == d) return;
  fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
  exit(1);
}


// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
  int c = next_c(json);
  while (isspace(c)) {
    c = next_c(json);
  }
  ungetc(c, json);
}


// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
  char buffer[129];
  int c = next_c(json);
  if (c != '"') {
    fprintf(stderr, "Error: Expected string on line %d.\n", line);
    exit(1);
  }
  c = next_c(json);
  int i = 0;
  while (c != '"') {
    if (i >= 128) {
      fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
      exit(1);
    }
    if (c == '\\') {
      fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
      exit(1);
    }
    if (c < 32 || c > 126) {
      fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
      exit(1);
    }
    buffer[i] = c;
    i += 1;
    c = next_c(json);
  }
  buffer[i] = 0;
  return strdup(buffer);
}

double next_number(FILE* json) {
  double value;
  fscanf(json, "%lf", &value);
  //printf("Value is: %lf\n", value);
  // Error check this..
  return value;
}

double* next_vector(FILE* json) {
  double* v = malloc(3*sizeof(double));
  expect_c(json, '[');
  skip_ws(json);
  v[0] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[1] = next_number(json);
  skip_ws(json);
  expect_c(json, ',');
  skip_ws(json);
  v[2] = next_number(json);
  skip_ws(json);
  expect_c(json, ']');
  return v;
}

//Error Test for if a value is not given
void read_scene(char* filename) {
    //Doesnt run with 4 lights without this print line....
    //HELP
    //printf("Not it is fixed, WHY.\n");
    //printf("Starting Read Scene Function.\n");
    int c;
    char current_object;

    FILE* json = fopen(filename, "r");

  if (json == NULL) {
    fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
    exit(1);
  }

  skip_ws(json);

  // Find the beginning of the list
  expect_c(json, '[');

  skip_ws(json);

  // Find the objects

  while (1) {
    c = fgetc(json);
    if (c == ']') {
       printf("%d", line);
      fprintf(stderr, "Error: This is the worst scene file EVER.\n");
      fclose(json);
      return;
    }
    if (c == '{') {
      skip_ws(json);

      // Parse the object
    char* key = next_string(json);
    if (strcmp(key, "type") != 0) {

        fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
        exit(1);
      }
    skip_ws(json);
    expect_c(json, ':');
    skip_ws(json);
    char* value = next_string(json);
    if (strcmp(value, "camera") == 0) {
        if(num_camera >1){
            fprintf(stderr, "Error: More than one camera object has been provided, will not render.");
            exit(1);
        }
        num_camera += 1;
        //Working with camera object
        //Put future values into camera object
        current_object = 'c';
        //printf("Working with a %camera\n", current_object);
      }
    else if (strcmp(value, "sphere") == 0) {
        //Working with Sphere object
        //Put future values into sphere object and then put into list
        total_objects += 1;
        current_object = 's';
        //printf("Current object: %c\n", current_object);
        obj_list[list_i].type = 's';
        //printf("Working with a %cphere\n", current_object);
      }
    else if (strcmp(value, "plane") == 0) {
        //Working with Plane object
        //Put future values into plane object and then put into list
        total_objects += 1;
        current_object = 'p';
        //printf("Current object: %c\n", current_object);
        obj_list[list_i].type = 'p';
        //printf("Working with a %clane\n", current_object);
      }
      else if (strcmp(value, "light") == 0) {
        total_objects += 1;
        current_object = 'l';
        //printf("Current object: %c\n", current_object);
        light_list[list_l].type = 'l';
        //printf("Working with a %clane\n", current_object);
      }
    else {
        fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
        exit(1);
      }
    skip_ws(json);
    int cam_vars = 0;
    while (1) {
	// , }

        c = next_c(json);
        if (c == '}') {
            // stop parsing this object
            break;
        }
        else if (c == ',') {
            // read another field
            skip_ws(json);
            char* key = next_string(json);
            skip_ws(json);
            expect_c(json, ':');
            skip_ws(json);
            if ((strcmp(key, "width") == 0) || (strcmp(key, "height") == 0) || (strcmp(key, "radius") == 0) || (strcmp(key, "theta") == 0)
                || (strcmp(key, "radial-a2") == 0) || (strcmp(key, "radial-a1") == 0) || (strcmp(key, "radial-a0") == 0) || (strcmp(key, "angular-a0") == 0)
                || (strcmp(key, "reflectivity") == 0) || (strcmp(key, "refractivity") == 0) || (strcmp(key, "ior") == 0)) {
                //Depending on which is key, put value into that object value
                double value = next_number(json);
                //printf("Single Value is: %lf.\n", value);

                if((strcmp(key, "width") == 0)){
                    main_camera.width = value;
                    cam_vars+=1;
                }

                if((strcmp(key, "height") == 0)){
                    main_camera.height = value;
                    cam_vars +=1;
                }

                if((strcmp(key, "radius") == 0)){
                    obj_list[list_i].radius = value;
                    //obj_vars +=1;
                }
                if((strcmp(key, "radial-a2") == 0)){
                    light_list[list_l].ra2 = value;
                    //obj_vars +=1;
                }
                if((strcmp(key, "radial-a1") == 0)){
                    light_list[list_l].ra1 = value;
                    //obj_vars +=1;
                }
                if((strcmp(key, "radial-a0") == 0)){
                    light_list[list_l].ra0 = value;
                    //obj_vars +=1;
                }
                if((strcmp(key, "angular-a0") == 0)){
                    light_list[list_l].aa0 = value;
                    //obj_vars +=1;
                }
                if((strcmp(key, "theta") == 0)){
                    light_list[list_l].theta = value;
                    //obj_vars +=1;
                }
                if((strcmp(key, "reflectivity") == 0)){
                    if(value < 0 || value > 1){
                        fprintf(stderr, "Reflectivity value is out of the range, 0.0-1.0\n");
                        exit(1);
                    }

                    obj_list[list_i].reflectivity = value;
                    //obj_vars +=1;
                    if(1 - obj_list[list_i].reflectivity - obj_list[list_i].refractivity < 0){
                        fprintf(stderr, "1 - reflectivity - refractivity is below 0.\n");
                        exit(1);
                    }
                }
                if((strcmp(key, "refractivity") == 0)){
                    if(value < 0 || value > 1){
                    fprintf(stderr, "Refractivity value is out of the range, 0.0-1.0\n");
                    exit(1);
                    }
                    obj_list[list_i].refractivity = value;
                    if(1 - obj_list[list_i].reflectivity - obj_list[list_i].refractivity < 0){
                        fprintf(stderr, "1 - reflectivity - refractivity is below 0.\n");
                        exit(1);
                    }
                    //obj_vars +=1;
                }
                if((strcmp(key, "ior") == 0)){
                    obj_list[list_i].index_of_refraction = value;
                    //obj_vars +=1;
                }

            }

            else if ((strcmp(key, "color") == 0) || (strcmp(key, "position") == 0) || (strcmp(key, "normal") == 0) || (strcmp(key, "direction") == 0) ||
                     (strcmp(key, "diffuse_color") == 0) || (strcmp(key, "specular_color") == 0)) {
                //Depending on which is key, put value into that object *value
                double* value = next_vector(json);
                //printf("Value is: %lf, %lf, %lf.\n", value[0], value[1], value[2]);

                if((strcmp(key, "color") == 0)){
                        //Add saving functionality for Diffuse and Specular Color
                    if(current_object != 'l'){
                        obj_list[list_i].color[0] = value[0];
                        obj_list[list_i].color[1] = value[1];
                        obj_list[list_i].color[2] = value[2];
                        obj_list[list_i].color[3] = '\0';
                    }
                    if(current_object == 'l'){
                        light_list[list_l].color[0] = value[0];
                        light_list[list_l].color[1] = value[1];
                        light_list[list_l].color[2] = value[2];
                        light_list[list_l].color[3] = '\0';
                    }

                }
                if((strcmp(key, "diffuse_color") == 0)){
                    obj_list[list_i].diffuse_color[0] = value[0];
                    obj_list[list_i].diffuse_color[1] = value[1];
                    obj_list[list_i].diffuse_color[2] = value[2];
                    obj_list[list_i].diffuse_color[3] = '\0';
                }
                if((strcmp(key, "specular_color") == 0)){
                    obj_list[list_i].specular_color[0] = value[0];
                    obj_list[list_i].specular_color[1] = value[1];
                    obj_list[list_i].specular_color[2] = value[2];
                    obj_list[list_i].specular_color[3] = '\0';
                }
                if((strcmp(key, "position") == 0)){
                        if(current_object != 'l'){
                            obj_list[list_i].position[0] = value[0];
                            obj_list[list_i].position[1] = value[1];
                            obj_list[list_i].position[2] = value[2];
                            obj_list[list_i].position[3] = '\0';
                            //obj_vars +=1;
                        }
                        if(current_object == 'l'){
                            light_list[list_l].position[0] = value[0];
                            light_list[list_l].position[1] = value[1];
                            light_list[list_l].position[2] = value[2];
                            light_list[list_l].position[3] = '\0';
                        }
                }

                if((strcmp(key, "normal") == 0)){
                    obj_list[list_i].normal[0] = value[0];
                    obj_list[list_i].normal[1] = value[1];
                    obj_list[list_i].normal[2] = value[2];
                    obj_list[list_i].normal[3] = '\0';
                    //obj_vars +=1;
                }

                if((strcmp(key, "direction") == 0)){
                    light_list[list_l].direction[0] = value[0];
                    light_list[list_l].direction[1] = value[1];
                    light_list[list_l].direction[2] = value[2];
                    light_list[list_l].direction[3] = '\0';
                    light_list[list_l].type = 's';
                }

            }
            else {
                fprintf(stderr, "Error: Unknown property, \"%s\", on line %d.\n",
                key, line);

            }
            skip_ws(json);
        }
        else {
            fprintf(stderr, "Error: Unexpected value on line %d\n", line);
            exit(1);
        }
      }
    skip_ws(json);
    c = next_c(json);
    if (c == ',') {
        //No Operation, another object is coming up
        //Iterate through list of objects by size of 1 object
        skip_ws(json);
        //printf("Before increment, list_i is: %d\n", list_i);
        if(current_object != 'c'){
            if(current_object != 'l'){
                //list_i += sizeof(scene_object);
                list_i++;
                //obj_vars = 0;
            }
            if(current_object == 'l'){
                //list_l += sizeof(scene_light);
                list_l++;
            }

            //printf("Incremented by %d\n", sizeof(scene_object));
        }
        if(current_object == 'c'){
            if(cam_vars != 2){
                fprintf(stderr, "ERROR: Your camera has either too little parameters or too many parameters.\n");
                exit(1);
            }
        cam_vars = 0;
        }
    }
    else if (c == ']') {
        //Iterated through all objects
        fclose(json);
        return;
    }
    else {
        fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
        exit(1);
      }
    }
  }
}

double plane_intersection(double* Ro, double* Rd, double* position, double* normal){
    normalize(normal);
    double a = normal[0];
    double b = normal[1];
    double c = normal[2];
    normalize(Rd);

    //D is the length of the shortest line segment from the origin to the plane
    //D is distance from origin/camera to the plane
    double x0 = position[0];
    double y0 = position[1];
    double z0 = position[2];

    double d = -(a*x0 + b*y0 + c*z0);

    double den = (a*Rd[0] + b*Rd[1] + c*Rd[2]);
    if(den == 0.0) return -1;
    double t = -(a*Ro[0] + b*Ro[1] + c*Ro[2] + d)/(a*Rd[0] + b*Rd[1] + c*Rd[2]);

    return t;
}

double sphere_intersection(double* Ro, double* Rd, double* position, double radius){
    //Center points of sphere
    double xc = position[0];
    double yc = position[1];
    double zc = position[2];
    normalize(Rd);

    double a = square(Rd[0])+square(Rd[1])+square(Rd[2]);
    double b = 2*(Rd[0]*(Ro[0]-xc) + Rd[1]*(Ro[1]-yc) + Rd[2]*(Ro[2]-zc));
    double c = (square((Ro[0]-xc)) + square((Ro[1]-yc)) + square((Ro[2]-zc)) - square(radius));

    double desc = (square(b) - 4*a*c);
    //If descriminate is negative, no intersection
    if(desc < 0.0) return INFINITY;

    double t0 = ((-b - sqrt(square(b) - 4.0*c*a))/(2.0*a));
    double t1 = ((-b + sqrt(square(b) - 4.0*c*a))/(2.0*a));

    //If descriminant is negative, imaginary number, dont return
    //If t0 is negative don't return it, return t1
    //Return t0 first
    //Return t1 next
    if(t0 > 0.0)return t0;
    if(t1 > 0.0) return t1;

    return -1;
}

void raycast(double num_width, double num_height){
    //printf("Starting Raycast function.\n");

    int x = 0;
    int y = 0;
    int i;

    //Where the camera is sitting
    double Ro[3] = {0.0, 0.0, 0.0};


    double width = main_camera.width;
    double height = main_camera.height;
    double N = num_width;
    double M = num_height;
    double pixel_width = width/N;
    double pixel_height = height/M;

    //Center of the view plane
    double center[2] = {0.0, 0.0};

    //c_xyz is center of the view plane

    //Distance from camera to view plane
    double p_z = 1;
    for(y=0; y<M; y+=1){

        double p_y = center[1] - height/2.0 + pixel_height*(y+0.5);

        for(x=0; x<N; x+=1){

            double p_x = center[0] - width/2.0 + pixel_width*(x+0.5);
            double Rd[3] = {p_x, p_y, p_z};
            //printf("Rd is: %lf %lf %lf\n", Rd[0], Rd[1], Rd[2]);
            normalize(Rd);

            double best_t = INFINITY;
            int closest_object_index = -1;
            for(i=0; i<=list_i; i++){
                //printf("i is: %d", i);
                double t = 0;
                if(obj_list[i].type == 's'){
                        t = sphere_intersection(Ro, Rd, obj_list[i].position, obj_list[i].radius);
                }
                if(obj_list[i].type == 'p'){
                        t = plane_intersection(Ro, Rd, obj_list[i].position, obj_list[i].normal);
                        //printf("Plane here\n");
                }
                //printf("Before copy the obj_list[j].type is: %c\n", obj_list[i].type);
                if(t > 0 && t < best_t){
                    best_t = t;
                    closest_object_index = i;

                    //printf("Closest_object.type: %c\n", closest_object.type);
                }
            }
            //printf("best_t is: %g.\n", best_t);
//After you loop through all the objects, I should have the closest object to the pixel I am looking at
            double color[3] = {0,0,0};

            color[0] = 0;
            color[1] = 0;
            color[2] = 0;

            double final_color[3] = {0,0,0};
            //If best_t is greater than 0 and not INFINITY, intersection is in view of camera
            if(best_t > 0 && best_t != INFINITY){
                //Painting Fucntion that illuminates the current pixel based upon the lights in the scene.
                memcpy(final_color,illuminate(3, color,closest_object_index,Ro,Rd, best_t), sizeof(double)*3);
                }
                int pos = (int)((M - y -1)*N +x);
                pixel_buffer[pos].r = (unsigned char)(clamp(final_color[0])*255);
                pixel_buffer[pos].g = (unsigned char)(clamp(final_color[1])*255);
                pixel_buffer[pos].b = (unsigned char)(clamp(final_color[2])*255);
        }
    }
}

double* illuminate(int recursion_depth, double* color, int closest_object_index, double Ro[], double Rd[], double best_t){
    //base case
    if(recursion_depth == 0 || best_t == INFINITY || closest_object_index < 0){
        double black[3] = {0,0,0};
        return black;
    }

    double Ro_new[3] = {0,0,0}; //Ro_new is vector from Light towards Object.  Current Pixel location of intersection.
    Ro_new[0] = (best_t * Rd[0]) + Ro[0];
    Ro_new[1] = (best_t * Rd[1]) + Ro[1];
    Ro_new[2] = (best_t * Rd[2]) + Ro[2];

    int j,k;
    for(j=0; j<=list_l; j++){ //Do Summation of Ambient + Diffuse + Emission Light

        //printf("Ro_new[0] is: %lf\n", Ro_new[0]);
        double Rd_new[3] = {0,0,0}; //Ray direction from object intersection position to Light position
        sub_vector(light_list[j].position, Ro_new, Rd_new);
        normalize(Rd_new);
        //printf("Ro_new is: %lf, %lf, %lf.\n", light_list[j].position[0], light_list[j].position[1], light_list[j].position[2]);
        double distance_light_to_pixel = sqrt(square(light_list[j].position[0]-Ro_new[0])+ square(light_list[j].position[1]-Ro_new[1])+ square(light_list[j].position[2]-Ro_new[2]));


        double best_t_shadow = INFINITY;
        int closest_shadow_index = -1;
        for(k=0; k<=list_i; k++){  //Iterate through objects to see what casts a shadow on pixel
            double t_shadow = 0;
            //if(compare_objects(obj_list[k], closest_object) == 0) continue;
            if(k == closest_object_index) continue;

            if(obj_list[k].type == 's'){
                    t_shadow = sphere_intersection(Ro_new, Rd_new, obj_list[k].position, obj_list[k].radius);
            }
            if(obj_list[k].type == 'p'){
                    t_shadow = plane_intersection(Ro_new, Rd_new, obj_list[k].position, obj_list[k].normal);
            }
            //if(best_t > square(vector_length(object_to_light))){
            //if(best_t > vector_length(object_to_light)){
            if(best_t > distance_light_to_pixel){ //Need to skip if object is past light
                //printf("best_t: %g.\n", best_t);
                //printf("distance from light to pixel: %g.\n", distance_light_to_pixel);
                continue;
            }
            if(t_shadow > 0.0 && t_shadow < best_t_shadow){
                best_t_shadow = t_shadow;
                closest_shadow_index = k;
                //printf("Closest_shadow_object.type: %c\n", closest_shadow_object.type);
            }
        }
        //if(closest_shadow_object.type == NULL){ //Color in the lighting for that pixel because there is no object casting a shadow on the closest object.
        //if(best_t_shadow == INFINITY){
        if(closest_shadow_index < 0){
            double closest_normal[3] = {0,0,0};
            if(obj_list[closest_object_index].type == 's') {
                    sub_vector(Ro_new, obj_list[closest_object_index].position, closest_normal); //Turn on for RINGS
                    //sub_vector(closest_object.position, Ro_new, closest_normal);

                    }
            if(obj_list[closest_object_index].type == 'p'){
                closest_normal[0] = obj_list[closest_object_index].normal[0];
                closest_normal[1] = obj_list[closest_object_index].normal[1];
                closest_normal[2] = obj_list[closest_object_index].normal[2];
            }

            normalize(closest_normal);
            //scale_vector(-1, closest_normal, closest_normal);

            double *vector_direction_to_light = Rd_new;  //vector from object intersection to Light position
            double reflection_of_light_vector[3] = {0,0,0};

            //vector is pointing from intersection point towards light, do I need to make it negative first and then dot product?
            vector_reflect(closest_normal, vector_direction_to_light, reflection_of_light_vector);
            //scale_vector(-1, reflection_of_light_vector, reflection_of_light_vector);
            //view_vector is the vector direction the camera sees an object at
            normalize(reflection_of_light_vector);

            double view_vector[3] = {0,0,0};
            view_vector[0] = 1*Rd[0];
            view_vector[1] = 1*Rd[1];
            view_vector[2] = 1*Rd[2];
            normalize(view_vector);

            double* current_diffuse = obj_list[closest_object_index].diffuse_color;
            double* current_specular = obj_list[closest_object_index].specular_color;
            double ns = 20; //property of diffuseness of object, will eventually be a property of objects

            normalize(closest_normal);
            normalize(vector_direction_to_light);

            double V_dot_R = dot_product(view_vector, reflection_of_light_vector);
            double N_dot_L = dot_product(closest_normal, vector_direction_to_light);


            //Setting values for reflection vector
            double rec_Ro[3] = {0,0,0};
            memcpy(rec_Ro,Ro_new, sizeof(double)*3);
            double rec_Rd[3] = {0,0,0};
            memcpy(rec_Rd,Rd_new, sizeof(double)*3);
            //vector_reflect(closest_normal, Rd_new, rec_Rd);
            //scale_vector(-1, rec_Rd, rec_Rd);
            double rec_best_t;
            rec_best_t = shoot_double(rec_Ro,rec_Rd);
            int rec_co_index = shoot_int(rec_Ro,rec_Rd);

            //printf("rec_best_t: %g, rec_co_index: %d.\n", rec_best_t, rec_co_index);


            color[0] += (1 - obj_list[closest_object_index].reflectivity - obj_list[closest_object_index].refractivity) *
                         f_rad(light_list[j], Ro_new) * f_ang(light_list[j],Rd_new,PI) *
                            (diffuse_contribution(0,current_diffuse,light_list[j],N_dot_L)
                             +
                            specular_contribution(0,current_specular,light_list[j],N_dot_L,V_dot_R,ns)
                             ) +
                        //Reflectivity Color
                        (obj_list[closest_object_index].reflectivity) * illuminate(recursion_depth-1,color, rec_co_index, rec_Ro, rec_Rd, rec_best_t)[0] +
                        //Refractivity Color
                        (obj_list[closest_object_index].refractivity) * 0
                        ;

            color[1] += (1 - obj_list[closest_object_index].reflectivity - obj_list[closest_object_index].refractivity) *
                            f_rad(light_list[j], Ro_new) * f_ang(light_list[j],Rd_new,PI) *
                            (diffuse_contribution(1,current_diffuse,light_list[j],N_dot_L)
                             +
                            specular_contribution(1,current_specular,light_list[j],N_dot_L,V_dot_R,ns)
                             ) +
                        //Reflectivity Color
                        (obj_list[closest_object_index].reflectivity) * illuminate(recursion_depth-1,color, rec_co_index, rec_Ro, rec_Rd, rec_best_t)[0] +
                        //Refractivity Color
                        (obj_list[closest_object_index].refractivity) * 0
                        ;

            color[2] += (1 - obj_list[closest_object_index].reflectivity - obj_list[closest_object_index].refractivity) *
                            f_rad(light_list[j], Ro_new) * f_ang(light_list[j],Rd_new,PI) *
                            (diffuse_contribution(2,current_diffuse,light_list[j],N_dot_L)
                             +
                            specular_contribution(2,current_specular,light_list[j],N_dot_L,V_dot_R,ns)
                             ) +
                        //Reflectivity Color
                        (obj_list[closest_object_index].reflectivity) * illuminate(recursion_depth-1,color, rec_co_index, rec_Ro, rec_Rd, rec_best_t)[0] +
                        //Refractivity Color
                        (obj_list[closest_object_index].refractivity) * 0
                        ;

        }//Not darkening shadows, just lighting up where there are no shadows
    }
    return color;
}

double shoot_double(double* Ro, double* Rd){
    normalize(Rd);
    double best_t = INFINITY;
    int closest_object_index = -1;
    int i;
    for(i=0; i<=list_i; i++){
        double t = 0;
        if(obj_list[i].type == 's'){
                t = sphere_intersection(Ro, Rd, obj_list[i].position, obj_list[i].radius);
        }
        if(obj_list[i].type == 'p'){
                t = plane_intersection(Ro, Rd, obj_list[i].position, obj_list[i].normal);
        }
        if(t > 0 && t < best_t){
            best_t = t;
            closest_object_index = i;
        }
    }
    return best_t;
}
int shoot_int(double* Ro, double* Rd){
    normalize(Rd);
    double best_t = INFINITY;
    int closest_object_index = -1;
    int i;
    for(i=0; i<=list_i; i++){
        double t = 0;
        if(obj_list[i].type == 's'){
                t = sphere_intersection(Ro, Rd, obj_list[i].position, obj_list[i].radius);
        }
        if(obj_list[i].type == 'p'){
                t = plane_intersection(Ro, Rd, obj_list[i].position, obj_list[i].normal);
        }
        if(t > 0 && t < best_t){
            best_t = t;
            closest_object_index = i;

        }
    }
    return closest_object_index;
}

int write(int w, int h, char *output_image){
    FILE *fp;
    char magic_number[2] = {'P', '6'};
    int width = w;
    int height = h;
    int j;

    fp = fopen(output_image, "wb");

    if(fp == 0){
        fprintf(stderr, "Error: Unable to create file for output image.\n");
        exit(1);
    }
    fwrite(magic_number, sizeof(magic_number), sizeof(magic_number)-1, fp);
    fprintf(fp,"\n%d %d", width, height);
    fprintf(fp,"\n%d", 255);

    fprintf(fp,"\n");
            for (j=0; j<width*height; j++){
                    fwrite(&pixel_buffer[j].r,1,1, fp);
                    fwrite(&pixel_buffer[j].g,1,1, fp);
                    fwrite(&pixel_buffer[j].b,1,1, fp);
            }
    return 0;
}

int main(int argc, char** argv) {
    if(argc != 5){
        fprintf(stderr, "Error: Not all arguments were provided or too many were given.\n");
        exit(1);
    }
    obj_list = malloc(sizeof(scene_object)*128);
    light_list = malloc(sizeof(scene_light)*128);
    //scene_object obj_list[128];
    //scene_light light_list[128];
    double N = (double)atoi(argv[1]);
    double M = (double)atoi(argv[2]);
    pixel_buffer = (pixels*)malloc(sizeof(pixels)*N*M);
    //memset(pixel_buffer, 255, 3*N*M);

    read_scene(argv[3]);
    raycast(N, M);
    write(N, M, argv[4]);

    //Free Memory I Allocate
    free(obj_list);
    free(light_list);
    free(pixel_buffer);

    //printf("Finished.\n");

    return 0;
}
