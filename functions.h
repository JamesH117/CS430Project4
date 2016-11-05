static inline double square(double v){
	return v*v;
}
static inline double power(double n, double r){
    double a = n;
    while(r!=1){
        a = a*n;
        r-=1;
    }
    return a;
}
static inline void normalize(double* v){
	double len = sqrt(square(v[0]) + square(v[1]) + square(v[2]));
	v[0] /= len;
	v[1] /= len;
	v[2] /= len;
}
static inline double vector_length(double* v){
    return sqrt(square(v[0])+square(v[1])+square(v[2]));

}
static inline void scale_vector(double scaler, double* vvector, double* answer){
    answer[0] = scaler * vvector[0];
    answer[1] = scaler * vvector[1];
    answer[2] = scaler * vvector[2];
}
static inline void add_vector(double* a, double* b, double* answer){
    answer[0] = a[0] + b[0];
    answer[1] = a[1] + b[1];
    answer[2] = a[2] + b[2];
}
static inline void sub_vector(double* a, double* b, double* answer){
    answer[0] = a[0] - b[0];
    answer[1] = a[1] - b[1];
    answer[2] = a[2] - b[2];
}
static inline int compare_objects(scene_object listo, scene_object pointero){
    if(listo.type != pointero.type) return -1;
    if(listo.diffuse_color[0] != pointero.diffuse_color[0] || listo.diffuse_color[1] != pointero.diffuse_color[1] || listo.diffuse_color[2] != pointero.diffuse_color[2]) return -1;
    if(listo.specular_color[0] != pointero.specular_color[0] || listo.specular_color[1] != pointero.specular_color[1] || listo.specular_color[2] != pointero.specular_color[2]) return -1;
    if(listo.position[0] != pointero.position[0] || listo.position[1] != pointero.position[1] || listo.position[2] != pointero.position[2]) return -1;
    if(listo.type == 'p'){
        if(listo.normal[0] != pointero.normal[0] || listo.normal[1] != pointero.normal[1] || listo.normal[2] != pointero.normal[2]) return -1;
    }
    if(listo.type == 's'){
        if(listo.radius != pointero.radius) return -1;
    }
    return 0;
}
static inline scene_object copy_object(scene_object coppy, scene_object original){
    //printf("orig type: %c\n", original.type);
    coppy.type = original.type;
    //printf("coppy type: %c\n", coppy.type);
    coppy.diffuse_color[0] = original.diffuse_color[0];
    coppy.diffuse_color[1] = original.diffuse_color[1];
    coppy.diffuse_color[2] = original.diffuse_color[2];

    coppy.specular_color[0] = original.specular_color[0];
    coppy.specular_color[1] = original.specular_color[1];
    coppy.specular_color[2] = original.specular_color[2];

    coppy.position[0] = original.position[0];
    coppy.position[1] = original.position[1];
    coppy.position[2] = original.position[1];

    if(original.type == 'p'){
        coppy.normal[0] = original.normal[0];
        coppy.normal[1] = original.normal[1];
        coppy.normal[2] = original.normal[2];
    }
    if(original.type == 's') coppy.radius = original.radius;
    return coppy;
}
static inline double clamp(double c){
    if(c>1)
        c = 1;
    if(c<0)
        c = 0;
        return c;
}
static inline double dot_product(double* a, double* b){
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}
static inline void vector_reflect(double* normal, double* light_vector, double* reflection){
    normalize(normal); //normalize normal just incase it is not already normalized
    normalize(light_vector);

    double a[3] = {0,0,0};
    /*
    scale_vector(2,normal,a);
    double num = dot_product(a,light_vector);
    scale_vector(num, normal, a);
    sub_vector(a,light_vector,reflection);
    */

    //scale_vector(1, light_vector, light_vector);

    double num = dot_product(normal, light_vector); //Dot product between normal and vector to reflect

    scale_vector(num,normal, a); //Scale the normal with that dot product
    scale_vector(2,a,a);        //Scale that new vector with 2
    sub_vector(light_vector,a,reflection);  //Subtract the new scaled vector from the vector to reflect
    //Reflection should now be the reflection of light_vector

    //free(a);//Free the memory used for this equation
}

//This Seems to be working with Spot Lights
static inline double f_ang(scene_light a, double* Rd_new, double PI){
    double* direction_from_light_to_object = malloc(sizeof(double)*3);
    //Rd_new is ray direction from object to light, want the inverse
    direction_from_light_to_object[0] = -1*Rd_new[0];
    direction_from_light_to_object[1] = -1*Rd_new[1];
    direction_from_light_to_object[2] = -1*Rd_new[2];

    if(a.type == 'l') return 1;

    double v0_v1 = ((a.direction[0]*direction_from_light_to_object[0]) +(a.direction[1]*direction_from_light_to_object[1]) + (a.direction[2]*direction_from_light_to_object[2]));

    if(v0_v1 < (cos(a.theta * PI/180))) return 0;

    free(direction_from_light_to_object);
    return pow(v0_v1, a.aa0);
}

static inline double f_rad(scene_light a, double* Ro_new){
    double* vector_to_light = malloc(sizeof(double)*3);
    double answer;

    //a.position - Ro_new = vector_to_light
    sub_vector(a.position, Ro_new, vector_to_light);
    //scale_vector(-1, vector_to_light, vector_to_light);
    double den = a.ra2*square(vector_length(vector_to_light)) + a.ra1*vector_length(vector_to_light) + a.ra0;
    double newden = 1*square(vector_length(vector_to_light));

    free(vector_to_light);

    if(den == 0.0) return 1/newden;
    return 1/den;

    free(vector_to_light);
}
static inline double diffuse_contribution(int index, double* obj_diff_color, scene_light light_obj, double N_dot_L){
    //normalize(closest_normal); //Normalize incase it is not already normalized
    //scale_vector(-1,vector_to_light,vector_to_light);
    if(N_dot_L <= 0) return 0;
    return obj_diff_color[index]*light_obj.color[index]*N_dot_L;
}
static inline double specular_contribution(int index, double* obj_spec_color, scene_light light_obj, double N_dot_L, double V_dot_R, double ns){
    //normalize(closest_normal);
    if(V_dot_R <= 0) return 0;
    if(N_dot_L <= 0) return 0;
    //double V_dot_R = dot_product(view_vector, reflection_of_light_vector);

    return obj_spec_color[index]*light_obj.color[index]*pow(V_dot_R,ns);
}
