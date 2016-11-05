typedef struct obj_camera{
    double width, height;
} obj_camera;

typedef struct scene_object{
        unsigned char type;
        double color[4];
        double diffuse_color[4];
        double specular_color[4];
        double position[4];
        double normal[4];
        double radius;
        double reflectivity;
        double refractivity;
        double index_of_refraction;

} scene_object;

//Maybe put lights as scene objects as well?
typedef struct scene_light{
        unsigned char type;
        double color[4];
        double position[4];
        double direction[4];
        double theta;
        double ra2;
        double ra1;
        double ra0;
        double aa0;

} scene_light;
typedef struct pixels{
    unsigned char r,g,b;
} pixels;

