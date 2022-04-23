#include <cmath>

#define G 6.6743e-11
#define THETA 1.0

class vector{
    public:
        double x, y;

        vector(){
            x = 0;
            y = 0;
        }

        vector(double x_, double y_){
            x = x_; 
            y = y_;
        }

        vector add(vector v){
            vector result;
            result.x = x + v.x;
            result.y = y + v.y;
            return result;
        }

        vector sub(vector v){
            vector result;
            result.x = x - v.x;
            result.y = y - v.y;
            return result;
        }

        vector mul(double coeff){
            vector result;
            result.x = x * coeff;
            result.y = y * coeff;
            return result;
        }

        double abs(){
            return sqrt(x * x + y * y);
        }

        void print(){
            printf("x: %le, y: %le\n", x, y);
        }
};


class Particle{
    public:
        double mass;
        vector position;
        vector velocity;
        int global_quad;

    void set_global_quad(int num_threads, double radius){
        if(num_threads == 4){
            if (position.x <= 0 && position.y >=0){
                global_quad = 0;
            } else if (position.x >= 0 && position.y >= 0){
                global_quad = 1;
            } else if (position.x <= 0 && position.y <= 0){
               global_quad = 2;
            } else {
                global_quad = 3;
            }
        }
        
    }
};

enum class Quad{
    QuadTopLeft,
    QuadTopRight,
    QuadBottomLeft,
    QuadBottomRight
};




class TreeNode{
    public:
        double mass = 0.0;
        vector mass_center;
        vector center;
        double radius;
        TreeNode* children[4] = {};
        int num_particles = 0;
        Particle* particle;

        TreeNode(double radius, double center_x, double center_y);
        ~TreeNode();

        bool is_external(Particle& p);
        Quad get_quad(double px, double py);
        void build_subtree(Quad quad);
        void add_particle(Particle& p);
        vector calculate_force(Particle& p);

        void print();
};