#include <cstdio>
#include <unistd.h>
#include <stdlib.h>
#include <chrono>
#include <omp.h>
#include "BHTree.h"

#define TIME_STEPSIZE 1.0

void write_csv(TreeNode* root, Particle* p, FILE* fp, int time_idx, int num_particles){
    for (int i = 0; i < num_particles; i++){
        if(root->is_external(p[i])){
            printf("id: %d is external\n", i);
            continue;
        }
        fprintf(fp, "%d,%d,%lf,%lf,%lf,%lf,%lf,%lf\n",
                time_idx, i, time_idx*TIME_STEPSIZE, p[i].position.x, p[i].position.y,
                p[i].velocity.x, p[i].velocity.y, p[i].mass);
    }
}

void get_nodeinfo(int num_threads, int index, double& x, double& y, double& rad){
    // printf("get_nodeinfo, origin index %d\n", index);
    while(num_threads>1){
        num_threads/=4;
        int i = index/num_threads;
        index -= num_threads*i;
        rad/=2;
        x += (2*(i%2)-1)*rad;
        y += (1-2*(i/2))*rad;
    }
    // printf("rad %le, x %le, y %le\n", rad, x,  y);
}

TreeNode** create_trees(int num_threads, double global_rad){
    TreeNode** t = new TreeNode*[num_threads];
    for (int i = 0; i < num_threads; i++){
        double x = 0; 
        double y = 0;
        double rad = global_rad;
        get_nodeinfo(num_threads, i, x, y, rad);

        t[i] = new TreeNode(rad, x, y);
    }
    return t;
}

TreeNode* merge_tree(TreeNode** nodes, int num_nodes){
    while(num_nodes>1){
        TreeNode** upper_nodes = new TreeNode*[num_nodes/4];
        num_nodes = num_nodes/4;

        for(int i = 0; i < num_nodes; i++){
            double x = 0;
            double y = 0;
            double rad = 0;
            for(int j = 0; j < 4; j++){
                x += nodes[i*4+j]->center.x;
                y += nodes[i*4+j]->center.y;
                rad += nodes[i*4+j]->radius;
            }
            upper_nodes[i] = new TreeNode(rad/2, x/4, y/4);
            for(int j = 0; j < 4; j++){
                if(nodes[i*4+j]->num_particles>0){
                    upper_nodes[i]->children[j] = nodes[i*4+j];
                }
            }
        }
        delete[] nodes;
        nodes = upper_nodes;
    }
    TreeNode* root = nodes[0];
    delete[] nodes;
    return root;
}

int main(int argc, char *argv[]){
    using namespace std::chrono;
    typedef std::chrono::high_resolution_clock Clock;
    typedef std::chrono::duration<double> dsec;

    int opt = 0;

    char* filename = nullptr;
    char* outfile = nullptr;
    int num_particles;
    int num_iterations = 5;
    int num_of_threads = 1;
    double universe_radius = 0.0;

    do {
        opt = getopt(argc, argv, "f:o:i:n:");
        switch (opt) {
        case 'f':
            filename = optarg;
            break;

        case 'o':
            outfile = optarg;
            break;
        
        case 'n':
            num_of_threads = atoi(optarg);
            break;

        // case 'p':
        //     prob = atof(optarg);
        //     break;

        case 'i':
            num_iterations = atoi(optarg);
            break;

        case -1:
            break;

        default:
            break;
        }
    } while (opt != -1);

    printf("%s\n", filename);

    FILE* fp = fopen(filename,"r");
    fscanf(fp, "%d\n", &num_particles);
    printf("num of particles: %d\n", num_particles);

    fscanf(fp, "%lf\n", & universe_radius);
    printf("universe radius: %le\n", universe_radius);

    Particle* particles = (Particle*) calloc(num_particles, sizeof(Particle));

    char* line;
    size_t len = 0;
    ssize_t read;
    int count = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        sscanf(line, "%lf %lf %lf %lf %lf %*s", &particles[count].position.x, &particles[count].position.y, 
                &particles[count].velocity.x,  &particles[count].velocity.y, &particles[count].mass);
        count++;
    }
    // for (int i= 0; i < num_particles; i++) {
    //     printf("px: %le, py: %le, vx: %le, vy: %le, mass: %le\n", 
    //             particles[i].position.x, particles[i].position.y, particles[i].velocity.x, particles[i].velocity.y, particles[i].mass);
    // }
    fclose(fp);

    fp = fopen(outfile, "w");
    fprintf(fp, "time_idx,body_idx,t,x,y,vx,vy,m\n");
        
    TreeNode* root = new TreeNode(universe_radius, 0, 0);
    write_csv(root, particles, fp, 0, num_particles);

    omp_set_num_threads(num_of_threads);
    auto compute_start = Clock::now();
    double compute_time = 0;
    double tree_time = 0;
    double force_time = 0;

    for (int iter = 1; iter <= num_iterations; iter++) {
        auto tree_start = Clock::now();

        TreeNode** nodes = create_trees(num_of_threads, universe_radius);

        #pragma omp parallel for default(shared) schedule(dynamic) 
        for (int i = 0; i < num_particles; i++){
            particles[i].set_global_quad(num_of_threads, universe_radius);
        }

        #pragma omp parallel
        {
            int tid = omp_get_thread_num(); // Scratch vectors allocated at startup 
            for (int i = 0; i < num_particles; i++){
                if(tid == particles[i].global_quad){
                    nodes[tid]->add_particle(particles[i]);
                }
            }
        }

        int tid = omp_get_thread_num();
        if (tid == 0) {
            root = merge_tree(nodes, num_of_threads);
        }
        tree_time += duration_cast<dsec>(Clock::now() - tree_start).count();

        // root->print();
        auto force_start = Clock::now();

        #pragma omp parallel for default(shared) schedule(dynamic) 
        for (int i = 0; i < num_particles; i++){
            vector force = root->calculate_force(particles[i]);
            // force.print();
            vector acc = force.mul(1.0 / particles[i].mass);
            // acc.print();
            vector delta_position = particles[i].velocity.mul(TIME_STEPSIZE)
                .add(acc.mul(0.5*TIME_STEPSIZE*TIME_STEPSIZE));
            // delta_position.print();
            particles[i].position = particles[i].position.add(delta_position);
            particles[i].velocity = particles[i].velocity.add(acc.mul(TIME_STEPSIZE));
        }
        force_time += duration_cast<dsec>(Clock::now() - force_start).count();
        write_csv(root, particles, fp, iter, num_particles);
        delete root;
        
    }
    compute_time += duration_cast<dsec>(Clock::now() - compute_start).count();
    printf("Computation Time: %lf.\n", tree_time + force_time);
    printf("Tree Time: %lf.\n", tree_time);
    printf("Force Time: %lf.\n", force_time);

    fclose(fp);

    free(particles);

    return 0;
}